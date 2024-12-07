# wfmash

_**a pangenome-scale aligner**_

[![build and test](https://github.com/ekg/wfmash/actions/workflows/test_on_push.yml/badge.svg)](https://github.com/ekg/wfmash/actions/workflows/test_on_push.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/wfmash)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6949373.svg)](https://doi.org/10.5281/zenodo.6949373)

`wfmash` is an aligner for pangenomes that combines efficient homology mapping with base-level alignment. It uses MashMap 3.5 to find approximate mappings between sequences, then applies WFA (Wave Front Alignment) to obtain base-level alignments.

`wfmash` is designed to make whole genome alignment easy. On a modest compute node, whole genome alignments of gigabase-scale genomes should take minutes to hours, depending on sequence divergence. It can handle high sequence divergence, with average nucleotide identity between input sequences as low as 70%.

`wfmash` is the key algorithm in [`pggb`](https://github.com/pangenome/pggb) (the PanGenome Graph Builder), where it is applied to make an all-to-all alignment of input genomes that defines the base structure of the pangenome graph. It can scale to support the all-to-all alignment of hundreds of human genomes.

## Process

By default, `wfmash` breaks query sequences into non-overlapping segments (default: 1kb) and maps them using MashMap. Consecutive mappings separated by less than the chain gap (default: 2kb) are merged. Mappings are limited to 50kb in length by default, which allows efficient base-level alignment using WFA. This length limit is important because WFA's computational complexity is quadratic in the number of differences between sequences, not their percent divergence - meaning longer sequences with the same divergence percentage require dramatically more compute time.

For longer sequences, use `-m/--approx-mapping` to get approximate mappings only, which allows working with much larger segment and mapping lengths.

## usage

`wfmash` has been developed to accelerate the alignment step in variation graph induction (the first step in the `seqwish` / `smoothxg` pipeline).
Suitable default settings are provided for this purpose.

Seven parameters shape the length, number, identity, and alignment divergence of the resulting mappings.

### mapping settings

These parameters affect the structure of the mappings:

* `-s[N], --segment-length=[N]` is the length of the mapping seed (default: `1k`). The best pairs of consecutive segment mappings are merged where separated by less than `-c[N], --chain-gap=[N]` bases.
* `-l[N], --block-length-min=[N]` requires seed mappings in a merged mapping to sum to more than the given length (default 5kb).
* `-p[%], --map-pct-id=[%]` is the percentage identity minimum in the _mapping_ step
* `-n[N], --n-secondary=[N]` is the maximum number of mappings (and alignments) to report for each segment above `--block-length-min` (the number of mappings for sequences shorter than the segment length is defined by `-S[N], --n-short-secondary=[N]`, and defaults to 1)

By default, we obtain base-level alignments by applying a high-order version of WFA to the mappings.
Various settings affect the behavior of the pairwise alignment, but in general the alignment parameters are adjusted based on expected divergence between the mapped subsequences.
Specifying `-m, --approx-map` lets us stop before alignment and obtain the approximate mappings (akin to `minimap2` without `-c`).

### all-to-all mapping

Together, these settings allow us to precisely define an alignment space to consider.
During all-to-all mapping, `-X` can additionally help us by removing self mappings from the reported set, and `-Y` extends this capability to prevent mapping between sequences with the same name prefix.
When working with large sequence collections we frequently use [PanSN](https://github.com/pangenome/PanSN-spec) naming convention and `-Y'#'` to specify that we want to group mappings by prefix, which in this context means genome or haplotype groupings.


## input indexing

`wfmash` requires a FASTA index (`.fai`) for its reference ("target"), and benefits if both reference and query are indexed.
We can build these indexes on BGZIP-indexed files, which we recommend due to their significantly smaller size.
To index your sequences, we suggest something like:

```sh
bgzip -@ 16 ref.fa
samtools faidx ref.fa.gz
```

Here, we apply `bgzip` (from `htslib`) to build a line-indexable gzip file, and then use `samtools` to generate the FASTA index, which is held in 2 files:

```sh
$ ls -l ref.fa.gz*
ref.fa.gz
ref.fa.gz.gzi
ref.fa.gz.fai
```

## examples

Map a set of query sequences against a reference genome:

```sh
wfmash reference.fa query.fa >aln.paf
```

For mapping longer sequences without alignment, use -m with larger segment and max length values:

```sh
wfmash -m -s 50k -P 500k reference.fa query.fa >mappings.paf
```

Self-mapping of sequences:

```sh
wfmash -X query.fa query.fa >aln.paf
```

Or just

```sh
wfmash query.fa >aln.paf
```

## sequence indexing

`wfmash` provides a progress log that estimates time to completion.

This depends on determining the total query sequence length.
To prevent lags when starting a mapping process, users should apply `samtools index` to index query and target FASTA sequences.
The `.fai` indexes are then used to quickly compute the sum of query lengths.


## Installation

### Static binaries

We provide [static builds of wfmash releases](https://github.com/waveygang/wfmash/releases) targeted at the `x86-64-v3` instruction set.

### Bioconda

`wfmash` recipes for Bioconda are available at https://anaconda.org/bioconda/wfmash.
To install the latest version using `Conda` execute:

``` bash
conda install -c bioconda wfmash
```

## Building from Source

The build process for `wfmash` is managed using `CMake`, providing various options to customize the build.

### Prerequisites

Before building `wfmash`, you need the following dependencies installed on your system:

- GCC (version 9.3.0 or higher) or a recent version of Clang/LLVM
- CMake
- Zlib
- GSL
- HTSlib
- LibLZMA
- BZip2
- Threads
- OpenMP

On Ubuntu >20.04, these dependencies can be installed with the following command:

```sh
sudo apt install build-essential cmake zlib1g-dev libgsl-dev libhts-dev liblzma-dev libbz2-dev
```

### Clone the Repository

Clone the `wfmash` repository:

```sh
git clone https://github.com/waveygang/wfmash.git
cd wfmash
```

### Build Options

`wfmash` provides several CMake options to customize the build process:

- `BUILD_STATIC` (default: `OFF`): Build a static binary.
- `BUILD_DEPS` (default: `OFF`): Build external dependencies (htslib, gsl, libdeflate) from source. Use this if system libraries are not available or you want to use specific versions. HTSlib will be built without curl support, which removes a warning for static compilation related to `dlopen`.
- `BUILD_RETARGETABLE` (default: `OFF`): Build a retargetable binary. When this option is enabled, the binary will not include machine-specific optimizations (`-march=native`).

These can be mixed and matched.

### Building with System Libraries

To build `wfmash` using system libraries:

```sh
cmake -H. -Bbuild && cmake --build build -- -j 8
```

This command will configure and build `wfmash` in the `build` directory, using as many cores as you specify with the `-j` option.

### Building with External Dependencies

If you need to build with external dependencies, use the `BUILD_DEPS` option:

```sh
cmake -H. -Bbuild -DBUILD_DEPS=ON && cmake --build build -- -j 8
```

This will download and build the necessary external dependencies.

### Building a Static Binary

To build a static binary, use the `BUILD_STATIC` option:

```sh
cmake -H. -Bbuild -DBUILD_STATIC=ON && cmake --build build -- -j 16
```

### Building a Retargetable Binary

To build a retargetable binary, use the `BUILD_RETARGETABLE` option:

```sh
cmake -H. -Bbuild -DBUILD_RETARGETABLE=ON && cmake --build build -- -j 8
```

This will configure the build without `-march=native`, allowing the binary to be run on different types of machines.

### Installing

After building, you can install `wfmash` using:

```sh
cmake --install build
```

This will install the `wfmash` binary and any required libraries to the default installation directory (typically `/usr/local/bin` for binaries).

#### Tests

To build and run tests, change to build directory and

```sh
ctest .
```

#### Notes for distribution

If you need to avoid machine-specific optimizations, use the `CMAKE_BUILD_TYPE=Generic` build type:

```shell
cmake -H. -Bbuild -D CMAKE_BUILD_TYPE=Generic && cmake --build build -- -j 8
```

The resulting binary should be compatible with all x86 processors.

#### Notes for debugging/plotting

To enable the functionality of emitting wavefront plots (in PNG format), tables (in TSV format), and timing information, add the `-DWFA_PNG_TSV_TIMING=ON` option:

```shell
cmake -H. -Bbuild -D CMAKE_BUILD_TYPE=Release -DWFA_PNG_TSV_TIMING=ON && cmake --build build -- -j 3
```

Note that this may make the tool a little bit slower.

### nix

If you have `nix`, you can install directly from the repository via:

```shell
nix profile install github:waveygang/wfmash
```

For local development, from the wfmash repo directory:

```shell
nix build .#wfmash
```

And you can install into your profile from the source repo with:

```shell
nix profile install .#wfmash
```

### guix

If you have `guix`:

```shell
guix build -f guix.scm
```

To build guix in a development container, see the instructions in the header of [guix.scm](./guix.scm).

#### Docker and Singularity images with nix

Nix is also able to build an Docker image, which can then be loaded by Docker and converted to a Singularity image.

```
nix build .#dockerImage
docker load < result
singularity build wfmash.sif docker-daemon://wfmash-docker:latest
```

This can be run with Singularity like this:

```
singularity run wfmash.sif $ARGS
```

Where `$ARGS` are your typical command line arguments to `wfmash`.

### Guix

#### installing via the guix-genomics git repository

First, clone the guix-genomics repository:

``` bash
git clone https://github.com/ekg/guix-genomics
```

And install the `wfmash` package to your default GUIX environment:

``` bash
GUIX_PACKAGE_PATH=. guix package -i wfmash
```

Now `wfmash` is available as a global binary installation.

#### installing via the guix-genomics channel

Add the following to your `~/.config/guix/channels.scm`:

``` scm
  (cons*
(channel
  (name 'guix-genomics)
  (url "https://github.com/ekg/guix-genomics.git")
  (branch "master"))
%default-channels)
```

First, pull all the packages, then install `wfmash` to your default GUIX environment:

``` bash
guix pull
guix package -i wfmash
```

If you want to build an environment only consisting of the `wfmash` binary, you can do:

``` bash
guix environment --ad-hoc wfmash
```

For more details about how to handle Guix channels, go to https://git.genenetwork.org/guix-bioinformatics/guix-bioinformatics.git.

## running wfmash on a cluster

When aligning a large number of very large sequences, one wants to distribute the calculations across a whole cluster.
This can be achieved by dividing the approximate mappings `.paf` into chunks of similar difficult alignment problems using [split_approx_mappings_in_chunks.py](scripts/split_approx_mappings_in_chunks.py).

### example

1. We restrict `wfmash` to its approximate mapping phase.

```sh
wfmash -m reference.fa query.fa > approximate_mappings.paf
```

2. We use the Python script to split the approximate mappings into chunks. A good approximation of the number of chunks is the number of nodes on your cluster. In the following, we assume a cluster with 5 nodes.

```python
python3 split_approx_mappings_in_chunks.py approximate_mappings.paf 5
```
This gives us:

```sh
ls
approximate_mappings.paf.chunk_0.paf
approximate_mappings.paf.chunk_1.paf
approximate_mappings.paf.chunk_2.paf
approximate_mappings.paf.chunk_3.paf
approximate_mappings.paf.chunk_4.paf
```

3. Dependent on your cluster workload manager, create a command line to submit 5 jobs to your cluster.

One example without specifying a workflow manager:

```sh
wfmash -i approximate_mappings.paf.chunk_0.paf reference.fa query.fa > approximate_mappings.paf.chunk_0.paf.aln.paf
```

The resulting `.paf` can be directly plugged into [seqwish](https://github.com/ekg/seqwish).

```sh
# list all base-level alignment PAFs
PAFS=$(ls *.aln.paf | tr '\n' ',')
# trim of the last ','
PAFS=${PAFS::-1}
seqwish -s reference.fa -p $PAFS -g seqwish.gfa
```

### use [nf-core/pangenome](https://github.com/nf-core/pangenome)

If you have `Nextflow` and `Docker` or `Singularity` available on your cluster, the lines above can become a one-liner:

```sh
nextflow run nf-core/pangenome -r dev --input references.fa --wfmash_only --wfmash_chunks 5
```

This emits a `results/wfmash` folder which stores all the `wfmash` output.


## <a name=“publications”></a>publications

- **Santiago Marco-Sola, Jordan M. Eizenga, Andrea Guarracino, Benedict Paten, Erik Garrison, and Miquel Moreto**. ["Optimal gap-affine alignment in O (s) space"](https://doi.org/10.1093/bioinformatics/btad074). *Bioinformatics*, 2023.

- **Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, and Antonio Espinosa** ["Fast gap-affine pairwise alignment using the wavefront algorithm"](https://doi.org/10.1093/bioinformatics/btaa777) *Bioinformatics*, 2020.

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.
