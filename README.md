# wfmash

_**a pangenome-scale aligner**_

[![build and test](https://github.com/ekg/wfmash/actions/workflows/test_on_push.yml/badge.svg)](https://github.com/ekg/wfmash/actions/workflows/test_on_push.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/wfmash)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6949373.svg)](https://doi.org/10.5281/zenodo.6949373)

`wfmash` is an aligner for pangenomes that combines efficient homology mapping with base-level alignment. It uses MashMap 3.5 to find approximate mappings between sequences, then applies WFA (Wave Front Alignment) to obtain base-level alignments. MashMap 3.5 employs minmers, a generalization of minimizers that provides unbiased Jaccard similarity estimation for improved mapping accuracy.

`wfmash` is designed to make whole genome alignment easy. On a modest compute node, whole genome alignments of gigabase-scale genomes should take minutes to hours, depending on sequence divergence. It can handle high sequence divergence, with average nucleotide identity between input sequences as low as 70%. By default, `wfmash` automatically determines an appropriate identity threshold based on the ANI (Average Nucleotide Identity) distribution of your input sequences, using the median (50th percentile) for optimal balance between coverage and alignment quality.

`wfmash` is the key algorithm in [`pggb`](https://github.com/pangenome/pggb) (the PanGenome Graph Builder), where it is applied to make an all-to-all alignment of input genomes that defines the base structure of the pangenome graph. It can scale to support the all-to-all alignment of hundreds of human genomes.

## Algorithm Overview

`wfmash` performs alignment in several stages:

1. **Mapping**: Query sequences are broken into segments based on window size (default: 1kb) and mapped using MashMap with minmer sketches. Minmers are a generalization of minimizers that select multiple smallest k-mers per window, enabling unbiased Jaccard similarity estimation.

2. **Chaining**: Consecutive mappings separated by less than the chain gap (default: 2kb) are merged into longer approximate mappings.

3. **Filtering**: Various filters can be applied:
   - L1 filtering requires a minimum number of sketch hits (default: 3)
   - Plane-sweep filtering removes overlapping mappings
   - Hypergeometric filtering assesses mapping significance

4. **Scaffolding** (optional): For large-scale alignments, scaffolding identifies syntenic regions:
   - Chains are merged with larger gaps (default: 100kb) to form scaffolds
   - Only chains with sufficient segments (default: 10) are considered
   - Mappings are retained if they fall within a maximum distance (default: 100kb) from scaffold anchors
   - This helps focus alignment on truly homologous regions while filtering out spurious matches

5. **Alignment**: Filtered mappings are aligned at base-level using WFA. Mappings are limited to 50kb by default because WFA's complexity is quadratic in the number of differences.

For approximate mapping only, use `-m/--approx-mapping` to skip the alignment stage, which allows working with much larger segment and mapping lengths.

## Usage

```
wfmash [target.fa] [query.fa] {OPTIONS}
```

### Basic Examples

Map query sequences against a reference:
```sh
wfmash reference.fa query.fa >aln.paf
```

All-vs-all alignment (map a set of sequences to themselves):
```sh
wfmash sequences.fa >aln.paf
```

Output only approximate mappings without base-level alignment:
```sh
wfmash -m reference.fa query.fa >mappings.paf
```

For PanSN-formatted all-vs-all mapping, exclude mappings within the same genome:
```sh
wfmash -Y '#' pangenome.fa >aln.paf
```

### Parameter Groups

#### Minmer Sketching
* `-k[INT], --kmer-size=[INT]` - k-mer size (default: 15)
* `-s[INT], --sketch-size=[INT]` - number of minmers per window (default: auto-calculated)
* `-w[INT], --window-size=[INT]` - window size for minmer selection (default: 1k)

#### Mapping Parameters
* `-m, --approx-mapping` - output mappings only, no alignment
* `-p[FLOAT|aniXX[+/-N]], --map-pct-id=[FLOAT|aniXX[+/-N]]` - minimum identity percentage or ANI preset (default: ani50)
  * Fixed percentage: `-p 85` sets 85% identity threshold
  * ANI presets: `-p ani25` uses 25th percentile, `-p ani50` uses median (default)
  * Adjustments: `-p ani50-10` uses median minus 10%, `-p ani75+5` uses 75th percentile plus 5%
* `-n[INT], --mappings=[INT]` - number of mappings per segment (default: 1)
* `-l[INT], --block-length=[INT]` - minimum mapping block length (default: 0, no minimum)
* `-c[INT], --chain-jump=[INT]` - maximum gap to chain mappings (default: 2k)
* `-P[INT], --max-length=[INT]` - maximum mapping length for alignment (default: 50k)
* `-N, --no-split` - map each sequence as a single block

#### Filtering Options
* `-f, --no-filter` - disable all filtering
* `-M, --no-merge` - keep fragment mappings separate
* `-o, --one-to-one` - report only best mapping per query/target pair
* `-H[INT], --l1-hits=[INT]` - minimum sketch hits for L1 filter (default: 3)
* `-F[FLOAT], --filter-freq=[FLOAT]` - filter high-frequency minimizers (default: 0.0002)
* `--hg-filter=[n,Δ,conf]` - hypergeometric filter parameters (default: 1.0,0.0,99.9)

#### Scaffolding Parameters (for synteny filtering)
* `-S[INT], --scaffold-mass=[INT]` - minimum segments per scaffold (default: 10)
* `-D[INT], --scaffold-dist=[INT]` - maximum distance from scaffold anchors (default: 100k)
* `-j[INT], --scaffold-jump=[INT]` - maximum gap for scaffold chaining (default: 100k)
* `--scaffold-out=[FILE]` - output scaffold chains to FILE
* `--scaffold-overlap=[FLOAT]` - overlap threshold for scaffold chain filtering (default: 0.5)

#### Selection Filters
* `-X, --self-maps` - include self-mappings
* `-Y[C], --group-prefix=[C]` - exclude mappings within groups by prefix delimiter
* `-L, --lower-triangular` - only map seq_i to seq_j if i>j
* `-T[pfx], --target-prefix=[pfx]` - only map to targets with prefix
* `-Q[pfxs], --query-prefix=[pfxs]` - only map queries with prefix(es)

#### Alignment Parameters
* `-g[m,go1,ge1,go2,ge2], --wfa-params=[m,go1,ge1,go2,ge2]` - WFA gap costs (default: 5,8,2,24,1)
* `-E[INT], --target-padding=[INT]` - bases to extend target region
* `-U[INT], --query-padding=[INT]` - bases to extend query region

#### Output Options
* `-a, --sam` - output in SAM format (default: PAF)
* `-d, --md-tag` - include MD tag in output

#### System Parameters
* `-t[INT], --threads=[INT]` - number of threads (default: 1)
* `-I[FILE], --read-index=[FILE]` - load pre-built index from FILE
* `-W[FILE], --write-index=[FILE]` - save index to FILE
* `-b[SIZE], --batch=[SIZE]` - target index batch size (default: 4G)


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

## Advanced Examples

### Mapping longer sequences without alignment
For long sequences where you only need approximate mappings:
```sh
wfmash -m -w 50k -P 500k reference.fa query.fa >mappings.paf
```

### Standard alignment with default parameters
For typical whole-genome alignment (default: ani50, -S 10):
```sh
wfmash reference.fa query.fa >aln.paf
```

### Higher identity threshold
For very similar sequences only (e.g., 95% identity):
```sh
wfmash -p 95 reference.fa query.fa >aln.paf
```

### Using ANI presets
Automatically determine identity threshold from data:
```sh
# Use median ANI for balanced sensitivity/specificity
wfmash -p ani50 reference.fa query.fa >aln.paf

# Use 75th percentile minus 5% for higher sensitivity
wfmash -p ani75-5 reference.fa query.fa >aln.paf
```

### Multiple mappings per segment
To explore alternative alignments:
```sh
wfmash -n 3 reference.fa query.fa >aln.paf
```

### Pangenome all-vs-all with scaffolding
For large-scale pangenome construction with synteny filtering:
```sh
wfmash -Y '#' -S 10 -j 200k --scaffold-out scaffolds.paf pangenome.fa >aln.paf
```

### One-to-one mapping
To get only the best mapping between each query-target pair:
```sh
wfmash -o reference.fa query.fa >aln.paf
```

## Scaffolding for Large-Scale Alignments

Scaffolding is a powerful feature for filtering alignments to focus on syntenic regions. It's particularly useful for:
- Whole-genome alignments
- Pangenome construction  
- Reducing noise in highly repetitive sequences

The scaffolding algorithm:
1. Merges chains with large gaps (up to `-j/--scaffold-jump`, default 100kb)
2. Filters for chains with sufficient support (≥ `-S/--scaffold-mass` segments, default 5)
3. Keeps only mappings within `-D/--scaffold-dist` (default 100kb) of scaffold anchors

This effectively identifies and preserves large-scale syntenic blocks while filtering out spurious matches.

## Sequence Indexing

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
git clone --recursive https://github.com/waveygang/wfmash.git
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

### Building with Vendored htslib

If your system doesn't have htslib installed (libhts-dev package), you can use the `VENDOR_HTSLIB` option to download and build htslib automatically:

```sh
cmake -H. -Bbuild -DVENDOR_HTSLIB=ON && cmake --build build -- -j 8
```

This option:
- Downloads htslib 1.20 from the official GitHub releases
- Builds it with minimal dependencies (no libcurl, S3, or GCS support)
- Links wfmash against the vendored htslib library
- Works with both shared and static builds (`-DBUILD_STATIC=ON`)

This is particularly useful for building on systems where htslib is not available through the package manager or when you need a specific version of htslib.

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

### Install with Guix in a local git repository

```shell
guix build -f guix.scm
```

To build guix in a development container, see the instructions in the header of [guix.scm](./guix.scm).
Note that our guix setup allows for static builds and specifying the target CPU architecture(!)

For example `--tune=native` builds for skylake on my laptop:

```shell
guix build -L . wfmash-gcc-static-git --without-tests=wfmash-gcc-static-git --tune=native
  guix build: tuning wfmash-gcc-static-git@0.21-HEAD.024dcaa for CPU skylake
  guix build: tuning gsl@2.8 for CPU skylake
...
```

To build for x86-64-v4 use `--tune=x86-64-v4`. A complete list can be found [here](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html).

See also the instructions in [guix.scm](./guix.scm).

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

wfmash is part of Guix:

```
guix package -A wfmash
wfmash  0.21.0  out     gnu/packages/bioinformatics.scm:24769:2
```

To compile a more recent version use the instructions in [guix.scm](./guix.scm).

#### installing via the guix-genomics git repository

*Note: this section is out of date.*

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


## Performance Optimization

### Default Settings Rationale

The default settings (`-p ani50 -S 10`) have been optimized based on extensive benchmarking:

- **ani50**: Uses the median ANI of input sequences, providing optimal balance between sensitivity and specificity
  - Reduces spurious inter-chromosomal mappings by up to 95%
  - Improves runtime by 2-8x compared to more permissive settings
  - Focuses on high-confidence regions (>99% identity in many cases)

- **-S 10**: Requires 10 segments for scaffold chains, filtering noise while preserving synteny
  - Reduces off-diagonal mappings by ~40%
  - Minimal impact on coverage for most use cases
  - Particularly effective for eukaryotic genomes with repetitive content

### Performance Tips

#### For Maximum Speed (Gene-space focus)
```sh
# Use higher identity threshold and more stringent scaffolding
wfmash -p ani75 -S 15 reference.fa query.fa >aln.paf
```
- 10x faster on divergent genomes
- Focuses on highly conserved regions
- Ideal for: gene annotation transfer, synteny analysis

#### For Maximum Sensitivity (Including repeats)
```sh
# Use lower identity threshold and relaxed scaffolding
wfmash -p ani25-5 -S 5 reference.fa query.fa >aln.paf
```
- Captures centromeric and highly repetitive regions
- ~90% query coverage on divergent genomes
- Ideal for: repeat analysis, centromere studies, comprehensive coverage

#### For Closely Related Genomes (>99% ANI)
```sh
# Default settings work well, or use fixed threshold
wfmash -p 99 reference.fa query.fa >aln.paf
```

### Memory Management

If you encounter memory issues:

1. **Reduce batch size**: Use `-b 500m` or `-b 100m` instead of default
2. **Reduce threads**: Memory usage scales with thread count
3. **Use scaffolding**: `-S 10` reduces the number of alignments computed

### Benchmarking Results

Based on real-world testing:

| Setting | Arabidopsis Runtime | Yeast Runtime | Memory Usage | Noise Level |
|---------|-------------------|---------------|--------------|-------------|
| ani25-5 -S 5 | 5:15 | 2:35 | 2.1 GB | High (17% inter-chr) |
| ani50 -S 10 | 0:38 (8.4x faster) | 1:19 (2x faster) | 1.5 GB | Low (1.7% inter-chr) |

## <a name="publications"></a>publications

- **Santiago Marco-Sola, Jordan M. Eizenga, Andrea Guarracino, Benedict Paten, Erik Garrison, and Miquel Moreto**. ["Optimal gap-affine alignment in O (s) space"](https://doi.org/10.1093/bioinformatics/btad074). *Bioinformatics*, 2023.

- **Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, and Antonio Espinosa** ["Fast gap-affine pairwise alignment using the wavefront algorithm"](https://doi.org/10.1093/bioinformatics/btaa777) *Bioinformatics*, 2020.

- **Bryce Kille, Erik Garrison, Todd J. Treangen, Adam M. Phillippy**. ["Minmers are a generalization of minimizers that enable unbiased local Jaccard estimation"](https://doi.org/10.1093/bioinformatics/btad512). *Bioinformatics*, 2023.

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.
