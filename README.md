# wfmash

_A DNA sequence read mapper based on mash distances and the wavefront alignment algorithm._

`wfmash` is a fork of [MashMap](https://github.com/marbl/MashMap) that implements base-level alignment using [WFA](https://github.com/Martinsos/WFA), via the [`wflign`](https://github.com/ekg/wflign) tiled wavefront global alignment algorithm.
It completes MashMap with a high-performance alignment module capable of computing base-level alignments for very large sequences.

## process

Each query sequence is broken into non-overlapping pieces defined by `-s[N], --segment-length=[N]`.
These segments are then mapped using MashMap's sliding minhash mapping algorithm and subsequent filtering steps.
To reduce memory, a temporary file is used to store initial mappings.
Each mapping location is then used as a target for alignment using the wavefront inception algorithm in `wflign`.

The resulting alignments always contain extended CIGARs in the `cg:Z:*` tag.
Approximate mapping (equivalent to `MashMap2` with a block length filter) can be obtained with `-m, --approx-map`.

Sketching, mapping, and alignment are all run in parallel using a configurable number of threads.
The number of threads must be set manually, using `-t`, and defaults to 1.

## usage

`wfmash` has been developed to accelerate the alignment step in variation graph induction (the first step in the `seqwish` / `smoothxg` pipeline).
Suitable default settings are provided for this purpose.

Seven parameters shape the length, number, identity, and alignment divergence of the resulting mappings.

### mapping settings

These parameters affect the structure of the alignments:

* `-s[N], --segment-length=[N]` is the length of the mapping seed (when `-N` is not set), consecutive segment mappings are merged and the merged mapping is aligned
* `-l[N], --block-length-min=[N]` defines a minimum length filter on our merged mappings
* `-N, --no-split` avoids splitting queries into segments, and instead maps them in their full length
* `-p[%], --map-pct-id=[%]` is the percentage identity minimum in the _mapping_ step
* `-n[N], --n-secondary=[N]` is the maximum number of mappings (and alignments) to report for each segment above `--block-length-min` (the number of mappings for sequences shorter than the segment length is defined by `-S[N], --n-short-secondary=[N]`, and defaults to 1)

### all-to-all mapping

Together, these settings allow us to precisely define an alignment space to consider.
During all-to-all mapping, `-X` can additionally help us by removing self mappings from the reported set, and `-Y` extends this capability to prevent mapping between sequences with the same name prefix.


## examples

Map a set of query sequences against a reference genome:

```sh
wfmash reference.fa query.fa >aln.paf
```

Setting a longer segment length to reduce spurious alignment:

```sh
wfmash -s 50000 reference.fa query.fa >aln.paf
```

Self-mapping of sequences:

```sh
wfmash -X query.fa query.fa >aln.paf
```

## sequence indexing

`wfmash` provides a progress log that estimates time to completion.

This depends on determining the total query sequence length.
To prevent lags when starting a mapping process, users should apply `samtools index` to index query and target FASTA sequences.
The `.fai` indexes are then used to quickly compute the sum of query lengths.


## installation

### building from source

The build is orchestrated with `cmake`. At least GCC version 9.3.0 is required for compilation. You can check your 
version via:

``` bash
gcc --version
g++ --version
```

Clone the `wfmash` git repository and build with:

```
sudo apt install libjemalloc-dev # ubuntu
sudo pacman -S jemalloc          # archlinux

git clone --recursive https://github.com/ekg/wfmash.git
cd wfmash
cmake -H. -Bbuild && cmake --build build -- -j 16
```

The `wfmash` binary will be in `build/bin`. To clean up, just remove the build directory.

### Bioconda

`wfmash` recipes for Bioconda are available at https://anaconda.org/bioconda/wfmash.
To install the latest version using `Conda` execute:

``` bash
conda install -c bioconda wfmash
```

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

Add the following to your ~/.config/guix/channels.scm:

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


## <a name=“publications”></a>publications

- **Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, and Antonio Espinosa** ["Fast gap-affine pairwise alignment using the wavefront algorithm"](https://doi.org/10.1093/bioinformatics/btaa777) *Bioinformatics*, 2020.

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.

- **Martin Šošić and Mile Šikić** ["Edlib: a C/C ++ library for fast, exact sequence alignment using edit distance"](https://doi.org/10.1093/bioinformatics/btw753), *Bioinformatics*, 2017.
