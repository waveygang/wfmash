# wfmash

_A DNA sequence read mapper based on mash distances and the wavefront alignment algorithm._

`wfmash` is a fork of [MashMap](https://github.com/marbl/MashMap) that implements base-level alignment using [WFA](https://github.com/Martinsos/WFA), via the [`wflign`](https://github.com/ekg/wflign) tiled wavefront global alignment algorithm.
It completes MashMap with a high-performance alignment module capable of computing base-level alignments for very large sequences.

## process

Each query sequence is broken into non-overlapping pieces defined by `-s[N], --segment-length=[N]`.
These segments are then mapped using MashMap's sliding minhash mapping algorithm and subsequent filtering steps.
To reduce memory, a temporary file is used to store initial mappings.
Each mapping location is then used as a target for alignment using WFA.

The resulting alignments always contain extended CIGARs in the `cg:Z:*` tag.
Approximate mapping (equivalent to `MashMap`) can be obtained with `-m, --approx-map`.

Mapping merging is disabled by default, as aligning merged approximate mappings with WFA under reasonable identity bounds can generate very long runtimes.
However, merging can be useful in some settings and is enabled with `-M, --merge-mappings`.

Sketching, mapping, and alignment are all run in parallel using a configurable number of threads.
The number of threads must be set manually, using `-t`, and defaults to 1.

## usage

`wfmash` has been developed to accelerate the alignment step in variation graph induction (the first step in the `seqwish` / `smoothxg` pipeline).
Suitable default settings are provided for this purpose.

Seven parameters shape the length, number, identity, and alignment divergence of the resulting mappings.

### mapping settings

The first three affect the structure of the mashmap2 mappings:

* `-s[N], --segment-length=[N]` is the length of the mapped and aligned segment (when `-N` is not set)
* `-N, --no-split` avoids splitting queries into segments, and instead maps them in their full length
* `-p[%], --map-pct-id=[%]` is the percentage identity minimum in the _mapping_ step
* `-n[N], --n-secondary=[N]` is the maximum number of mappings (and alignments) to report for each segment above `segment-length` (the number of mappings for sequences shorter than the segment length is defined by `-S[N], --n-short-secondary=[N]`, and defaults to 1)

### alignment settings

The last four essential parameters control the WFA alignment process and filter its output.

WF-min and WF-diff prune unlikely solutions from the set in consideration:

* `-l[N], --wf-min=[N]` the number of wavefronts is required to trigger reduction
* `-d[N], --wf-diff=[N]` prune wavefronts whose are more than WF-diff cells (on the diagonal) behind the max wavefront

The exact WFA may be computed if desired, which requires more time and memory but is equivalent to affine Needleman-Wunsch.
(Note that WFA already has adaptive features due to its formulation.)

* `-e, --exact-wfa` compute the exact WFA, don't use adaptive wavefront reduction

An alignment identity filter can be used to remove very low-quality alignments:

* `-a[%], --align-pct-id=[%]` defines the minimum percentage identity alignment to report from the _alignment_step

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

The build is orchestrated with cmake:

```
cmake -H. -Bbuild && cmake --build build -- -j 16
```

The `wfmash` binary will be in `build/bin`.
To clean up, just remove the build directory.

## <a name=“publications”></a>publications

- **Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, and Antonio Espinosa** ["Fast gap-affine pairwise alignment using the wavefront algorithm"](https://doi.org/10.1093/bioinformatics/btaa777) *Bioinformatics*, 2020.

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.

- **Martin Šošić and Mile Šikić** ["Edlib: a C/C ++ library for fast, exact sequence alignment using edit distance"](https://doi.org/10.1093/bioinformatics/btw753), *Bioinformatics*, 2017.
