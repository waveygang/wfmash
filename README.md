# wfmash

_A DNA sequence read mapper based on mash distances and the wavefront alignment algorithm._

`wfmash` is a fork of [MashMap](https://github.com/marbl/MashMap) that implements base-level alignment using the wavefront alignment algorithm [WFA](https://github.com/smarco/WFA).
It completes an alignment module in MashMap and extends it to enable multithreaded operation.
A single command-line interface simplfies usage.
The [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) output format is harmonized and made equivalent to that in [minimap2](https://github.com/lh3/minimap2), and has been validated as input to [seqwish](https://github.com/ekg/seqwish).
`wfmash` is also a fork of `edyeet`, which uses edlib to obtain an edit-distance based alignment, which is fast but may not be appropriate for many biological applications.

## process

Each query sequence is broken into pieces defined by `-s[N], --segment-length=[N]`.
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

Four parameters shape the length, number, identity, and alignment divergence of the resulting mappings:

* `-s[N], --segment-length=[N]` is the length of the mapped and aligned segment
* `-p[%], --map-pct-id=[%]` is the percentage identity minimum in the _mapping_ step
* `-n[N], --n-secondary=[N]` is the maximum number of mappings (and alignments) to report for each segment
* `-a[N], --align-max-score=[N]` defines the maximum score per bp in the WFA _alignment_ step (higher is more divergent)

Together, these settings allow us to precisely define an alignment space to consider.
During all-to-all mapping, `-X` can additionally help us by removing self mappings from the reported set.

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

## installation

Follow [`INSTALL.txt`](INSTALL.txt) to compile and install wfmash.

## <a name=“publications”></a>publications

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.
