# edyeet

`edyeet` is a fork of [MashMap](https://github.com/marbl/MashMap) that implements base-level alignment using [edlib](https://github.com/Martinsos/edlib).
It completes an alignment module in MashMap and extends it to enable multithreaded operation.
A single command-line interface simplfies usage.
The [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) output format is harmonized and made equivalent to that in [minimap2](https://github.com/lh3/minimap2), and has been validated as input to [seqwish](https://github.com/ekg/seqwish).

## usage

`edyeet`'s purpose is to accelerate the alignment step in variation graph induction.
Suitable default settings are provided for this purpose.
Mapping merging is disabled by default, as aligning merged approximate mappings with edlib under reasonable identity bounds can generate very long runtimes.

Four parameters shape the length, number, and identity of the resulting mappings:

* `-s[N], --segment-length=[N]` is the length of the mapped and aligned segment
* `-p[%], --map-pct-id=[%]` is the percentage identity minimum in the _mapping_ step
* `-n[N], --n-secondary=[N]` is the maximum number of mappings (and alignments) to report for each segment
* `-a[%], --align-pct-id=[%]` defines the minimum percentage identity allowed in the _alignment_ step

Together, these settings allow us to precisely define an alignment space to consider.
During all-to-all mapping, `-X` can additionally help us by removing self mappings from the reported set.
The number of threads must be set manually, using `-t`, and defaults to 1.

## examples

Map a set of query sequences against a reference genome:

```sh
edyeet reference.fa query.fa >aln.paf
```

Setting a longer segment length to reduce spurious alignment:

```sh
edyeet -s 50000 reference.fa query.fa >aln.paf
```

Self-mapping of sequences:

```sh
edyeet -X query.fa query.fa >aln.paf
```

## installation

Follow [`INSTALL.txt`](INSTALL.txt) to compile and install edyeet.

## <a name=“publications”></a>publications

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.
