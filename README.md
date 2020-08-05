# edyeet

------------

`edyeet` is a fork of [MashMap](https://github.com/marbl/MashMap) that implements base-level alignment using [edlib](https://github.com/Martinsos/edlib).
It completes an alignment module in MashMap and extends it to enable multithreaded operation.
A single command-line interface simplfies usage.
The [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md) output format is harmonized and made equivalent to that in [minimap2](https://github.com/lh3/minimap2), and has been validated as input to [seqwish](https://github.com/ekg/seqwish).

## Installation

Follow [`INSTALL.txt`](INSTALL.txt) to compile and install edyeet.

## Usage

* Map set of query sequences against a reference genome:
  ```sh
  edyeet reference.fna query.fa >aln.paf
  ```

## <a name=“publications”></a>Publications

- **Chirag Jain, Sergey Koren, Alexander Dilthey, Adam M. Phillippy, and Srinivas Aluru**. ["A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps"](https://doi.org/10.1093/bioinformatics/bty597). *Bioinformatics (ECCB issue)*, 2018.

- **Chirag Jain, Alexander Dilthey, Sergey Koren, Srinivas Aluru, and Adam M. Phillippy**. ["A fast approximate algorithm for mapping long reads to large reference databases."](https://link.springer.com/chapter/10.1007/978-3-319-56970-3_5) In *International Conference on Research in Computational Molecular Biology*, Springer, Cham, 2017.
