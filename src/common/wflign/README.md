# wflign

_the we-flyin WFA-guided ultralong sequence aligner_

## what

`wflign` derives a global alignment for pairs of long sequences using very little memory and time.

Rather than directly aligning the sequences, `wflign` aligns them to each other in small pieces, typically _s_ = 1kbp.
It computes the global alignment with [WFλ](https://github.com/ekg/wflambda), a version of the affine-gapped reduced wavefront algorithm (WFA) based on an abstract match function λ.

WFλ is computed over a matrix representing the global alignment between the sequences at _s_-bp resolution.
Each cell in the matrix corresponds to the alignment of a specific pair of _s_-length fragment from the target and query.

To determine if the cell is a match (and guide WFλ), we compute the alignment using a fast pairwise alignment algorithm, indicating a match when we find an alignment within the given edit distance threshold.
Finally, we use the WFλ traceback to determine and emit the set of base-level alignments on the optimal path corresponding to our _s_-bp matrix cells.

## why

It isn't easy to align long sequences.
Generally, the time and memory required to compute a pairwise alignment increase quadratically with query and target length.
It can be impractical to directly align very long sequences, and we often require heuristic methods to determine regions of the input that are likely to contain a chain of matches that we can align with an expensive pairwise algorithm.

WFA suggests an efficient way to reduce the amount of computation required to obtain the optimal alignment between two sequences.
It reduces the cost to be quadratic in the score (penalty) of the optimal alignment.
Thus, it can be very cheap when sequences are highly similar and the optimal alignment score is low.
However, regions of the alignment of low quality can increase the memory and runtime costs of the algorithm.

To avoid such limitations, _wflign_ applies WFA at a high level, using it to guide the alignment process, and keeping the largest alignment problem size small.
Our implementation thus requires only the memory to align a single segment and 1/_s_ the memory of the full WFA.
This approach is more flexible than using a fixed-width band (it effectively has a variable bandwidth), requires no heuristic seeding step, and can benefit from parallel exploration of the wavefront.

## how

To build:

```
cmake -H. -Bbuild && cmake --build build -- -j 16
```

It can then be run with `wflign`.

## limitations

`wflign` assumes that the sequences are of suitable identity and oriented in the same direction.
The provided executable is intended for testing.

## todo

This is a work in progress.
To complete:

- Correctly compute traceback to bound PAF output (requires caching alignments).
- Explore recursive evaluation to find SVs (direct repeats and inversions) in unaligned regions of the query and target.
- Parallel wavefront expansion.
