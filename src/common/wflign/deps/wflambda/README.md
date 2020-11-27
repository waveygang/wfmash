# WFλ

## What is WFλ?

The wavefront alignment (WFA) algorithm is an exact gap-affine algorithm that takes advantage of homologous regions between sequences to accelerate the alignment process.

WFλ is an adaptation of the WFA implementation that supports arbitrary comparison functions between the pattern and target vector.
No sequence is provided to the algorithm.
Instead, it takes the alignment matrix dimension and a callback function that determines if the given matrix `(v, h)` index is a match or mismatch.

## Why does WFλ exist?

WFλ's abstract alignment process can be used to align very long sequences in low memory.
Rather than aligning the sequences directly, we align them in small parts, guided by the WFA algorithm's wavefront progression.
In this application, the callback function compares a (potentially-overlapping) small pair of subsequences of the query and target defined by the index pair `(v, h)`.
If alignment is successful within some bounds, then the callback returns true and WFλ considers the cell a match.
Otherwise, it is not a match.
Thus, the wavefront progression of WFλ bounds our evaluation of the full alignment matrix.
The total maximum memory is bounded by the cost of a single alignment problem and the cost of maintaining the WFA wavefronts.

## Building

Either make or cmake can be used to build wflambda.

```
cmake -H. -Bbuild && cmake --build build -- -j 16
```

The library can be included in another cmake project like this:

```
add_subdirectory(deps/wflambda EXCLUDE_FROM_ALL)
target_link_libraries(your_exe wflambda) # or target_link_libraries(your_exe wflambda)
```

## License

WFλ is distributed under MIT licence.

## Citation

**Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
