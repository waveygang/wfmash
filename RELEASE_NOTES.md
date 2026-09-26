Release notes mirrored from the upstream GitHub releases of
[waveygang/wfmash](https://github.com/waveygang/wfmash/releases); see also
[contributors](https://github.com/waveygang/wfmash/graphs/contributors) and
[commits](https://github.com/waveygang/wfmash/commits/master).

## ChangeLog ongoing

## ChangeLog v0.24.3 (pre)

* Added a `--legacy` switch that restores the pre-0.15 (mashmap v3.1.1-era) parameter defaults in one flag (k-mer 19, 5 kb window, 25 kb block floor, 20 kb chain gap, one mapping per segment, fixed 90% identity, no 50 kb mapping cap, scaffold filtering off, 0.001% high-frequency k-mer filter, and an alignment identity cutoff of 0.8 x `-p`); explicit options always win. Documented in the README with the default-vs-legacy table and measured equivalence (8-strain yeast all-vs-all: 620/621 old pairs and 99.8% of old chain bases recovered, 2,243 vs 1,622 mappings - close but not byte-identical by design)
* Tracked the release version in a top-level `VERSION` file (0.24.2) and use it for display
* Fixed the failing GitHub Actions build (`build and test`): both jobs died at `Install wgatools` before compiling anything. The unpinned `cargo install --git wgatools` re-resolves dependencies and now pulls hts-sys 2.2.x, whose vendored htslib renamed `bam1_core_t.isize` to `isize_` (plus usize/u64 signature changes), breaking compilation of rust-htslib 0.44.1 (E0609 no field `isize`, E0308 mismatched types); wgatools' own Cargo.lock pins the compatible hts-sys 2.1.1, so `cargo install --locked` builds it cleanly (verified for wgatools v1.1.0 and pafcheck v0.1.0). Also pinned the runner images (`ubuntu-22.04` as upstream does, `macos-15`) because `ubuntu-latest`/`macos-latest` rolled to toolchains that cannot build this dependency at all. Bumped `actions/checkout` from v2 to v5 to clear the Node.js 20 deprecation warnings on current runners

## ChangeLog vlocal (agc-short-names worktree)

* Built the `agc-short-names` branch (AGC-archive input, PRs #400/#401) with libagc v2.1 for the workshop AGC experiment; see `cwl/agc.md` in the workshop tree for results and the two small build fixes needed (kernel headers for raduls, missing `#include <limits>` in wflign rkmh.cpp)

## ChangeLog v0.24.2 (20260130)

Buildable source tarball: [wfmash-v0.24.2.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.24.2/wfmash-v0.24.2.tar.gz)

### Bug Fixes

* **Fix input-seeds parameter conflict**: Changed `-E` to `-K` to avoid conflict with existing parameters (#392)
* **Fix stdout mapping truncation**: Use persistent output streams to prevent truncated output (#377)
* **Fix OOM handler deadlock**: Handler now exits immediately instead of potentially hanging (#387)
* **Fix ANI logging**: Correctly displays ani50-2 offset

### New Features

* **External seed input support**: New `-K/--input-seeds` parameter for FastGA integration with two-level plane sweep filtering (#391)
* **Scaffold overlap threshold**: New `--scaffold-overlap` parameter (default 0.5) (#385)

### Improvements

* **Optimized defaults**: Changed default ANI parameter to `ani50-2` based on performance analysis (#388)
* **Better OOM handling**: Aborts cleanly with helpful error message (#384)
* **Comprehensive dependency vendoring**: Added `VENDOR_EVERYTHING` build option with static linking support (#375)

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.24.1...v0.24.2

## ChangeLog v0.24.1 (20250729)

Buildable source tarball: [wfmash-v0.24.1.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.24.1/wfmash-v0.24.1.tar.gz)

These changes specifically work on ensuring that the updates to the mapping system are memory safe. There were some cases of uninitialized variables being used in computation that may lead to instability in the previously released version.

## ChangeLog v0.24.0 (20250728)

Buildable source tarball: [wfmash-v0.24.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.24.0/wfmash-v0.24.0.tar.gz)

## wfmash v0.24.0 Release

This release brings significant memory optimizations, improved mapping scaffolding capabilities, and enhanced ANI-based identity estimation.

### Major Improvements

#### Memory Optimization
* **Drastically reduced memory usage** during mapping phase (~66% reduction)
* Optimized alignment phase with on-demand record loading
* Clean memory separation between mapping and alignment phases
* Optimized ANI sketching phase memory consumption
* New compact mapping structures for better memory efficiency

#### Mapping Scaffolding
* **New 2D distance graph scaffolding algorithm** for improved syntenic block detection
* Enhanced scaffold filtering with plane sweep optimization
* Configurable minimum scaffold length (default: 5kb)
* Support for scaffold mapping output via new options
* Better handling of boundary mappings for improved leniency

#### ANI-based Identity Estimation
* **New ANI preset system** (ani25, ani25-5, etc.) for automatic identity threshold selection
* Automatic identity estimation with `-p auto`
* Streaming MinHash implementation for efficient ANI computation
* Parallel ANI estimation with TaskFlow
* Per-group identity calculations with better CPU utilization

### New Features

#### Build System
* Added `VENDOR_HTSLIB` CMake option for building without system htslib
* Updated WFA2-lib submodule integration
* Improved build optimization flags

#### Command Line Interface
* Redesigned CLI parameters for better usability
* Changed sketch parameter to `-s` (was `-S`)
* Changed window-size parameter to `-w` (was `segLength`)
* Updated default overlap threshold from 1.0 to 0.95
* Minimum L1 hits now defaults to 3 (configurable with `-H`)
* Map sparsification parameter for controlling mapping density

#### Performance
* ~25% speedup for small genomes through optimized reverseComplement function
* Per-group mutexes for better parallel scaling
* Thread-local reader functions for improved I/O performance
* Progress reporting for all pipeline phases

### Bug Fixes
* Fixed critical bug in MinHash sketch computation for groups
* Resolved `stoi` conversion errors with invalid records
* Fixed type conversion overflow in parameter handling
* Corrected mapping merge logic for query and reference spans
* Fixed boundary mapping criteria for better edge case handling

### Technical Details
* Maintained chain identity information for `ch:Z:` tag
* Helper function to merge adjacent CIGAR operations
* Improved sequence loading patterns for ANI estimation
* Better progress reporting throughout all phases

### Contributors
Thanks to all contributors who made this release possible, with special mentions to those who worked on memory optimization, scaffolding improvements, and the ANI estimation system.

### What's Changed
* Save mappings at the boundaries by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/350
* Created a submodule for WFA2-lib to reflect credits by @pjotrp in https://github.com/waveygang/wfmash/pull/351
* Improved ReverseComplement function to use a lookup table by @pjotrp in https://github.com/waveygang/wfmash/pull/352
* Update WFA2-lib subproject to latest commit 49c255df by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/353
* Sparsify again by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/360
* Fix type conversion in handy_parameter function to use int64_t for co… by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/363
* Lower density minmers by @ekg in https://github.com/waveygang/wfmash/pull/362
* fix: apply query padding consistently in total work calculation by @ekg in https://github.com/waveygang/wfmash/pull/364
* Merge adjacent CIGAR operations by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/365
* Resolved skipping invalid record: stoi by @unavailable-2374 in https://github.com/waveygang/wfmash/pull/367
* Low memory mapper by @ekg in https://github.com/waveygang/wfmash/pull/369
* improve readme to match current version by @ekg in https://github.com/waveygang/wfmash/pull/370
* -p from ANI estimate by @ekg in https://github.com/waveygang/wfmash/pull/371
* Add VENDOR_HTSLIB CMake option for building without system htslib by @ekg in https://github.com/waveygang/wfmash/pull/372
* squeeze down memory by @ekg in https://github.com/waveygang/wfmash/pull/373

### New Contributors
* @unavailable-2374 made their first contribution in https://github.com/waveygang/wfmash/pull/367

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.23.0...v0.24.0

## ChangeLog v0.23.0 (20250414)

Buildable source tarball: [wfmash-v0.23.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.23.0/wfmash-v0.23.0.tar.gz)

### What's Changed
* fix: Update script to generate git version with additional source path by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/328
* Mapping memory cleanup by @ekg in https://github.com/waveygang/wfmash/pull/329
* Optimize mapping process for memory efficiency by @ekg in https://github.com/waveygang/wfmash/pull/330
* Alignment badness by @ekg in https://github.com/waveygang/wfmash/pull/331
* Fix index management when target and queries are in different files by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/332
* Thread safe faidx by @ekg in https://github.com/waveygang/wfmash/pull/333
* Head/tail patching + new alignment penalties by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/335
* Indicators progress logging by @ekg in https://github.com/waveygang/wfmash/pull/334
* Add option to disable alignment patching at chain boundaries by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/338
* query padding by @ekg in https://github.com/waveygang/wfmash/pull/339
* feat: Update WFA scoring parameters to [5,8,2,24,1] by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/340
* Apply target padding everywhere, while query padding only at the ends by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/341
* Add `--progress-bar` option by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/342
* fix: change batch size type from uint64_t to int64_t for consistency by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/343
* guix: fix profiler to run with two tests by @pjotrp in https://github.com/waveygang/wfmash/pull/336
* remove unused code by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/344
* refactor: update align progress meter to use `shared_ptr` for consistency by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/345
* refactor `mergeMappingsInRange` by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/346
* Fix progress again by @ekg in https://github.com/waveygang/wfmash/pull/348

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.22.0...v0.23.0

## ChangeLog v0.22.0 (20250305)

Buildable source tarball: [wfmash-v0.22.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.22.0/wfmash-v0.22.0.tar.gz)

### Overview

wfmash v0.22.0 represents a significant evolution in our pangenome-scale alignment approach, featuring fundamental algorithmic changes to mapping and alignment processes, a cleaner command line interface. It presents a major reworking of key parts of the mapping and alignment pipeline to improve reliability, sensitivity, and accuracy. This version introduces mutual-best-buddy based mapping chaining, smaller segment sizes for greater SV breakpoint sensitivity, scaffold-mapping based filtering to detect large homology regions, long mapping splitting, (which allows) direct biWFA alignment, improved memory management with batch indexing, and a completely rewritten TaskFlow-based execution model that delivers substantial performance improvements.

### Alignment Engine Improvements

#### Direct biWFA Integration

The most substantial change in this release is the transition from WFlign to biWFA as the default alignment algorithm. Previously, wfmash used a complex hierarchical approach with WFlign and several intermediate alignment steps, which sometimes led to inconsistent results and performance issues. The new direct biWFA implementation provides several benefits:

* Simpler, more reliable alignment process with fewer intermediate steps
* More consistent alignment results across diverse sequence types
* Improved handling of complex structural variations
* Better performance through vectorized code for lower divergence sequences

The alignment engine now uses a direct approach to apply biWFA to mappings found in the initial MashMap phase, resulting in cleaner alignments while reducing computational overhead.

#### Target Padding

We've implemented target padding around mapping boundaries to improve alignment quality. When requested through the new `-E/--target-padding` parameter, wfmash extends the reference sequence by a specified amount on both sides of the mapped region before alignment. This helps capture sequence context that might be missed with exact boundaries and ensures more complete alignments at sequence edges.

To maintain valid alignments, we employ coordinate swizzling during the alignment process, ensuring that any indels resulting from target padding are always placed at alignment boundaries, even though we're using global alignment with biWFA.

### Mapping and Filtration Improvements

#### Fundamental Changes to Mapping and Chaining

One of the most significant architectural changes in this version is the complete rewrite of the mapping and chaining logic. These changes fundamentally alter how wfmash detects and represents sequence homology:

* Segment length defaults to 1kb (versus previous 5kb), allowing detection of much finer-grained homology
* Chain gap parameter now defaults to 2kb (versus previous 30kb), providing more precise control over what constitutes a chain
* Chain selection now uses a more sophisticated distance-based metric that considers both reference and query coordinates
* Chains are now properly tracked with unique identifiers, positions, and lengths throughout the process
* Chain information is preserved in output formats with `chain:i:id.pos.length` tags
* Improved merging logic respects maximum mapping length while maintaining chain integrity
* Better handling of divergent regions with intelligent chain splitting and statistics computation
* Support for merging chains on either forward or reverse strand with proper coordinate handling

These changes collectively result in much more accurate mapping chains that better represent the underlying biological reality, especially for sequences with complex evolutionary histories.

#### Scaffold-Based Mapping

A major advancement in this release is the introduction of scaffold-based mapping, which substantially improves our ability to detect and represent structural variations. The scaffolding process works by:

1. Creating "super chains" from maximally merged mappings with aggressive gap parameters
2. Using a rotated coordinate system to efficiently filter mappings based on their relationship to these scaffolds
3. Employing a plane sweep algorithm with interval trees to identify mappings that fall within scaffold envelopes

This approach preserves mappings that contribute to larger structural patterns while filtering out spurious alignments. Users can control this process with the new `-S/--scaffolding` parameter which accepts gap size, minimum length, and maximum deviation values.

#### Improved Chaining

The chaining logic has been completely rewritten to better handle complex genomic arrangements. The new approach:

* Focuses on finding optimal chain pairs based on precise distance metrics
* Uses shorter segment lengths (default 1kb versus previous 5kb) to capture finer-grained homology
* Implements more intelligent and flexible chain gap parameters (default 2kb vs. previous 30kb)
* Respects maximum mapping length constraints when merging chains

These changes allow wfmash to detect smaller structural variants directly in the mapping phase, rather than relying on the more complex and sometimes error-prone WFlign approach used in previous versions.

### Memory Optimization and Scalability

#### Batch Indexing

This release introduces batch indexing for reference sequences, a critical feature for working with very large genomes on memory-constrained systems. With the new `-b/--batch` parameter, wfmash:

* Partitions reference sequences into batches based on specified size
* Builds and processes indices for each batch independently
* Combines mapping results across all batches for consistent output

This enables processing of reference collections that would otherwise exceed available memory, such as large mammalian pangenome projects spanning terabases of sequence.

#### Parallel Index Building

Index construction is now fully parallelized, with significant performance improvements:

* Parallel k-mer frequency counting across all sequences
* Thread-local processing during index construction
* Improved synchronization between indexing threads
* Better memory utilization through optimized data structures

The result is much faster index building, particularly for large reference collections, with reduced peak memory usage.

### Performance Enhancements

#### TaskFlow-Based Execution Model

The execution engine has been completely rewritten using the TaskFlow library, replacing the previous atomic queue-based system. This new model provides:

* More efficient task scheduling with explicit dependencies
* Better load balancing across all available threads
* Improved pipeline parallelism for multi-stage processing
* Reduced thread contention and synchronization overhead

The TaskFlow implementation manages the entire workflow from reading input files to writing output, with appropriate parallelization at each stage.

#### Enhanced FASTA I/O Performance

FASTA input/output operations now benefit from several optimizations:

* Integrated thread pooling for FASTA reading with optimized BGZF queue sizes
* Zero-copy sequence view processing with string_view
* Batch reading of input mapping files
* Efficient memory handling with custom allocators

These improvements significantly reduce I/O bottlenecks, especially for highly compressed reference files.

### Command Line Interface Improvements

The command-line interface has been thoroughly reorganized for better usability:

* More logical grouping of related parameters
* Short options for commonly used parameters (e.g., `-g` for alignment scoring)
* Clearer parameter names that better reflect their functionality
* More descriptive help text with improved formatting
* Simplified parameter handling with sensible defaults

New parameters include:

* `-S/--scaffolding` for controlling scaffold-based mapping
* `-E/--target-padding` for reference sequence padding
* `-b/--batch` for controlling batch size in indexing
* `-g/--wfa-params` for alignment scoring configuration

Default values have been carefully tuned based on extensive testing across diverse genome types, providing good out-of-the-box performance for most use cases.

### Output and Reporting Enhancements

#### Chain Information in PAF Output

The PAF output format now includes a `chain:i` field that exposes detailed information about mapping chains:

* Chain ID to identify mappings that belong together
* Position within the chain (1-based)
* Total length of the chain

This makes it easier to track related mappings and understand the structure of complex alignments, particularly when processing outputs with downstream tools.

#### Enhanced Progress Reporting

Progress reporting has been significantly improved:

* More accurate time estimates during long-running operations
* Detailed statistics about reference and query datasets
* Information about sequence groups and average sizes
* Clear reporting of filtering parameters and their effects
* Better error messages with more context and recovery options

These improvements make wfmash more informative during execution and help troubleshoot potential issues.

### Other Notable Changes

* HTSlib thread pooling integration for better performance with compressed files
* Improved error handling and validation for all input parameters
* More robust handling of sequence naming and ID management
* Enhanced detection and processing of overlapping mappings
* Optimization of minimum hits calculation for better sensitivity
* Support for chain field exposed in PAF output for downstream processing
* Improved hypergeometric filtering with configurable parameters

### Conclusion

wfmash v0.22.0 represents a substantial step forward in pangenome alignment capability, with fundamental improvements to core algorithms, significantly enhanced performance, and better memory efficiency. These changes enable more accurate alignment of complex genomic regions while making the tool more accessible for large-scale projects on diverse computing environments.

### What's Changed
* Fix query end position by @ekg in https://github.com/waveygang/wfmash/pull/271
* Smooth chain by @ekg in https://github.com/waveygang/wfmash/pull/272
* smooth the introduction of max mapping length parameter by @ASLeonard in https://github.com/waveygang/wfmash/pull/269
* fix: invalid paf produced for some patch alignments by @kdm9 in https://github.com/waveygang/wfmash/pull/274
* Fix uint64_t underflow by @bkille in https://github.com/waveygang/wfmash/pull/280
* Map chunk query by @ekg in https://github.com/waveygang/wfmash/pull/277
* Super basic hypergeometric filter by @ekg in https://github.com/waveygang/wfmash/pull/284
* Subindexes in one file and remove frequent kmer filtering by @ekg in https://github.com/waveygang/wfmash/pull/282
* Parallel filter by @ekg in https://github.com/waveygang/wfmash/pull/285
* Log indexing by @ekg in https://github.com/waveygang/wfmash/pull/287
* biWFA it by @ekg in https://github.com/waveygang/wfmash/pull/288
* Flagtastrophe by @ekg in https://github.com/waveygang/wfmash/pull/290
* Tweak hits by @ekg in https://github.com/waveygang/wfmash/pull/292
* Freq filter yes by @ekg in https://github.com/waveygang/wfmash/pull/293
* Fix minmer filt by @ekg in https://github.com/waveygang/wfmash/pull/294
* Let me overlap by @ekg in https://github.com/waveygang/wfmash/pull/295
* Full precision cli by @ekg in https://github.com/waveygang/wfmash/pull/296
* Parallelize k-mer frequency counting and index building by @ekg in https://github.com/waveygang/wfmash/pull/297
* Logging target subset count and average size by @ekg in https://github.com/waveygang/wfmash/pull/298
* Reenable writing index only by @ASLeonard in https://github.com/waveygang/wfmash/pull/300
* Prefilter mappings to save memory in batched mapping by @ekg in https://github.com/waveygang/wfmash/pull/301
* Guix static build for wfmash. Fixing tests and working on guix instructions/scripts by @pjotrp in https://github.com/waveygang/wfmash/pull/302
* wfmash now builds with clang by @pjotrp in https://github.com/waveygang/wfmash/pull/304
* New guix build targets and shells. Fixes libasan and adds profiling. by @pjotrp in https://github.com/waveygang/wfmash/pull/305
* map: update mapping selection logic based on merge and split parameters by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/306
* apply mapping filter before it is too late by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/307
* merged mappings has to be <= map_mapping_length by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/308
* Revert on disabling tests and adding updated regression output by @pjotrp in https://github.com/waveygang/wfmash/pull/309
* update regression outputs by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/310
* Fix tests by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/313
* Allow big P when mappings are already computed by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/314
* fix MD cigar when using biwfa by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/315
* pad the target when aligning by @ekg in https://github.com/waveygang/wfmash/pull/312
* Group fix by @ekg in https://github.com/waveygang/wfmash/pull/317
* update GitHub Actions workflow to use ubuntu-latest by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/318
* Fix wfa params by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/320
* Scaffold mapping by @ekg in https://github.com/waveygang/wfmash/pull/319
* Scaffold maxi by @ekg in https://github.com/waveygang/wfmash/pull/321
* update WFA2-lib by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/322
* Taskflow map by @ekg in https://github.com/waveygang/wfmash/pull/323
* Taskflow align by @ekg in https://github.com/waveygang/wfmash/pull/324
* Multi index again by @ekg in https://github.com/waveygang/wfmash/pull/326
* Remove unused atomic_queue implementation and includes from aligner by @ekg in https://github.com/waveygang/wfmash/pull/327
* Expose `chain:i` field in the aligned PAF output by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/325

### New Contributors
* @ASLeonard made their first contribution in https://github.com/waveygang/wfmash/pull/269
* @kdm9 made their first contribution in https://github.com/waveygang/wfmash/pull/274
* @pjotrp made their first contribution in https://github.com/waveygang/wfmash/pull/302

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.21.0...v0.22.0

## ChangeLog v0.21.0 (20240902)

Buildable source tarball: [wfmash-v0.21.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.21.0/wfmash-v0.21.0.tar.gz)

Previously, settings that might make runtime slightly better when aligning pangenomes hurt performance in comparative genomics contexts. Updates related to mashmap3 and alignment have made us much more robust to defaults that are more sensitive.

In this release, we're setting a bunch of defaults which have become standard in testing:

* Default minimum mapping identity reduced from 90% to 70%.
* Set maximum mapping length to 50k by default (previously unlimited).
* Changed block length default from 5x segment length to 3x segment length.
* Set default chain gap to 30kb (previously was 6x segment length, up to 30k).
* Reduced default segment length from 5k to 1k.
* Changed default kmer size from 19 to 15.
* Modified wflign to run on all fragments except very small ones (less than 1000 bp).
* Changed filtering logic to use Euclidean distance as an absolute cutoff instead of axis-weighted Euclidean distance, while still ranking based on axis-weighted distance.

These should tend to make wfmash more sensitive at the edges of its performance envelope with minimal costs for easy, low-divergence pangenome alignment problems.

## ChangeLog v0.20.0 (20240826)

Buildable source tarball: [wfmash-v0.20.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.20.0/wfmash-v0.20.0.tar.gz)

### Major Changes

1. **New Global Alignment Approach**:
* Replaced the previous head and tail patching with a comprehensive global alignment strategy.
* Implemented `erode_head` and `erode_tail` functions to remove small, potentially spurious matches at alignment boundaries.
* The alignment now aims to include the entire query sequence, crucial when using the `-P` option for chunking mappings.
* This change ensures continuity across the entire sequence, especially important when mappings are broken into smaller pieces for easier alignment.
* Switched from a semi-global approach (pinned at one end) to a fully global alignment, improving accuracy across the entire sequence length.

2. **Improved Chaining Algorithm**:
* Introduced an axis-weighted Euclidean distance function for more accurate chaining of mappings.
* This new function helps break mappings when encountering large indels, which can be computationally expensive to align.
* Improves detection of large structural variations directly from the mapping stage.
* Reduces spurious chaining in satellite repetitive sequences by considering the diagonal nature of true matches.
* The weighting maintains the original chain gap threshold for on-diagonal matches while effectively shortening the allowed distance for off-diagonal matches.

3. **Mapping and Alignment Improvements**:
* Modified the logic for determining cuttable positions in long alignments to avoid breaking alignments in the middle of structural variations (SVs).
* Adjusted the merging of consecutive mappings to be more selective, prioritizing the preservation of potential SV signals.
* Enhanced the handling of complex genomic structures by improving coordination between mapping and alignment stages.

4. **Performance Optimization**:
* Temporarily disabled multithreaded FASTA input processing due to thread safety issues with the samtools faidx reader.
* This change addresses memory efficiency concerns and prevents potential errors in multi-threaded environments.
* Future updates may reintroduce multi-threaded processing with improved memory management.
* Optimized the mapping process when not splitting sequences.
* Improved efficiency of long mapping handling, particularly when max mapping length is set to infinity.

5. **Default Changes**:
* Changed the default maximum mapping length (`-P/--max-mapping-length`) to infinity, allowing for longer continuous alignments when appropriate.

### Minor Improvements and Bug Fixes

* Enhanced error handling and validation throughout the alignment process.
* Improved coordinate calculations, especially in edge cases involving sequence boundaries and large structural variations.
* Added additional PAF output fields, including a chain identifier for merged mappings.
* Adjusted parameters for more robust alignment in complex regions.

This release significantly improves wfmash's efficiency when handle complex genomic structures (e.g. centromeres) and large-scale variations, particularly when using the `-P` option to chunk mappings for more efficient alignment. While this option has been left unset by default, we do strongly recommend exploring it if you find your alignment times are very slow. A good setting in testing has been `-P50k`.

## ChangeLog v0.19.0 (20240819)

Buildable source tarball: [wfmash-v0.18.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.19.0/wfmash-v0.19.0.tar.gz)

### What's Changed
* It patches by @ekg in https://github.com/waveygang/wfmash/pull/264

This release fixes a bunch of small issues with previous updates to the mapping merging and splitting logic.

The main update should improve mapping coverage by correctly calculating the block length of the mapping based on the pre-split mapping. We also correctly organize cuts to be in regions without SVs.

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.18.0...v0.19.0

## ChangeLog v0.18.0 (20240805)

Buildable source tarball: [wfmash-v0.18.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.18.0/wfmash-v0.18.0.tar.gz)

Improving mapping in complex regions, debugging recursive patching, and other fun.

1. Recursive Inversion Patching:
* Implemented recursive patching for inversions, completing the "multipatch" functionality.
* This allows for more accurate alignment of complex genomic regions with inversions.

2. SAM Output for Multipatch Alignments:
* Added support for SAM output format for multipatch alignments.
* Ensures consistent representation of complex alignments across different output formats.

3. Orientation-Consistent Alignments:
* Improved alignment consistency across all orientations of reference-query pairs.
* Enhances reliability and reproducibility of alignment results.

4. Optimized Inversion Patching:
* Implemented a bound on the maximum score for inverted patches.
* Allows for early termination of alignment when the inverted patch is worse than the forward alignment.

5. Dynamic Multi-Producer Alignment Module:
* Rewrote the alignment module to support multiple producers filling the work queue.
* Dynamically handles memory issues, improving efficiency and scalability.

6. Overlap Filtering in Plane Sweep Algorithm:
* Implemented an overlap filter to prevent keeping suboptimal mappings.
* New CLI option: `-O, --overlap-threshold <F>`
* Allows setting the fraction F for dropping mappings overlapping with higher scoring mappings.
* Default value is 0.5.

7. Long Mapping Fragmentation:
* Enabled breaking of long mappings into smaller fragments at junction points.
* Junctions are defined by four consecutive segments, allowing for more precise breakpoint detection around structural variations.
* New CLI option: `-P, --max-mapping-length <N>`
* Sets the maximum length of a single mapping before breaking.
* Default value is 1M (1 million bases).

8. Improved Handling of Satellite Sequences:
* The combination of overlap filtering, mapping fragmentation, and recursive patching significantly improves wfmash's ability to handle satellite sequences.
* These changes address common performance issues and mapping problems associated with highly repetitive regions.
* Users should expect better accuracy and efficiency when aligning genomes with abundant satellite sequences.

9. Performance Improvements:
* Various optimizations and code refactoring for better overall performance.

10. Bug Fixes and Minor Enhancements:
* Multiple bug fixes and small improvements throughout the codebase.

This release significantly enhances wfmash's ability to handle complex genomic structures, including challenging satellite sequences. It improves output consistency and optimizes performance for large-scale alignments. The new features and CLI options provide more accurate and detailed alignment information, particularly for regions with inversions, structural variations, and repetitive elements, while offering users greater control over the alignment process. These improvements make wfmash more robust and efficient for a wider range of genomic analyses, especially those involving highly repetitive or complex regions.

### What's Changed
* speed up inversion patching by @ekg in https://github.com/waveygang/wfmash/pull/260
* Break mappings by @ekg in https://github.com/waveygang/wfmash/pull/261

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.17.0...v0.18.0

## ChangeLog v0.17.0 (20240720)

Buildable source tarball: [wfmash-v0.17.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.17.0/wfmash-v0.17.0.tar.gz)

This release introduces multipatch alignment capabilities, significantly enhancing wfmash's ability to handle complex genomic structures, particularly inversions and other rearrangements. Multipatching refers to a process in which the initial wflign traceback is patched, we determine that an inverted orientation of the patch is preferable (as introduced in v0.16.0), and (in v0.17.0) we now attempt multiple patching steps to span the gap. Key improvements include:

**Multipatch Alignment**:
* Implemented a progressive alignment approach that can detect and align multiple patches, including inversions, within a single alignment region.
* Added a new tag `patch:Z:true` to indicate multipatch alignments in the output.
* Introduced an `inv:Z:true/false` tag to specify whether a patch is inverted.

**Alignment Refinements**:
* Implemented trimming of alignments to remove leading and trailing indels, improving alignment quality.
* Added bounds detection for alignments to better handle partial matches.
* Increased the default chain gap to 6x segment length or 30k, allowing for detection of larger variants.

**Output Enhancements**:
* Modified the output format to clearly distinguish multipatch alignments.
* Improved logging and debugging output for better insight into the alignment process.

**Code Improvements**:
* Enhanced the `alignment_t` class with new accessors for query and target begin/end positions.
* Implemented pruning of overlapping patches to avoid redundant alignments.
* Refactored several core functions for better modularity and readability.

**Build System**:
* Added libdeflate as a dependency in the Guix build configuration.

This release significantly improves wfmash's ability to handle complex genomic alignments, particularly those involving local inversions and other structural variations. The multipatch approach allows for a more complete representation of genomic relationships in challenging regions than is available in other methods.

Happy aligning with enhanced structural variation breakpoint resolution! 🧬🔍🧮

### What's Changed
* add deflate to guix.scm by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/258
* Multi-patch by @ekg in https://github.com/waveygang/wfmash/pull/259

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.16.0...v0.17.0

## ChangeLog v0.16.0 (20240702)

Buildable source tarball: [wfmash-v0.16.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.16.0/wfmash-v0.16.0.tar.gz)

The primary enhancement in this release is the implementation of inversion detection during the alignment patching process. This feature significantly improves the alignment accuracy for sequences containing inversions.

#### How it works:

1. **Patching Process**: During the wflign high-level trace patching, the algorithm identifies regions that do not align well in the forward orientation.

2. **Reverse Complement Alignment**: For these poorly aligned regions, the algorithm attempts an alignment with the reverse complement of the sequence.

3. **Score Comparison**: The algorithm compares the alignment scores of the forward and reverse complement alignments.

4. **Selection**: If the reverse complement alignment produces a better score, it is selected for that region.

5. **Output**: Reverse complement alignments are reported with an additional SAM tag `rc:Z:true`.

#### Key Components:

* New parameter `wflign_min_inv_patch_len`: Sets the minimum length of an inverted patch to be considered (default: 23).
* `calculate_alignment_score` function: Computes alignment scores based on the CIGAR string and penalties.
* Modified `do_wfa_patch_alignment` function: Now handles both forward and reverse complement alignments.
* Updated `write_merged_alignment` function: Processes and outputs reverse complement alignments.

This feature allows wfmash to accurately align sequences with inversions, improving its utility for complex genomic comparisons.

### Other Significant Changes

1. **MashMap Index Support**:
* Implemented creation and usage of MashMap indexes for faster repeat mapping.
* New CLI options: `--mm-index`, `--create-index-only`, `--overwrite-mm-index`.

2. **Memory Optimization**:
* Improved memory usage in the `Sketch` class.

3. **Kmer Size Calculation**:
* Modified to handle edge cases with high-identity alignments.

4. **Alignment Class Improvements**:
* Enhanced `alignment_t` class with proper copy and move semantics.

5. **Index File Handling**:
* Improved reading and writing processes with parameter validation.

### Detailed Log of Changes

#### src/align/include/align_parameters.hpp
* Added `wflign_min_inv_patch_len` parameter to `Parameters` struct.

#### src/align/include/computeAlignments.hpp
* Integrated `wflign_min_inv_patch_len` into `WFlign` constructor call.

#### src/common/wflign/src/wflign.cpp and wflign.hpp
* Added `min_inversion_length` to `WFlign` constructor and member variables.
* Modified `minhash_kmer_size` calculation for edge cases.

#### src/common/wflign/src/wflign_alignment.cpp and wflign_alignment.hpp
* Implemented copy/move constructors and assignment operators for `alignment_t`.
* Added `calculate_alignment_score` function.

#### src/common/wflign/src/wflign_patch.cpp and wflign_patch.hpp
* Modified `do_wfa_patch_alignment` for reverse complement handling.
* Updated `write_merged_alignment` for reverse complement output.
* Refined patching process for bidirectional alignment consideration.

#### src/interface/parse_args.hpp
* Added CLI options for MashMap indexing and `wflign_min_inv_patch_len`.

#### src/map/include/map_parameters.hpp
* Added parameters for MashMap indexing support.

#### src/map/include/parseCmdArgs.hpp
* Updated parsing for new MashMap indexing options.

#### src/map/include/winSketch.hpp
* Implemented MashMap index functions (create, read, write).
* Added CLI-index file parameter validation.
* Optimized `Sketch` class memory usage.

## ChangeLog v0.15.0 (20240617)

Buildable Source Tarball: [wfmash-v0.15.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.15.0/wfmash-v0.15.0.tar.gz)

Initial experiments in our all-to-all alignment of the draft [vertebrate genomes project](https://vertebrategenomesproject.org/) demonstrated that we were not generating end-to-end alignments for many mashmap3 homology pairs at 70\% ANI (`wfmash -m -p 70`). Exploration showed that our attempts at automatically tuning alignment parameters based on mashmap estimated identity simply didn't work. The parameter settings we used meant that optimal wflign alignments were often `I*D*`, or "fully indel-ed", leading to no insight into the homology between the pairs even when internally WFA segments did match.

To avoid this "gotcha" and ensure we obtain an alignment, we set the softest wflign parameters possible to maintain the inequality `match < gap-extend < mismatch < gap-open`: match=0 mismatch=2, gap-open=3, gap-extend=1. We also use 0,3,4,2,24,1 for our WFA patching parameters, matching minimap2's asm20 setting. These changes lead to a major improvement in runtime and memory usage during alignment. In WFA, where everything is order of score or score*score, smaller scores mean lower memory and faster runtime.

We also ran into portability issues. The biggest improvement was to bring back static builds with options to enable generic compatibility with many recent x86 systems. This will allow direct distribution of binaries in these releases.

We also hit some very weird software bugs that led us to drop jemalloc. It was causing very strange problems (like IOT like invalid instruction errors, signal 9 allocation errors with 5% RAM usage, etc.) and offers no obvious performance advantage in wfmash's current setup, mentioning here because it was a very tricky bug to resolve.

### New Features and Enhancements

#### Breaking Changes

`wfmash` now requires the query FASTA sequence to be bgzipped and samtools faidx indexed as well as the target sequence. This lets us basically be able to randomly access the query which improves performance in parallel and high-performance computing settings because we don't have to spool through very big query files if we're only aligning a very small part of them.

#### Publications
* Added a new citation for the biWFA algorithm:
* **Santiago Marco-Sola, Jordan M. Eizenga, Andrea Guarracino, Benedict Paten, Erik Garrison, and Miquel Moreto**. ["Optimal gap-affine alignment in O (s) space"](https://doi.org/10.1093/bioinformatics/btad074). *Bioinformatics*, 2023.

#### Build System
* **Configurable Build Options**: Introduced new CMake options to make the build process more flexible:
* `BUILD_STATIC`: Option to build a static binary.
* `BUILD_DEPS`: Option to build external dependencies (htslib, gsl, libdeflate) from source.
* `BUILD_RETARGETABLE`: Option to build a retargetable binary without machine-specific optimizations.
* **Static Compilation**: Improved support for static compilation, including the ability to build static binaries and handle external dependencies more flexibly.
* **OpenMP Support**: Added OpenMP support for parallel processing.
* **Improved Documentation**: Updated the README to provide detailed instructions for building from source, including static and retargetable binaries.

#### Performance and Optimization
* **Optimized Compilation Flags**: Adjusted compilation flags for better performance and compatibility across different systems.
* **Memory Management**: Improved memory management by reducing the number of sketches kept in memory during large alignments.
* **Query Sequence Handling**: Enhanced the handling of query sequences to support random access, reducing memory usage and improving performance.

#### Bug Fixes
* **Memory Access Errors**: Fixed potential memory access errors by adding bounds checks for sequence indices.
* **Thread Safety**: Ensured thread safety by using a single `faidx_t` object for sequence fetching, shared among multiple threads.
* **Alignment Filtering**: Disabled low-identity filtering by default to ensure all alignments are kept for post-processing.

#### Miscellaneous
* **Nix and Guix Support**: Added support for building wfmash using Nix and Guix, including Docker image generation.
* **Test Cases**: Added a script to generate test cases for wflign, facilitating easier testing and validation.

### Detailed Changes

#### Commit Highlights
* **Commit 577c3de**: Added biWFA citation to the README.
* **Commit 1d142d9**: Merged changes for Stampede3 build configuration.
* **Commit d55cfe7**: Made the build configurable and documented how to use the new options.
* **Commit 18e33b0**: Fixed the path for libdeflate in the CMake configuration.
* **Commit 9ff0452**: Merged updates for scoring parameter optimizations.
* **Commit 609082b**: Updated build to use Clang and removed jemalloc dependency.
* **Commit e6f1824**: Restored micromamba/anaconda support.
* **Commit debeff7**: Debugged build on TACC's Stampede3 cluster.
* **Commit 75a6631**: Improved build process for Stampede3 cluster.
* **Commit 719381c**: Avoided `-march=native` for broader compatibility.
* **Commit 081213c**: Fixed memory management issues in alignment code.
* **Commit fb4c6d0**: Used generic modern optimizations, avoiding processor-specific flags.
* **Commit c04088e**: Ensured zero-termination of sequence data fetched with `faidx_fetch_seq64`.
* **Commit ea76722**: Added validation for mashmap input rows.
* **Commit fedad55**: Reduced the number of sketches kept in memory during large alignments.
* **Commit a7aa342**: Improved queue behavior and memory management in alignment code.
* **Commit 581b364**: Disabled low-identity filtering by default.
* **Commit 3f1f7af**: Corrected documentation of queues.
* **Commit acd7fdc**: Updated atomic queue definition for better single-producer multi-consumer behavior.
* **Commit 2b91145**: Avoided deadlock on empty input files.
* **Commit edcd281**: Used a single `faidx_t` object for sequence fetching to save memory.
* **Commit a04908b**: Fixed scoring parameters for diverse alignment problems.
* **Commit bb0d43d**: Merged updates for forcibly using biWFA alignment.
* **Commit 3e434f8**: Added a script to create test cases for wflign.
* **Commit cc89be6**: Added option to force global biWFA alignment.
* **Commit 35194d8**: Merged updates for random access to queries during alignment.
* **Commit c738c1d**: Stopped sorting the input mapping file for better performance.
* **Commit 70e896a**: Removed redundant query sequence processing.
* **Commit 6fcddc2**: Enabled random access of query subsequences in alignment.
* **Commit d9a0880**: Limited to one query sequence file for simplicity.
* **Commit 729d9d7**: Merged updates for static build reimplementation.
* **Commit 8477337**: Corrected debugging build with PNG and TSV support.
* **Commit efc8f04**: Added libdeflate as a dependency.
* **Commit deb1472**: Updated minimum CMake versions.
* **Commit d1588e6**: Described static compilation options in the README.
* **Commit 2713899**: Defaulted to non-static builds.
* **Commit a96919e**: Reimplemented static builds.
* **Commit 899e154**: Bumped Nix build configuration.
* **Commit 7376468**: Reverted removal of `flake.lock`.
* **Commit 526995a**: Removed `flake.lock`.
* **Commit f94b7ee**: Updated Nix build configuration.
* **Commit e2df9c8**: Moved to Nix flake.
* **Commit cbedc8f**: Locked the Nix flake.
* **Commit b0d0ada**: Added Nix flake configuration.

Happy whole-genome-aligning! 🔬🧬📊

## ChangeLog v0.14.0 (20240529)

Buildable Source Tarball: [wfmash-v0.14.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.14.0/wfmash-v0.14.0.tar.gz)

This release provides support for subsetting the queries which are used in addition to the target subsetting. A list of queries can be offered. (We still work with only a single target though.) The idea is that this will make it possible for us to subdivide the all-versus-all alignment matrix and run many small jobs where multiple queries are aligned against a single target. However, running all queries against one target would be computationally infeasible, because there might be many hundreds of thousands of queries. There are some other bug fixes and updates as well, but the main difference that triggers a release is the change in the command line API.

### changelog

#### Query filtering and specification improvements
* Added support for specifying a comma-delimited list of query name prefixes to filter queries with the `-Q`/`--query-prefix` option.
* Added `-A`/`--query-list` option to specify a file containing a list of query sequence names to use.
* Updated internal sequence iteration and counting logic to properly apply the new query filtering options.

#### Target filtering option name changes
* Renamed target prefix filtering option from `-P`/`--target-prefix` to `-T`/`--target-prefix` for consistency.
* Renamed target list filtering option from `-A`/`--target-list` to `-R`/`--target-list`.

#### All-to-all alignment script improvements
* Updated `scripts/all2all_jobs.py` to:
* Support grouping by genome, haplotype, or contig.
* Allow specifying different grouping levels for target and query sequences.
* Directly generate wfmash command lines.
* Added `scripts/make_source_targball.sh` to generate a source tarball for releases.

#### Build and testing updates
* Added back `rt` library to CMake configuration.
* Updated CI tests to run on the `main` branch.
* Adjusted CI test cases for the subset of the LPA dataset.

#### Bug fixes
* Fixed a `heap-use-after-free` error in `wflign_affine_wavefront()`.

## ChangeLog v0.13.1 (20240430)

Buildable Source Tarball: [wfmash-v0.13.1.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.13.1/wfmash-v0.13.1.tar.gz)

### What's Changed
* drop timing from the output by default by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/231
* fix "alignment block length" in PAF output by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/235
* improve MAC compatibility by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/236
* Make chain gap dynamic by @bkille in https://github.com/waveygang/wfmash/pull/234
* for the chain_gap, do not go above 20k by default by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/238

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.13.0...v0.13.1

## ChangeLog v0.13.0 (20240323)

Buildable Source Tarball: [wfmash-v0.13.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.13.0/wfmash-v0.13.0.tar.gz)

### What's Changed
* Do not allow a segment to chain with itself by @bkille in https://github.com/waveygang/wfmash/pull/225
* Convex penalties for the alignment patching by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/229

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.12.6...v0.13.0

## ChangeLog v0.12.6 (20240205)

Buildable Source Tarball: [wfmash-v0.12.6.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.12.6/wfmash-v0.12.6.tar.gz)

### What's Changed
* Add error when the index is empty  by @bkille in https://github.com/waveygang/wfmash/pull/219
* Don't include endpoint in L2 range by @bkille in https://github.com/waveygang/wfmash/pull/220
* make the script more robust by @AndreaGuarracino in https://github.com/waveygang/wfmash/pull/222

**Full Changelog**: https://github.com/waveygang/wfmash/compare/v0.8.4...v0.12.6

## ChangeLog v0.12.5 (20240124)

Buildable Source Tarball: [wfmash-v0.12.5.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.12.5/wfmash-v0.12.5.tar.gz)

This closes an integration process that pulled in minmer-based mashmap3. The CLI has been updated to make it easy to run large all-to-all alignment jobs in parts on a cluster using the `-P[pfx], --target-prefix=[pfx]` and target list options. We have included a new `--one-to-one` mapping option that ensures symmetry between mapping query and target ranges (they must each be the best option).

## ChangeLog v0.10.5 (20230802)

Buildable Source Tarball: [wfmash-v0.10.5.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.10.5/wfmash-v0.10.5.tar.gz)

This release fixes the issue of getting N exact mappings when less than N mappings are requested. Moreover, now it is possible to specify a prefix, or a list of sequences, to filter target sequences. Finally, `wfmash` building with `cmake` is improved.

This introduces:
* fix the build to use cmake everywhere #180 #181 by @ekg
* avoid exact mapping duplicates #182 by @ekg
* subset targets by list or prefix #183 by @ekg

## ChangeLog v0.10.4 (20230627)

Buildable Source Tarball: [wfmash-v0.10.4.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.10.4/wfmash-v0.10.4.tar.gz)

This release strongly improves the runtime performance of the alignments without affecting their accuracy (in all alignments tested). We apply an heuristic BiWFlign/BiWFλ and the mapping chaining is less aggressive, improving also the runtime during the mapping phase with the price of getting a bit more mappings.

This introduces:
* make gcc 11.4.0 and 13.1.0 happy #177 by @AndreaGuarracino
* update WFA2-lib #179 by @AndreaGuarracino
* Heuristic BiWFλ and less aggressive mapping chaining #176 by @AndreaGuarracino

## ChangeLog v0.10.3 (20230419)

Buildable Source Tarball: [wfmash-v0.10.3.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.10.3/wfmash-v0.10.3.tar.gz)

This release reimplements the bidirectional `wfλ`, with a very tiny, but consistent, performance improvement.

This introduces:
* Fix versioning #170 by @AndreaGuarracino
* Refactored BiWFλ (BiWFlambda) #169 by @smarco @AndreaGuarracino

## ChangeLog v0.10.2 (20230113)

Buildable Source Tarball: [wfmash-v0.10.2.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.10.2/wfmash-v0.10.2.tar.gz)

Minor change. This release is used to update the conda release.

This introduces:
* executable scripts #167 by @AndreaGuarracino

## ChangeLog v0.10.1 (20221217)

Buildable Source Tarball: [wfmash-v0.10.1.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.10.1/wfmash-v0.10.1.tar.gz)

Minor changes. This release is used to update the conda release.

This introduces:
* always emit long version format #164 by @AndreaGuarracino
* add `-O/--invert-filtering` parameter #165 by @AndreaGuarracino

## ChangeLog v0.10.0 (20221027)

Buildable Source Tarball: [wfmash-v0.10.0.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.10.0/wfmash-v0.10.0.tar.gz)

The main change is the introduction of the bidirectional `wfλ`, which allows a full exploration of the high-order alignment matrix without blowing up with the memory. This leads to more and better anchors that tend to simplify the alignments during the patching phase, therefore even improving the runtime.

This introduces (in reverse chronological order):
* WFA2-lib update #148 #149 by @smarco @AndreaGuarracino
* bug fixes in reduction #151 #152 #155 by @AndreaGuarracino
* implementation of `BiWFλ`, bidirectional flavor of `wfλ` #153 by @AndreaGuarracino
* better head and tail patchings #158 by @AndreaGuarracino
* fix trace plot for short/medium sequences #159 by @AndreaGuarracino
* tune parameter for high divergence #160 by @AndreaGuarracino

## ChangeLog v0.9.2 (20220801)

Fixes several bugs (including a memory leak), updates WFA2-lib and documentation.

## ChangeLog v0.9.1 (20220703)

Buildable Source Tarball: [wfmash-v0.9.1.tar.gz](https://github.com/waveygang/wfmash/releases/download/v0.9.1/wfmash-v0.9.1.tar.gz)

Major changes in this release are an improvement of mapping merging and updates to the base level alignment code which allow for exact patching of large gaps in the alignment traceback.

Fixes and updates:

* signed 32-bit overflow on large plant chromosomes
* avoid overfiltering due to equal scores on multiple targets
* mapping min-hash sparsification (hash mapping struct and keep hashes below a threshold) to sample from very large all-vs-all alignment problems
* simplify mapping merging and cap the merge distance at the wavefront reduction distance used in wflign
* reduce wavefronts while following the diagonal, to handle wflign of sequences with very different lengths without alignment failure
* universal minimizer spaced seeds: ALeS seeds now work with unbiased universal minimizer scheme
* integration of WFA2-lib and application of biWFA for `O(s)`-memory patching of large gaps in alignment traceback

## ChangeLog v0.8.2 (20220413)

Buildable Source Tarball: [wfmash-v0.8.2.tar.gz](https://github.com/ekg/wfmash/releases/download/v0.8.2/wfmash-v0.8.2.tar.gz)

This introduces:

* updates in how `wfmash` is compiled/built to ensure greater inter-system compatibility;
* adaptive penalties for the alignment, with more permissive wflambda/wflign parameters.

## ChangeLog v0.8.1 (20220328)

**Buildable Source Tarball:** [wfmash-v0.8.1.tar.gz](https://github.com/ekg/wfmash/releases/download/v0.8.1/wfmash-v0.8.1.tar.gz)

This introduces:
* fixed a bug in mapping filtering for short sequences;
* default segment size (`-s`) at 10 kbps;
* fixed alignment penalties regardless of the requested mapping identity (`-p`): this strongly reduces the runtime and lead to much more compressed representations of the alignments between sequences.

## ChangeLog v0.8.0 (20220323)

**Buildable Source Tarball:** [wfmash-v0.8.0.tar.gz](https://github.com/ekg/wfmash/releases/download/v0.8.0/wfmash-v0.8.0.tar.gz)

`wfmash` is now substantially better at mapping and alignment at very high sequence divergences. This involves many changes relative to v0.7.0.

### mashmap3

The mapping module has been largely rewritten to allow for mappings to span large structural variation. We now apply multiple merging passes in 2D over the query/target mapping matrix (mashmap2 used a 1D approach in the query). The first unites mappings found within 2x the segment length (`wfmash -s`). Subsequently, multiple rounds of greedy merging and plane-sweep filtering merge the closest mappings on a near diagonal within a given chaining gap (`wfmash -c`). We finally filter the mappings at 5x segment length (`wfmash -l`) rather than 3x in previous releases.

The updated mapping merging also allows us to make a sparser first mapping step, as segment mapping drop-outs can be spanned using this approach. This allows us to use relatively sparse minimizer selection, which reduces the number of candidate (usually erroneous) mappings to consider.

We have also applied [world minimizers, which are unbiased](https://www.biorxiv.org/content/10.1101/2022.01.14.476226v1.abstract) and faster to compute than window minimizers. To ensure efficient performance, we implement a much stronger filter on repetitive minimizers, filtering out the top 0.5% of most-frequent minimizers, which is now configurable with `wfmash -H`.

### divergence-adaptive wflign

This release also features improvements to the base level alignment that are essential for sensitive alignment at high divergence. We now rest more heavily on the wflign matrix, which leads to a more complete exploration of alignment possibilities. Alignment parameters---such as dynamic programming scoring (for WFA), maximum sketch distance to evaluate a local alignment, and max allowed alignment score---are now set based on a function of the mashmap-based identity.

### testing approach

To develop this release, we tested on sequence collections with up to 30% divergence. We ensured that adjustments worked on a series of test cases drawn from humans, yeast, e. coli, potato, and fish, including a [scale-up test to an all-vs-all alignment of 45 fish assemblies](https://twitter.com/erikgarrison/status/1497502046668435457).

### user considerations

In contrast to previous versions, `wfmash` v0.8.0 is less sensitive to particular segment length settings. The meaning of `-p`, or the minimum pairwise identity of the mappings, is also somewhat softened, because mappings can now span very large gaps, up to `--chain-gap` which defaults to 100x the segment length. Very long segment lengths of 50-100kb are probably less necessary, and we're seeing good performance at 5kb to 20kb segment lengths.

The increase in the minimum mapping length filter (from 3x to 5x segment length) reflects increased sensitivity and also potential errors caused by these changes.

An additional concern is that users seeking to map against extremely repetitive sequences may need to set `-H` lower. Increasing `-s` can also span gaps caused by repeats and derive alignments for them. Alignments that focus strongly on repetitive regions may still need special parameter tuning. The default settings are now focused on obtaining reasonable homology maps for pangenome and pan-clade alignment problems.

### visualization of wflign alignment matrix

Parameter tuning was assisted with visualizations of the wflign (high-order, over 256bp `wfmash -W`-length segments) alignment matrix. These show regions compared using kmer jaccards in gray, attempted successful alignments in green, and blue for failed alignments.

Two 1Mbp regions of yeast genomes:

![image](https://user-images.githubusercontent.com/145425/159750315-0115bf3c-9ced-4bb7-bd4e-87e44619ec42.png)

... and the full alignment matrix (pafplot):

![image](https://user-images.githubusercontent.com/145425/159750340-4a0bb88c-683c-452b-a47b-3c099b41380b.png)

Two fish chromosomes at ~25% pairwise divergence in aligned regions (`wfmash -p 70 -s 20k`).

![image](https://user-images.githubusercontent.com/145425/159750448-60a848fe-fec9-496e-b959-abf1235f567d.png)

And a few alignments through human lipoprotein A (LPA):

![image](https://user-images.githubusercontent.com/145425/159750828-929b3c81-d303-4555-bcd6-1458c180f403.png)

![image](https://user-images.githubusercontent.com/145425/159750858-578d3f5a-a7d8-4df5-96f1-a98ade16c56c.png)

![image](https://user-images.githubusercontent.com/145425/159750880-aae2ac8f-7072-4773-a69e-30382cb49637.png)

## ChangeLog v0.7.0 (20210909)

**Buildable Source Tarball:** [wfmash-v0.7.0.tar.gz](https://github.com/ekg/wfmash/releases/download/v0.7.0/wfmash-v0.7.0.tar.gz)

This release introduces a huge amount of updates:
* the mapping parameters (window size and kmer size) are adaptive with respect to the requested segment identity;
* the alignment parameters (mismatch/gap penalties and the max mash distance heuristic) are adaptive with respect to the estimated identity for each mapping region;
* WFA was updated to the last WFAv2, which includes important memory usage optimizations;
* wflign / wflambda are upgraded to WFAv2, leading to a strong reduction of the memory usage;
* alignment accuracy is improved during the patching, thanks to the reduced memory usage of the new WFAv2;
* robin-hood structures are applied, to improve runtimes;
* matches and (part of the) mis-matches are cached in wflign, improving the runtime by paying little memory overhead;
* pure-WFA alignment is performed for short sequences (and short mapping regions in long sequences);
* ends-free WFA for head/tail patching, replacing edlib;
* fixed a reduction bug in the WFA library;
* input PAF from other aligners are supported;

## ChangeLog v0.6.1 (20210726)

**Buildable Source Tarball:** [wfmash-v0.6.1.tar.gz](https://github.com/ekg/wfmash/releases/download/v0.6.1/wfmash-v0.6.1.tar.gz)

This (little) release includes:

* handy parameters (#89);
* a buildable source tarball;
* a little compiling fix.

## ChangeLog v0.6 (20210716)

Here, we sparsify the wflign problem, and then patch through the gaps using a low-memory version of WFA (cheers @smarco !)

## ChangeLog v0.5 (20210517)

A number of changes in wfmash have completed the alignment patching in wflign, rendering it stable and memory-thrifty enough to safely apply to large genomes. The mapping in general has also been improved by targeting a smaller windowSize parameter, and capping it at 256 to not generate confusion when mapping large segments.

## ChangeLog v0.4 (20210422)

With this version, alignments are patched with WFA (for unaligned regions where the short axis is up to 8kb) and edlib (for very short unaligned stretches). Edlib in semiglobal mode is used to patch up the heads and tails of the alignments. Previous versions have significant dropouts in alignments, but with these changes the issue is largely resolved.

## ChangeLog v0.3.1 (20210115)

This point release updates wflign to emit a single merged alignment for each mapping. The output is compact and ready for eventual adaptation to SAM output.

## ChangeLog v0.3 (20210111)

`wfmash` is now sync'ed with `edyeet` and an update to `wflign` lets us use WFA to obtain base-level alignment with affine gap costs. This is more biologically plausible than edit-distance based alignment provided by edilb.

Alignment runtime increases by 2-3x, depending on divergence rate given by `-p[%], --map-pct-id=[%]`, with higher thresholds experiencing lower relative slowdown.

`wfmash` uses both wavefronts and mash distance (locality sensitive hashing) in two contexts. For mapping, it uses MashMap2's algorithm. For base-level alignment, it uses `wflign`, which is WFλ with λ = WFA guided heuristically with mash distance.

## ChangeLog v0.2 (20201019)

Other people can compile it.

And it's generally improved since v0.1.

## ChangeLog v0.1 (20200923)

I'm still seeing problems with memory allocation in WFA, but this commit does produce apparently valid alignments!

So it might be helpful to others, and it's reasonably safe to use :stuck_out_tongue_winking_eye:

Happy aligning!

### mappingFilter skip-prefix grouping fix (merged into `main`)

* Fixed the reference-group subrange scan in the mapping filter
  (`src/map/include/mappingFilter.hpp`, active with `--skip-prefix`):
  the scan treated reference names *without* the prefix delimiter as
  one hardcoded shared group (a `group()` lambda returning 0), instead
  of asking the sequence index for each sequence's actual group.
  Since `getRefGroup()` returns the per-sequence `groupId` computed at
  index time, the hardcoded-0 shortcut could mis-group delimiter-less
  names and extend filter subranges across reference groups, corrupting
  the per-group top-`n` selection.  The scan now captures the first
  mapping's group and compares candidates against
  `idManager.getRefGroup()`, as before the regression.

### Low-divergence engine (merged into `main`)

* One binary, two engines: the **0.24 engine runs by default**;
  `wfmash --engine low-divergence ...` opts in to the 0.14 lineage.
* New `--engine` CLI concept: `src/engine/NAME/` hosts complete,
  self-contained engines with their own parameter parsing;
  `wfmash --engine NAME ...` dispatches the raw argv to the engine
  before the mainline parser runs.
* First engine: `low-divergence` — the wfmash 0.14 lineage (mashmap
  v3.1.1-era mapping plus the wflign aligner), now vendored from the
  maintained **v0.14.1 branch head (9b2a7388)** — the revision pggb
  master pins — into `src/engine/low-divergence/` with `lde_`-prefixed
  files and `lde_`-prefixed namespaces.  Mapping and alignment output
  are byte-identical to the wfmash-0.14.1-0.9b2a7388 binary on yeast
  all-vs-all (1,622 mappings), MHC map+align at `-p 95` and at the
  default cutoff, with clean exits; at workshop scale (623 alignments
  over 15 primate assemblies) the v0.14.1 branch and the 7bf8988
  snapshot the engine originally vendored produce byte-identical
  output, so the CWL regression suite is unaffected.
* The v0.14.1 re-vendor carries the branch's fixes on top of the
  snapshot: WFA2 aligner reuse across records (`thread_local`), lazier
  traceback handling in wflign patching, the k-mer overfiltering fix
  and the unaligned-load fix.  `agc_index.hpp` is vendored again
  (`lde_agc_index.hpp` — v0.14.1's seqiter includes it); without
  `WFMASH_HAVE_AGC` it compiles down to the `is_agc_file()` stub.
* WFA2-lib is vendored into the engine
  (`src/engine/low-divergence/deps/lde_WFA2-lib`, at the v0.14.1
  revision — it now carries the new `cigar_utils` API; the renamed
  symbol list grew from 430 to 437) with every linked C symbol renamed
  `lde_` and the C++ binding namespace renamed `wfa`→`lde_wfa`.  The
  mainline WFA2 revision differs and produces different alignment
  boundaries, and sharing one copy caused symbol/ODR collisions
  between the two engines (which also corrupted the mainline aligner
  through merged inline symbols).
* The engine no longer depends on mainline helpers
  (`wfmash::handy_parameter`, `wfmash::is_a_number`); it carries its
  own copies.  A `lde_standalone` test target builds the engine alone,
  with no mainline sources linked.
* `lde_main` no longer falls off the end without a return statement:
  legal for `main()` (implicit `return 0`), undefined behaviour after
  the rename — GCC turned the cleanup path into an `ud2` trap, crashing
  every engine alignment run after it had written its complete,
  correct output.
* Portability fix: added missing `<cmath>` includes in the vendored
  wflign code (GCC 14).
