# wfmash v0.24.0 Release

This release brings significant memory optimizations, improved mapping scaffolding capabilities, and enhanced ANI-based identity estimation.

## Major Improvements

### Memory Optimization
- **Drastically reduced memory usage** during mapping phase (~66% reduction)
- Optimized alignment phase with on-demand record loading
- Clean memory separation between mapping and alignment phases
- Optimized ANI sketching phase memory consumption
- New compact mapping structures for better memory efficiency

### Mapping Scaffolding
- **New 2D distance graph scaffolding algorithm** for improved syntenic block detection
- Enhanced scaffold filtering with plane sweep optimization
- Configurable minimum scaffold length (default: 5kb)
- Support for scaffold mapping output via new options
- Better handling of boundary mappings for improved leniency

### ANI-based Identity Estimation
- **New ANI preset system** (ani25, ani25-5, etc.) for automatic identity threshold selection
- Automatic identity estimation with `-p auto`
- Streaming MinHash implementation for efficient ANI computation
- Parallel ANI estimation with TaskFlow
- Per-group identity calculations with better CPU utilization

## New Features

### Build System
- Added `VENDOR_HTSLIB` CMake option for building without system htslib
- Updated WFA2-lib submodule integration
- Improved build optimization flags

### Command Line Interface
- Redesigned CLI parameters for better usability
- Changed sketch parameter to `-s` (was `-S`)
- Changed window-size parameter to `-w` (was `segLength`)
- Updated default overlap threshold from 1.0 to 0.95
- Minimum L1 hits now defaults to 3 (configurable with `-H`)
- Map sparsification parameter for controlling mapping density

### Performance
- ~25% speedup for small genomes through optimized reverseComplement function
- Per-group mutexes for better parallel scaling
- Thread-local reader functions for improved I/O performance
- Progress reporting for all pipeline phases

## Bug Fixes
- Fixed critical bug in MinHash sketch computation for groups
- Resolved `stoi` conversion errors with invalid records
- Fixed type conversion overflow in parameter handling
- Corrected mapping merge logic for query and reference spans
- Fixed boundary mapping criteria for better edge case handling

## Technical Details
- Maintained chain identity information for `ch:Z:` tag
- Helper function to merge adjacent CIGAR operations
- Improved sequence loading patterns for ANI estimation
- Better progress reporting throughout all phases

## Contributors
Thanks to all contributors who made this release possible, with special mentions to those who worked on memory optimization, scaffolding improvements, and the ANI estimation system.