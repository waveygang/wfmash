# Stampede3 Mapping Truncation Fix Summary

## Root Cause
The issue ONLY occurs on Stampede3 (and likely other HPC systems with parallel file systems like Lustre/GPFS). The problem is **concurrent file I/O on parallel file systems**.

## Why It Fails on Stampede3
1. **Parallel file systems have different consistency models** - Thread-local FAIDX readers see inconsistent file states
2. **High metadata latency** - FAIDX index operations are very slow on Lustre/GPFS
3. **Lock contention** - Multiple threads accessing the same file causes issues
4. **Race conditions** - File handles and metadata aren't properly synchronized across threads

## Why Other Approaches Failed
1. **Nested subflows** - Taskflow can't properly track dynamically created nested tasks
2. **Flattened structure** - Still has concurrent file I/O issues on parallel file systems
3. **Thread-local readers** - Lifetime issues + parallel file system inconsistency
4. **wait_for_all()** - Doesn't fix the underlying file I/O race conditions

## The Solution That Works
**Pre-load all sequences** (from high-concurrency-shuo branch)
- Sequential file I/O eliminates all race conditions
- No concurrent access to parallel file system
- Simple task structure without nested subflows
- Works reliably on all systems including Stampede3

## Why Pre-loading Is Necessary for HPC
- It's NOT about the task structure - it's about the file system
- Parallel file systems are optimized for large sequential I/O, not random access
- Thread-local readers + parallel file system = recipe for disaster
- Pre-loading is standard practice in HPC for this exact reason

## Clean Implementation
Remove the semaphore from Shuo's branch - it's unnecessary. The pre-loading itself fixes the issue.

## Memory Impact
- Acceptable for most use cases (3GB for human genome)
- Can add batching for extremely large datasets if needed
- Much better than missing mappings!

## Conclusion
This is a **system-specific issue** caused by parallel file system behavior. Pre-loading sequences is the correct solution for HPC environments. The nested subflow issues were a red herring - the real problem was always concurrent file I/O on Stampede3's parallel file system.