# Mapping Truncation Issue - Current Status and Direction

## Problem Summary
wfmash is experiencing mapping output truncation where not all query sequences produce output, even though the mapping process appears to run. This manifests as missing queries in the output PAF file.

## Repository State
- **Current branch**: `fix-mapping-fundamental` (branched from main)
- **Previous attempt branch**: `fix-mapping-truncation-clean` (pushed but ineffective)
- **Other relevant branches**: 
  - `high-concurrency-shuo` - Pre-loads all queries, works but uses high memory and still has some reliability issues
  - `work-stealing-subflows` - Failed attempt with work-stealing pattern

## Key Findings

### Root Cause Analysis
The issue stems from improper task completion tracking in nested taskflow subflows:

1. **Current Structure** (in main):
```cpp
subset_flow (tf::Taskflow)
  └── processQueries_task (tf::Subflow& sf)  
       └── for each query: sf.emplace(...)
            └── query_task (tf::Subflow& query_sf)  // NESTED SUBFLOW
                 └── for each fragment: query_sf.emplace(...)
                 └── query_sf.join()  // Attempts to wait for fragments
```

2. **The Problem**:
- `executor.run(*subset_flow).wait()` only waits for top-level tasks
- Nested subflows created dynamically aren't properly tracked
- Evidence: Code has a secondary wait loop `while (executor.num_topologies() > 0)` suggesting the developers knew `wait()` wasn't sufficient
- Thread-local FAIDX readers can be destroyed while still in use

### Attempted Fixes

#### 1. Flat Taskflow Approach (`fix-mapping-truncation-clean`)
- **Changes**: Removed nested subflows, made tasks flat, used static thread-local readers
- **Result**: FAILED - Still missing queries, 10-20x performance regression
- **Why it failed**: Unknown - possibly the flat structure at query level is too granular

#### 2. Shuo's Pre-loading Approach (`high-concurrency-shuo`)
- **Changes**: Loads all query sequences into memory upfront, uses semaphore control
- **Result**: WORKS but high memory usage and still has occasional missing queries
- **Why it (mostly) works**: Eliminates concurrent file I/O, simpler task structure

#### 3. Work-Stealing Approach (`work-stealing-subflows`)
- **Changes**: Fixed worker pool with work-stealing
- **Result**: FAILED - Still has issues
- **Why it failed**: Still had nested structures and thread-local reader issues

## Current Understanding

### What We Know Works
- Pre-loading sequences (eliminates file I/O races)
- Simpler task structures (fewer nested levels)
- Synchronous index building

### What Doesn't Work
- Completely flat structure at query level (massive performance hit)
- Thread-local readers without careful lifecycle management
- Relying on `executor.run().wait()` with nested subflows

### The Nested Subflow Question
We need to understand if we're using taskflow's nested subflows correctly:
- Is `query_sf.join()` sufficient to wait for fragment tasks?
- Does `executor.run().wait()` properly track dynamically created nested subflows?
- Are we missing some taskflow pattern for proper completion tracking?

## Proposed Direction

### Option 1: Fix Nested Subflows Properly
- Research proper taskflow nested subflow patterns
- Ensure all dynamically created tasks are tracked
- Add explicit synchronization if needed
- Keep the efficient nested structure but fix completion tracking

### Option 2: Hybrid Approach
- Keep subflows for organization but simplify:
  - Process queries in batches (not all at once like Shuo)
  - Use a queue + fixed workers pattern
  - Ensure each batch completes fully before next

### Option 3: Replace Problematic Parts
- Keep taskflow for top-level orchestration
- Replace the nested query/fragment processing with simpler threading
- Use std::async or thread pool for query processing

## Code Locations

### Key Files
- `/home/erik/wfmash/src/map/include/computeMap.hpp` - Main mapping logic
- Line 329: `mapQuery()` function start
- Line 510-542: processQueries_task with nested subflows
- Line 756: `executor.run(*subset_flow).wait()` call
- Line 767-771: Secondary wait loop checking `executor.num_topologies()`

### Test Setup
```bash
# Build
cd /home/erik/wfmash
mkdir -p build && cd build
cmake .. && make -j4

# Test command (should produce 7 unique queries in output)
./build/bin/wfmash -m -t 8 cerevisiae.chrV.fa.gz cerevisiae.chrV.fa.gz > test.paf 2>&1
cut -f1 test.paf | sort -u | wc -l  # Should be 7

# Test file
cerevisiae.chrV.fa.gz - contains 7 sequences
```

## Next Steps

1. **Understand Taskflow Better**: Research proper nested subflow usage patterns in C++ Taskflow
2. **Minimal Reproducer**: Create a simple test case that demonstrates the truncation
3. **Instrumentation**: Add logging to track exactly which queries start/complete
4. **Consider Alternatives**: If taskflow's nested subflows are fundamentally problematic, consider partial replacement

## Important Notes

- Performance is critical - the flat approach's 10-20x slowdown is unacceptable
- Memory usage matters - Shuo's full pre-loading isn't ideal for large datasets  
- The issue is intermittent/timing-dependent, suggesting race conditions
- Even Shuo's "working" solution still has reliability issues

## Questions to Answer

1. Why does the secondary wait loop `while (executor.num_topologies() > 0)` exist?
2. What exactly does `query_sf.join()` wait for?
3. How does taskflow track dynamically created tasks in nested subflows?
4. Is there a taskflow pattern we're missing for this use case?

## Session State
- Branch: `fix-mapping-fundamental`
- Ready to implement a solution once we understand the proper taskflow patterns
- All analysis documents in repository root (*.md files)