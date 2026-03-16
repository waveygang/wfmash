# Taskflow Nested Subflow Fix Analysis

## Problem Identification

After thorough research and code analysis, I've identified the root cause of the mapping truncation issue:

### Core Issues

1. **Improper Nested Subflow Synchronization**
   - Current structure has nested subflows: `subset_flow` → `processQueries_task` (sf) → `query_task` (query_sf)
   - `executor.run(*subset_flow).wait()` only waits for top-level tasks to be scheduled, not for nested subflow completion
   - The secondary wait loop `while (executor.num_topologies() > 0)` suggests developers knew `.wait()` wasn't sufficient

2. **Thread-Local Reader Lifetime Problem**
   - Thread-local FAIDX readers can be destroyed while tasks are still using them
   - Thread-local destructors run when thread exits, not when all tasks complete
   - This causes use-after-free issues when tasks access destroyed readers

3. **Lambda Capture Hazards**
   - Nested lambdas capture references that may become invalid
   - Parent task completion doesn't guarantee nested task completion
   - Race conditions between task completion and reference validity

## Research Findings

From Taskflow documentation:
- `subflow.join()` only waits for tasks within that specific subflow
- `executor.run().wait()` doesn't properly track dynamically created nested subflows
- `wait_for_all()` is the correct method to wait for ALL submitted tasks
- Nested subflows should be avoided when possible for better synchronization

## Proposed Solution

### Approach 1: Flatten the Subflow Structure
Remove one level of nesting by processing queries directly in the first subflow without creating nested query subflows.

### Approach 2: Use Proper Synchronization
Replace `executor.run().wait()` with proper synchronization that accounts for all dynamic tasks.

### Approach 3: Fix Reader Lifetime Management
Move readers out of thread-local storage to ensure they remain valid for task lifetime.

## Implementation Strategy

1. **Immediate Fix**: Replace the nested subflow pattern with a flatter structure
2. **Synchronization Fix**: Use `wait_for_all()` or maintain topology count properly
3. **Reader Fix**: Create readers with proper lifetime management
4. **Testing**: Verify all queries are processed without truncation