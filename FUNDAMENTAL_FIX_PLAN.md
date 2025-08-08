# Fundamental Fix Plan: True Flat Taskflow

## The Core Problem

The current structure has queries being processed INSIDE a subflow task:
```
subset_flow
  └── processQueries_task (creates subflow)
       └── for each query: creates tasks in subflow
            └── query_task (creates another subflow!)  
                 └── fragment tasks
```

`executor.run().wait()` only waits for the top-level tasks, not the dynamically created subtasks.

## The Solution: True Flat Structure

```cpp
// Build index first (synchronously)
buildIndex();

// Create taskflow with one task per query at TOP LEVEL
tf::Taskflow taskflow;
for (const auto& queryName : querySequenceNames) {
    taskflow.emplace([queryName]() {
        // Load sequence
        // Process ALL fragments inline (no subtasks)
        // Output results
    });
}

// This will actually wait for all queries
executor.run(taskflow).wait();

// Clean up
delete refSketch;
```

## Key Changes

1. **NO processQueries_task** - eliminate this wrapper entirely
2. **Query tasks at TOP LEVEL** - each query is a direct child of the main taskflow
3. **NO nested subflows** - zero use of tf::Subflow
4. **Synchronous index building** - build index before creating taskflow
5. **Inline fragment processing** - process fragments in a loop, not as tasks

## Why This Will Work

- `executor.run(taskflow).wait()` will track ALL query tasks because they're at the top level
- No nested subflows means no ambiguity about what needs to complete
- Every query MUST complete before `wait()` returns
- The thread pool (8 threads) naturally limits concurrency

## Implementation Steps

1. Move index building out of taskflow (make it synchronous)
2. Move query iteration to top level (not inside processQueries_task)
3. Create query tasks directly in main taskflow
4. Remove ALL tf::Subflow usage
5. Process fragments inline within each query task
6. Remove the `executor.num_topologies() > 0` wait loop (won't be needed)

## Testing

After implementation, verify:
1. All queries produce output
2. No truncation with multiple queries
3. Performance is reasonable (should be faster without subflow overhead)
4. Works with both single and multiple target subsets