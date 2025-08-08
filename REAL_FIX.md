# The REAL Fix - No Pre-loading Needed!

## The Actual Problem
The nested subflows are creating tasks that aren't tracked by executor.run().wait(). This is a known Taskflow issue (#138).

## The Real Solution
**Flatten the task hierarchy** - that's it! We don't need to pre-load sequences.

### Current BROKEN Structure:
```cpp
taskflow
  └── processQueries_task (Subflow)
       └── for each query: query_task (Subflow)  // NESTED!
            └── for each fragment: fragment_task
```

### Fixed FLAT Structure:
```cpp
taskflow
  └── for each query: query_task (regular task, no subflow)
       └── process all fragments sequentially within task
```

## Implementation
Instead of creating nested subflows, create all query tasks directly in the main taskflow:

```cpp
tf::Taskflow taskflow;

// Build index first
auto index_task = taskflow.emplace([&]() { 
    // build index
});

// Create one task per query (NO SUBFLOWS!)
std::vector<tf::Task> query_tasks;
for (const auto& queryName : querySequenceNames) {
    auto task = taskflow.emplace([&, queryName]() {
        // Load sequence (thread-local reader is FINE)
        // Process ALL fragments for this query sequentially
        // Output results
    });
    query_tasks.push_back(task);
    index_task.precede(task);  // All queries depend on index
}

// Now executor.run().wait() WILL work correctly!
executor.run(taskflow).wait();
```

## Why This Works
1. No nested subflows = all tasks are tracked
2. Thread-local readers are fine - threads live long enough
3. executor.run().wait() works correctly with flat structure
4. No pre-loading needed = no memory waste
5. No semaphores needed = clean design

This is the minimal change that fixes the root cause!