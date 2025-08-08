# Clean Solution for Mapping Truncation Issue

## The Real Problem
After extensive research, the core issue is clear:
1. **Nested subflows in Taskflow are fundamentally problematic** - executor.run().wait() doesn't properly track dynamically created nested tasks
2. **Thread-local resources with dynamic task creation is dangerous** - lifetime issues
3. **The semaphore in Shuo's branch is a workaround, not a fix**

## What Actually Works (from Shuo's branch)
1. **Pre-loading all sequences** - Eliminates concurrent file I/O races
2. **Simpler task structure** - Less nesting = fewer synchronization issues
3. **Sequential I/O followed by parallel processing** - Clean separation of concerns

## Clean Solution Design

### Architecture
```
1. Pre-load Phase (Sequential)
   - Load all query sequences into memory once
   - Use single reader, no concurrency issues
   - Store sequences in simple vector

2. Processing Phase (Parallel)
   - Create simple taskflow with one task per query
   - Each query task processes its fragments sequentially
   - No nested subflows, just simple parallel tasks
   
3. Output Phase
   - Simple mutex-protected output
   - No complex synchronization needed
```

### Key Principles
1. **NO nested subflows** - They're the root of all evil here
2. **NO thread-local storage** - Too many lifetime issues  
3. **NO semaphores** - If you need them, your design is wrong
4. **Simple is better** - One level of parallelism is enough

### Implementation Steps
1. Pre-load all sequences into a vector of structs
2. Create a simple taskflow with one task per query
3. Each task processes all fragments for that query
4. Use executor.run().wait() which WILL work with flat structure
5. No need for wait_for_all() or num_topologies() loops

### Why This Will Work
- Flat task structure = proper synchronization
- Pre-loaded sequences = no file I/O races
- Simple parallelism = predictable behavior
- This is essentially what Shuo's branch does, minus the unnecessary complexity

### Memory Impact
- For cerevisiae.chrV.fa.gz: ~3.8MB total sequence data
- For human genome: ~3GB (acceptable on modern systems)
- Can add batching later if needed for huge datasets