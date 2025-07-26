# Scaffold Logging Simplification

## Summary
Simplified the scaffold progress tracking by removing complex thread spawning, atomic counters, and separate progress meters. Replaced with a simple logging message.

## Problems with Previous Implementation
1. **Deadlock risk**: Complex thread synchronization with atomic counters
2. **Overhead**: Separate monitoring thread checking progress every 10ms
3. **Complexity**: Multiple shared pointers, atomic operations, and memory ordering
4. **Unnecessary**: Scaffold filtering is typically fast enough not to need detailed progress

## Changes Made

### 1. In `computeMap.hpp`
- **Removed**: Thread spawning for scaffold progress monitoring
- **Removed**: Atomic counters (`scaffold_total_work`, `scaffold_completed_work`)
- **Removed**: Separate progress meter for scaffolding
- **Added**: Simple log message: `[wfmash::mashmap] Scaffolding mappings...`

### 2. In `mappingFilter.hpp`
- **Modified**: `filterByScaffolds` function signature to make progress tracking parameters optional with default nullptr
- **Removed**: All `fetch_add` calls to update scaffold work counters
- **Removed**: Progress tracking from the parallel distance computation

## Benefits
1. **No deadlock risk**: No threads or synchronization needed
2. **Simpler code**: Easier to maintain and debug
3. **Same functionality**: Scaffold filtering still works exactly the same
4. **Better performance**: No overhead from progress tracking

## Testing
Confirmed that:
- Program compiles successfully
- Scaffold filtering still works (`--scaffold-out` produces output)
- No deadlocks or hangs
- Simple message "Scaffolding mappings..." is displayed

## Example Output
```
[wfmash::mashmap] mapping  [100.0% complete, 4029392/4029392 units, 0s elapsed]
[wfmash::mashmap] mapping  completed.
[wfmash::mashmap] Scaffolding mappings...
[wfmash::mashmap] Mapped query in 0.4s, results saved to: /dev/stdout
```