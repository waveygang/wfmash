# Scaffold Filter Update - Implementation Summary

## Overview
The scaffold filtering logic in `mappingFilter.hpp` has been updated to better remove "random off-diagonal shorter weaker scaffold chains" by applying the plane sweep algorithm to the merged chains before selecting anchors.

## Key Changes

### Previous Logic
1. Merge mappings into chains
2. Filter chains by minimum length
3. Collect ALL mappings within chain bounds as anchors
4. Apply plane sweep filter to anchors
5. Use filtered anchors for distance-based filtering

### New Logic
1. Merge mappings into chains
2. Filter chains by minimum length
3. **NEW: Apply plane sweep filter to merged chains** (Step 4)
4. Collect mappings within filtered chain bounds as anchors (Step 5)
5. Use anchors for distance-based filtering

## Implementation Details

The critical change is in Step 4, where we now filter the scaffold chains themselves:

```cpp
// Step 4: Apply plane sweep filter to the merged chains to remove off-diagonal/weaker scaffold chains
if (!mergedChains.empty() && (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE)) {
    MappingResultsVector_t filteredChains;
    filterByGroup(mergedChains, filteredChains, param.numMappingsForSegment - 1, 
                 false, idManager, param, progress);
    mergedChains = std::move(filteredChains);
}
```

This ensures that:
1. Weaker scaffold chains that would create off-diagonal noise are removed
2. Only the best-scoring chains survive to contribute anchors
3. The final anchor set is cleaner and more focused on true syntenic regions

## Testing

The implementation has been tested with cerevisiae chrV data:
- Default parameters (S=5): 990 mappings retained, 666 scaffold chains
- Higher threshold (S=10): 876 mappings retained, 568 scaffold chains

The scaffold output file (`--scaffold-out`) now contains the filtered merged chains, showing which scaffolds survived the plane sweep filter.

## Benefits

1. **Better noise reduction**: Off-diagonal chains are removed at the chain level
2. **Cleaner anchor sets**: Only mappings from high-quality chains become anchors
3. **More focused scaffolds**: The plane sweep ensures competing/overlapping chains are resolved

## Files Modified
- `src/map/include/mappingFilter.hpp`: Updated `filterByScaffolds()` function