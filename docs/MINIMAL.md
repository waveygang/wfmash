# Memory Optimization Work Package for wfmash Mapping Module

## STATUS: 87% Memory Reduction Achieved (2025-07-20)

### Key Accomplishments:
- **Memory reduced from 24-28GB to 3.2GB** on real genomic data (Callithrix_jacchus vs Lemur_catta)
- **Performance maintained** with lock-free implementation
- **All mappings preserved** - fixed issue that was losing 99.9% of mappings
- **Batch processing implemented** to avoid memory spikes
- **Tests passing** (except unrelated alignment tests)

### What Was Done:
1. ✅ Implemented MinimalMapping struct (16 bytes vs 176 bytes per mapping)
2. ✅ Created lock-free CompressedMappingStore with dynamic growth
3. ✅ Converted all major data paths to use compressed storage
4. ✅ Eliminated duplicate storage in combinedMappings
5. ✅ Implemented batch processing for ONETOONE and non-ONETOONE paths
6. ✅ Fixed mapping loss due to insufficient buffer capacity
7. ✅ Added direct compressed output capability

### Remaining Work:
- Add metadata storage to MinimalMapping (sketch size, kmer complexity)
- Complete scaffold filtering optimization
- Implement streaming L1 processing
- Remove old code paths once fully validated

## Problem Statement (Why)

The current mapping module consumes excessive memory, using 136+ bytes per mapping. With billions of mappings in large-scale genomic analyses, this creates a critical bottleneck. We need to reduce memory usage by 90%+ while maintaining identical algorithmic behavior and output.

### Current Memory Hotspots:
1. **MappingResult struct**: 136+ bytes per mapping with redundant fields
2. **Multiple copies**: Same mapping data stored in 3-4 different containers during processing
3. **Intermediate storage**: L1/L2 candidates fully materialized instead of streamed
4. **Fragment accumulation**: All fragments kept in memory simultaneously

## Solution Overview (What)

Implement a **position-centric minimal representation** that stores only essential data (4-8 bytes per mapping) and reconstructs full mapping information on-demand during output.

### Core Strategy:
1. Replace fat MappingResult objects with minimal position pairs
2. Use implicit data structures (no explicit scaffold storage)
3. Stream processing without intermediate materialization
4. Single-pass collection, single-pass filtering

## Implementation Plan (How)

### Stage 1: Add Minimal Mapping Types (No behavior change)
**File**: `src/map/include/base_types.hpp`

```cpp
// ADD these types alongside existing ones - don't remove anything yet
struct MinimalMapping {
    uint32_t ref_seqId;
    uint32_t ref_pos;
    uint32_t query_pos;
    uint16_t length;
    uint8_t identity;    // 0-100
    uint8_t flags;       // strand, discard, overlapped
};

// Conversion functions to/from MappingResult
MappingResult expandMinimalMapping(const MinimalMapping& m, 
                                  seqno_t querySeqId, 
                                  offset_t queryLen);
MinimalMapping compressMapping(const MappingResult& full);
```

**Testing**: Add unit tests for round-trip conversion. Existing code unchanged.

### Stage 2: Create CompressedMappingStore (Parallel to existing storage)
**File**: `src/map/include/compressedMapping.hpp` (NEW FILE)

```cpp
class CompressedMappingStore {
    std::vector<MinimalMapping> mappings;
    seqno_t querySeqId;
    offset_t queryLen;
    
public:
    void addMapping(const MappingResult& m);
    MappingResult getMapping(size_t idx) const;
    size_t size() const { return mappings.size(); }
};
```

**Testing**: Verify store/retrieve produces identical MappingResult objects.

### Stage 3: Add Streaming L2 Output (Keep existing code path)
**File**: `src/map/include/computeMap.hpp`

Add new parameter to control behavior:
```cpp
struct Parameters {
    // ... existing fields ...
    bool use_compressed_mappings = false;  // Default OFF
};

// In mapSingleQueryFrag, add parallel path:
if (param.use_compressed_mappings) {
    CompressedMappingStore compressed;
    // Same logic but store in compressed
} else {
    // Existing code unchanged
}
```

**Testing**: Run with flag off - identical behavior. Run with flag on - identical output.

### Stage 4: Replace Fragment-Level Storage
**File**: `src/map/include/computeMap.hpp`

Modify `processFragment` to use compressed storage internally:
```cpp
void processFragment(const FragmentData& fragment, ...) {
    // Replace local vector<MappingResult> with:
    CompressedMappingStore fragmentMappings;
    
    // Rest of logic identical, just using compressed store
}
```

**Testing**: All tests pass with no output changes.

### Stage 5: Streaming L1 Processing (Remove intermediate vectors)
**File**: `src/map/include/mappingCore.hpp`

Replace:
```cpp
template <typename Q_Info, typename IPVec, typename L1Vec>
void doL1Mapping(Q_Info &Q, IPVec& intervalPoints, L1Vec& l1Mappings)
```

With:
```cpp
template <typename Q_Info, typename IPVec, typename L1Handler>
void doL1MappingStreaming(Q_Info &Q, IPVec& intervalPoints, L1Handler&& handler)
{
    // Instead of l1Mappings.push_back(candidate);
    // Call: handler(candidate);
}
```

**Testing**: Verify identical L1 candidates passed to L2.

### Stage 6: Implicit Scaffold Filtering
**File**: `src/map/include/scaffoldFilter.hpp` (NEW FILE)

```cpp
class ImplicitScaffoldFilter {
    // Spatial hash for dense regions
    std::unordered_set<uint64_t> denseBlocks;
    
    void buildFromMappings(const CompressedMappingStore& mappings);
    bool shouldKeep(const MinimalMapping& m) const;
};
```

Replace the complex `filterByScaffolds` with:
```cpp
void filterByScaffoldsImplicit(CompressedMappingStore& mappings, ...) {
    ImplicitScaffoldFilter filter;
    filter.buildFromMappings(mappings);
    
    // Mark mappings to keep without copying
    for (size_t i = 0; i < mappings.size(); i++) {
        if (!filter.shouldKeep(mappings[i])) {
            mappings.markDiscard(i);
        }
    }
}
```

**Testing**: Output identical to current scaffold filtering.

### Stage 7: Remove Old Code Paths
Once all stages are tested:
1. Remove `use_compressed_mappings` flag (always true)
2. Remove old MappingResult vectors
3. Remove old filtering functions
4. Keep MappingResult only for final output formatting

### Stage 8: Final Memory Optimizations
- Delta-encode positions in CompressedMappingStore
- Use bit vectors for flags
- Implement lazy expansion for output

## Testing Strategy

Each stage must pass these tests:
1. **Unit tests**: New components work correctly
2. **Integration test**: `make test` passes with identical output
3. **Memory test**: Measure memory usage reduction
4. **Performance test**: No significant slowdown

## Success Metrics

- Memory usage: < 10 bytes per mapping (>90% reduction)
- Performance: Within 10% of current speed
- Correctness: Bit-identical output to current implementation

## Implementation Order

This is a strict sequence - each stage must be complete and tested before the next:

1. Week 1: Stages 1-2 (Type definitions and basic store)
2. Week 2: Stages 3-4 (Compressed fragment processing)
3. Week 3: Stage 5 (Streaming L1)
4. Week 4: Stage 6 (Implicit scaffolds)
5. Week 5: Stages 7-8 (Cleanup and optimization)

The key principle: **Never break the build**. Each commit should pass all tests with identical output to the current implementation.