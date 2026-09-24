# CLAUDE.md - Memory Optimization Implementation Guide

## Project Overview
This document tracks the implementation of memory-efficient mapping structures in wfmash to replace the large MappingResult struct (240+ bytes) with a compact MinimalMapping struct (32 bytes).

## Stable Baseline Revision
**Git SHA**: 8de04614196f34ca4189ccf61814e2d82ddaae58
**Date**: 2025-07-22
**Description**: Last stable revision before memory optimization changes. Use this for baseline generation and comparison.

## Implementation Strategy: "Rip and Replace"
The strategy is to perform a direct swap of the primary mapping object, then use compiler errors as a precise to-do list for adapting the surrounding code.

### Core Principles
1. **Never expand compressed mappings unless absolutely necessary** - Only expand for complex algorithms or final output
2. **Use expand/compress wrappers** - Don't rewrite complex algorithms; wrap them with conversion functions
3. **Trust the compiler** - Compilation errors after the initial change are your roadmap

## Key Data Structures

### MinimalMapping (32 bytes)
```cpp
struct MinimalMapping {
    uint32_t ref_seqId;         // 4 bytes
    uint32_t ref_pos;           // 4 bytes  
    uint32_t query_pos;         // 4 bytes
    uint32_t length;            // 4 bytes
    uint32_t n_merged;          // 4 bytes
    uint32_t conservedSketches; // 4 bytes
    uint16_t identity;          // 2 bytes (scaled 0-10000)
    uint8_t  flags;             // 1 byte (strand, discard, overlapped)
    uint8_t  kmerComplexity;    // 1 byte (scaled 0-100)
};
```

### Conversion Functions
- `compressMapping(const FullMapping&) -> MinimalMapping`
- `expandMapping(const MinimalMapping&, seqno_t querySeqId, offset_t queryLen) -> FullMapping`

## Pipeline Stages & Testing

### 1. Raw Mappings (L2 Mapping)
- **Function**: `doL2Mapping` in `computeMap.hpp`
- **Test**: `wfmash -m -t 8 -M -f none -S 0 cerevisiae.chrV.fa.gz > baseline.raw.paf`
- **Purpose**: Tests initial mapping generation and compression

### 2. Merged Mappings (Union-Find)
- **Function**: `mergeMappingsInRange` in `mappingFilter.hpp`
- **Test**: `wfmash -m -t 8 -f none -S 0 cerevisiae.chrV.fa.gz > baseline.merged.paf`
- **Purpose**: Tests the complex union-find merge algorithm
- **Note**: Requires full mapping data; must use expand-execute-compress pattern

### 3. Scaffold Filtering
- **Function**: `filterByScaffolds` in `mappingFilter.hpp`
- **Test**: Part of full pipeline test
- **Purpose**: Identifies large syntenic blocks
- **Note**: Complex two-tier filtering; requires full mapping data

### 4. Full Pipeline
- **Test**: `wfmash -m -t 8 cerevisiae.chrV.fa.gz > baseline.full.paf`
- **Purpose**: Final acceptance test of all components

## Implementation Checklist

### Phase 1: Core Structure Changes
1. Create `src/map/include/compact_mapping.hpp` with MinimalMapping and conversion functions
2. Modify `base_types.hpp`: rename MappingResult to FullMapping, redefine MappingResultsVector_t
3. Update `computeMap.hpp`: compress mappings at creation
4. Update `mappingOutput.hpp`: expand mappings for reporting

### Phase 2: Fix Compilation Errors
For each error where MinimalMapping lacks required fields:
1. Create temporary FullMappingResultsVector_t at function entry
2. Expand all MinimalMappings to FullMappings
3. Run existing algorithm unchanged on full mappings
4. Compress results back to MinimalMappings

### Key Functions Requiring Wrappers
- `mergeMappingsInRange` - Complex union-find algorithm
- `filterByScaffolds` - Two-tier syntenic filtering
- `filterByGroup` - Plane-sweep filtering
- Any sorting functions using fields not in MinimalMapping

## Testing Protocol

### Before Any Changes
Generate all baselines:
```bash
# Raw mappings (no merging, no filtering)
wfmash -m -t 8 -M -f none -S 0 cerevisiae.chrV.fa.gz > baseline.raw.paf

# Merged mappings (merging only, no filtering)
wfmash -m -t 8 -f none -S 0 cerevisiae.chrV.fa.gz > baseline.merged.paf

# Full pipeline (all processing)
wfmash -m -t 8 cerevisiae.chrV.fa.gz > baseline.full.paf
```

### After Each Change
1. Generate current outputs with same commands
2. Check in order: raw → merged → full
3. First failure indicates which stage broke

## Memory Impact
- Original MappingResult: 240+ bytes per mapping
- MinimalMapping: 32 bytes per mapping
- Expected reduction: ~87% memory usage for mapping storage

## Critical Warnings
1. Do NOT modify the union-find merge algorithm - it's too complex
2. Do NOT modify the scaffold filtering logic - it uses two-tier filtering
3. Always pass querySeqId and queryLen when expanding mappings
4. The ScaffoldIndex mentioned in MAP_COMPACT.md is a separate optimization - ignore for now

## Status Updates
- 2025-07-22: Initial documentation created, stable baseline recorded