# Union-Find Mapping Merge Algorithm
## Comprehensive Technical Documentation

### Abstract

This document provides a complete technical specification of the sophisticated union-find based mapping merge algorithm used in wfmash v0.23.0-40-g408cbbed. The algorithm employs a combination of geometric distance scoring, union-find data structures, and chain fragmentation to merge collinear genomic mappings while preserving accuracy and maintaining computational efficiency.

### 1. Overview

The `mergeMappingsInRange` algorithm is designed to merge collinear sequence mappings that likely represent the same underlying genomic alignment split across multiple segments. The core insight is that genuine split mappings exhibit geometric consistency in both query and reference coordinate spaces, with bounded distances between consecutive segments.

**Key Properties:**
- **Input**: Vector of `MappingResult` objects representing individual mapping segments
- **Output**: Vector of maximally merged `MappingResult` objects with improved contiguity
- **Complexity**: O(n log n + n·α(n)) where α is the inverse Ackermann function
- **Geometric Model**: 2D coordinate space optimization with strand-aware distance metrics

### 2. Data Structures

#### 2.1 MappingResult Structure
```cpp
struct MappingResult {
    seqno_t refSeqId, querySeqId;           // Sequence identifiers
    offset_t refStartPos, refEndPos;        // Reference coordinates  
    offset_t queryStartPos, queryEndPos;    // Query coordinates
    strnd strand;                           // Alignment orientation
    
    // Chain tracking fields
    int64_t splitMappingId;                 // Unique mapping identifier
    double chainPairScore;                  // Distance score to best pair
    int64_t chainPairId;                    // ID of best pairing candidate
    int discard;                            // Discard flag (0/1)
    
    // Quality metrics
    float nucIdentity;                      // Nucleotide identity percentage
    int conservedSketches, sketchSize;      // k-mer sketch statistics
    float kmerComplexity;                   // Sequence complexity measure
    int n_merged;                           // Number of segments merged
};
```

#### 2.2 DisjointSets (Union-Find) Data Structure
The algorithm employs a lock-free union-find data structure with path compression and union by rank:

```cpp
class DisjointSets {
    using Aint = __uint128_t;  // 128-bit integer for parent+rank storage
    static const Aint parentMask = ~0ULL;        // Lower 64 bits = parent ID
    static const Aint rankMask = parentMask << 64; // Upper 64 bits = rank
    
public:
    uint64_t find(uint64_t id);     // Find root with path compression
    uint64_t unite(uint64_t a, uint64_t b);  // Union by rank
    bool same(uint64_t a, uint64_t b);       // Equivalence test
};
```

**Storage Optimization**: Each 128-bit entry encodes both parent pointer (64 bits) and rank (64 bits), enabling efficient atomic operations.

### 3. Algorithm Structure

The algorithm consists of five major phases:

1. **Preprocessing & Initialization**
2. **Geometric Chain Discovery** 
3. **Union-Find Merging**
4. **Maximal Fragment Construction**
5. **Chain Post-Processing**

### 4. Phase 1: Preprocessing & Initialization

#### 4.1 Input Validation
```cpp
if (!param.split || readMappings.size() < 2) return readMappings;
```
Early termination for trivial cases or when splitting is disabled.

#### 4.2 Lexicographic Sorting
```cpp
std::sort(readMappings.begin(), readMappings.end(),
    [](const MappingResult &a, const MappingResult &b) {
        return std::tie(a.refSeqId, a.strand, a.queryStartPos, a.refStartPos)
             < std::tie(b.refSeqId, b.strand, b.queryStartPos, b.refStartPos);
    });
```

**Sorting Key Priority:**
1. `refSeqId` - Group mappings by reference sequence
2. `strand` - Separate forward/reverse orientations  
3. `queryStartPos` - Order by query coordinate
4. `refStartPos` - Break ties with reference coordinate

#### 4.3 Unique Identifier Assignment
```cpp
for (auto it = readMappings.begin(); it != readMappings.end(); ++it) {
    it->splitMappingId = std::distance(readMappings.begin(), it);
    it->discard = 0;
    it->chainPairScore = std::numeric_limits<double>::max();
    it->chainPairId = std::numeric_limits<int64_t>::min();
}
```

**Field Initialization:**
- `splitMappingId`: Sequential index for union-find operations
- `discard`: Reset discard flag 
- `chainPairScore`: Initialize to infinity (no pair found)
- `chainPairId`: Initialize to sentinel value

### 5. Phase 2: Geometric Chain Discovery

#### 5.1 Group Processing
The algorithm processes mappings in groups defined by `(refSeqId, strand)` tuples:

```cpp
auto group_begin = readMappings.begin();
while (group_begin != readMappings.end()) {
    auto group_end = std::find_if_not(group_begin, readMappings.end(),
        [refSeqId = group_begin->refSeqId, strand = group_begin->strand](const MappingResult &m) {
            return m.refSeqId == refSeqId && m.strand == strand;
        });
    // Process group...
}
```

#### 5.2 Chain Pair Discovery Algorithm

For each mapping **M₁** in a group, the algorithm searches for the optimal subsequent mapping **M₂**:

##### 5.2.1 Search Space Pruning
```cpp
auto end_it2 = std::upper_bound(it + 1, group_end, it->queryEndPos + max_dist,
    [](offset_t val, const MappingResult &m) {
        return val < m.queryStartPos;
    });
```

**Binary Search Optimization**: Only consider mappings M₂ where:
```
M₂.queryStartPos ≤ M₁.queryEndPos + max_dist
```

This leverages the lexicographic sort order to prune the search space from O(n) to O(log n).

##### 5.2.2 Geometric Distance Scoring

For each candidate pair (M₁, M₂), the algorithm computes:

**Query Space Distance:**
```cpp
int64_t query_dist = it2->queryStartPos - it->queryEndPos;
```

**Reference Space Distance (Strand-Aware):**
```cpp
int64_t ref_dist = (it->strand == strnd::FWD) ? 
                   it2->refStartPos - it->refEndPos : 
                   it->refStartPos - it2->refEndPos;
```

**Geometric Constraints:**
1. **Query Gap Constraint**: `query_dist ≤ max_dist`
2. **Reference Gap Constraint**: `-segLength/5 ≤ ref_dist ≤ max_dist`
3. **Non-Overlapping**: `query_dist > 0` (implicit from `queryStartPos` comparison)

The negative bound (`-segLength/5`) allows for small overlaps, accounting for alignment uncertainty.

##### 5.2.3 Euclidean Distance Scoring
```cpp
double dist_sq = static_cast<double>(query_dist) * query_dist + 
                 static_cast<double>(ref_dist) * ref_dist;
double max_dist_sq = static_cast<double>(max_dist) * max_dist;

if (dist_sq < max_dist_sq && dist_sq < best_score && dist_sq < it2->chainPairScore) {
    best_it2 = it2;
    best_score = dist_sq;
}
```

**Scoring Function**: Euclidean distance in (query_gap, ref_gap) space
- **Geometric Interpretation**: Measures deviation from perfect collinearity
- **Optimality**: Each mapping can only be the best successor for one predecessor
- **Competition**: Multiple predecessors compete for the same successor

##### 5.2.4 Best Pair Assignment
```cpp
if (best_it2 != group_end) {
    best_it2->chainPairScore = best_score;
    best_it2->chainPairId = it->splitMappingId;
}
```

**Pairing Semantics:**
- Each mapping stores pointer to its **best predecessor**
- Forms a reverse-directed forest structure
- Multiple mappings may point to the same predecessor

### 6. Phase 3: Union-Find Merging

#### 6.1 Union-Find Initialization
```cpp
std::vector<dsets::DisjointSets::Aint> ufv(readMappings.size());
auto disjoint_sets = dsets::DisjointSets(ufv.data(), ufv.size());
```

**Initial State**: Each mapping forms its own singleton set with `parent[i] = i`.

#### 6.2 Chain Union Operations
```cpp
for (auto it = group_begin; it != group_end; ++it) {
    if (it->chainPairScore != std::numeric_limits<double>::max()) {
        disjoint_sets.unite(it->splitMappingId, it->chainPairId);
    }
}
```

**Union Strategy:**
- Union mapping with its best predecessor
- Creates connected components representing chains
- **Transitivity**: If A→B and B→C, then A and C are in the same chain

#### 6.3 Representative Assignment
```cpp
for (auto &mapping : readMappings) {
    mapping.splitMappingId = disjoint_sets.find(mapping.splitMappingId);
}
```

**Path Compression**: All mappings in a chain get the same representative ID.

### 7. Phase 4: Maximal Fragment Construction

#### 7.1 Chain-Based Grouping
```cpp
std::sort(readMappings.begin(), readMappings.end(),
    [](const MappingResult &a, const MappingResult &b) {
        return std::tie(a.splitMappingId, a.queryStartPos, a.refStartPos)
             < std::tie(b.splitMappingId, b.queryStartPos, b.refStartPos);
    });
```

**Secondary Sort**: Within each chain, order by genomic coordinates.

#### 7.2 Fragment Generation Algorithm

For each chain, the algorithm creates maximal fragments respecting length constraints:

```cpp
// Process chain into fragments
auto fragment_start = it;
while (fragment_start != it_end) {
    MappingResult mergedMapping = *fragment_start;
    auto fragment_end = fragment_start;
    auto next = std::next(fragment_end);
   
    offset_t current_length = fragment_start->queryEndPos - fragment_start->queryStartPos;
    
    while (next != it_end) {
        offset_t potential_length = next->queryEndPos - fragment_start->queryStartPos;
        if (potential_length > param.max_mapping_length) break;
        fragment_end = next;
        current_length = potential_length;
        next = std::next(next);
    }
    // Merge fragment...
}
```

**Greedy Extension**: Extend each fragment maximally while respecting `max_mapping_length`.

#### 7.3 Fragment Merging Mathematics

**Coordinate Merging:**
```cpp
mergedMapping.queryStartPos = fragment_start->queryStartPos;
mergedMapping.queryEndPos = fragment_end->queryEndPos;

if (mergedMapping.strand == strnd::FWD) {
    mergedMapping.refStartPos = fragment_start->refStartPos;
    mergedMapping.refEndPos = fragment_end->refEndPos;
} else {
    mergedMapping.refStartPos = fragment_end->refStartPos;
    mergedMapping.refEndPos = fragment_start->refEndPos;
}
```

**Strand-Aware Reference Merging**: Reverse strand mappings require coordinate reversal.

**Quality Metrics Aggregation:**
```cpp
double totalNucIdentity = 0.0;
double totalKmerComplexity = 0.0;
int totalConservedSketches = 0;
int totalSketchSize = 0;
int fragment_size = 0;

for (auto subIt = fragment_start; subIt != std::next(fragment_end); ++subIt) {
    totalNucIdentity += subIt->nucIdentity;
    totalKmerComplexity += subIt->kmerComplexity;
    totalConservedSketches += subIt->conservedSketches;
    totalSketchSize += subIt->sketchSize;
    fragment_size++;
}

mergedMapping.n_merged = fragment_size;
mergedMapping.nucIdentity = totalNucIdentity / fragment_size;
mergedMapping.kmerComplexity = totalKmerComplexity / fragment_size;
mergedMapping.conservedSketches = totalConservedSketches;
mergedMapping.sketchSize = totalSketchSize;
```

**Aggregation Rules:**
- **Identity & Complexity**: Arithmetic mean
- **Sketch Counts**: Sum (additive evidence)
- **Merge Count**: Fragment size

### 8. Phase 5: Chain Post-Processing

#### 8.1 Chain Split Processing
```cpp
processChainWithSplits(it, it_end, param);
```

The `processChainWithSplits` function applies several refinement operations:

##### 8.1.1 Gap Closure Analysis
```cpp
std::vector<bool> is_cuttable(std::distance(begin, end), true);

for (auto it = std::next(begin); it != end; ++it) {
    auto prev = std::prev(it);
    if (it->queryStartPos - prev->queryEndPos > param.segLength / 5 ||
        it->refStartPos - prev->refEndPos > param.segLength / 5) {
        is_cuttable[std::distance(begin, prev)] = false;
        is_cuttable[std::distance(begin, it)] = false;
    }
}
```

**Cuttability Rules**: Positions with large gaps (> segLength/5) cannot be used for fragmentation cuts.

##### 8.1.2 Consecutive Mapping Adjustment
```cpp
void adjustConsecutiveMappings(begin_mapping, end_mapping, threshold) {
    for (auto it = std::next(begin_mapping); it != end_mapping; ++it) {
        auto& prev = *std::prev(it);
        auto& curr = *it;
        
        int query_gap = curr.queryStartPos - prev.queryEndPos;
        int ref_gap = curr.refStartPos - prev.refEndPos;
        
        if (query_gap > 0 && ref_gap > 0 && query_gap <= threshold && ref_gap <= threshold) {
            int query_mid = (prev.queryEndPos + curr.queryStartPos) / 2;
            int ref_mid = (prev.refEndPos + curr.refStartPos) / 2;
            
            prev.queryEndPos = query_mid;
            prev.refEndPos = ref_mid;
            curr.queryStartPos = query_mid;
            curr.refStartPos = ref_mid;
        }
    }
}
```

**Gap Closure**: Small gaps (≤ threshold) are closed by extending boundaries to the midpoint.

#### 8.2 Fragment Length Enforcement
```cpp
auto fragment_start = begin;
offset_t accumulate_length = 0;

for (auto it = begin; it != end; ++it) {
    accumulate_length += it->queryEndPos - it->queryStartPos;
    
    if (accumulate_length >= param.max_mapping_length && is_cuttable[std::distance(begin, it)]) {
        processMappingFragment(fragment_start, std::next(it));
        fragment_start = std::next(it);
        accumulate_length = 0;
    }
}
```

**Length Constraint Enforcement**: Chains exceeding `max_mapping_length` are cut at cuttable positions.

#### 8.3 Discard Flag Management
```cpp
readMappings.erase(
    std::remove_if(readMappings.begin(), readMappings.end(), 
                   [](const MappingResult& e) { return e.discard == 1; }),
    readMappings.end()
);
```

**Cleanup**: Mappings marked for discard during post-processing are removed.

### 9. Complexity Analysis

#### 9.1 Time Complexity

**Phase 1 (Preprocessing)**: O(n log n)
- Lexicographic sorting dominates

**Phase 2 (Chain Discovery)**: O(n log n)  
- For each mapping: O(log n) binary search + O(k) distance calculations
- Expected k << n due to geometric constraints

**Phase 3 (Union-Find)**: O(n α(n))
- n union operations with path compression
- α(n) is inverse Ackermann function, effectively constant

**Phase 4 (Fragment Construction)**: O(n log n)
- Secondary sort + linear fragment processing

**Phase 5 (Post-Processing)**: O(n)
- Linear chain processing operations

**Total**: O(n log n)

#### 9.2 Space Complexity

**Union-Find Storage**: O(n) × 16 bytes = O(16n)
**Chain Pair Storage**: O(1) per mapping (embedded fields)
**Auxiliary Structures**: O(n)

**Total**: O(n)

### 10. Geometric Properties & Correctness

#### 10.1 Collinearity Preservation
The algorithm preserves geometric collinearity through:
1. **Euclidean Distance Minimization**: Favors pairs with minimal coordinate deviation
2. **Strand Consistency**: Separate processing for forward/reverse orientations  
3. **Gap Constraints**: Bounded allowable gaps prevent spurious connections

#### 10.2 Optimality Properties
- **Local Optimality**: Each mapping chooses the best available predecessor
- **Global Structure**: Union-find ensures transitivity of chaining relationships
- **Length Constraints**: Fragments respect maximum mapping length bounds

### 11. Parameter Dependencies

#### 11.1 Critical Parameters
- `max_dist`: Maximum allowable gap between chain segments
- `param.segLength`: Segment length (affects overlap tolerance)
- `param.max_mapping_length`: Maximum fragment length
- `param.split`: Enable/disable splitting functionality

#### 11.2 Parameter Interactions
- **Gap Tolerance**: Reference gap tolerance = `max_dist` ± `segLength/5`
- **Cuttability Threshold**: `segLength/5` determines non-cuttable positions
- **Fragment Length**: Trade-off between contiguity and length constraints

### 12. Implementation Notes

#### 12.1 Numerical Stability
- **Double Precision**: Distance calculations use `double` to avoid overflow
- **Integer Overflow**: Gap distances use `int64_t` for large genomes
- **Score Initialization**: Infinity values prevent uninitialized comparisons

#### 12.2 Memory Management
- **In-Place Operations**: Most operations modify existing vectors
- **Temporary Allocation**: Union-find requires O(n) temporary storage
- **Fragment Generation**: Creates new merged mappings, original chains remain

### 13. Conclusion

The union-find mapping merge algorithm represents a sophisticated approach to the genomic sequence alignment merging problem. By combining geometric optimization with efficient data structures, it achieves optimal time complexity while preserving alignment accuracy. The algorithm's success depends on careful parameter tuning and geometric constraint enforcement to distinguish genuine split mappings from spurious alignments.

The algorithm's robustness comes from its multi-phase approach: geometric discovery identifies candidate chains, union-find ensures global consistency, and post-processing enforces biological constraints. This design makes it suitable for large-scale genomic alignment applications where both accuracy and performance are critical.