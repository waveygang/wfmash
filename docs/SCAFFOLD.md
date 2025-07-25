# Scaffold Filtering in wfmash

## Overview

Scaffold filtering is a two-level hierarchical approach used by wfmash to identify regions of synteny between sequences, particularly useful for divergent species comparisons. It solves the problem of distinguishing true homologous regions from spurious matches when comparing distantly related genomes.

## The Two-Level Process

### Level 1: Local Chain Merging
- **Input**: Individual segment mappings (default 1kb segments)
- **Gap parameter**: `chain_gap` (default 2kb, set with `-c`)
- **Process**: Merge mappings that are within `chain_gap` distance in both query and reference
- **Output**: Local chains representing continuous alignable regions

### Level 2: Scaffold Construction
- **Input**: Local chains from Level 1
- **Gap parameter**: `scaffold_gap` (default 100kb, first value in `-S`)
- **Process**: Merge chains with much larger gaps to identify syntenic blocks
- **Output**: Large-scale scaffold regions representing areas of conserved synteny

### Filtering Step
- **Parameters**: 
  - `scaffold_min_length` (default 5kb, second value in `-S`)
  - `scaffold_max_deviation` (default 100kb, third value in `-S`)
- **Process**: 
  1. Keep all chains that are part of scaffolds ≥ `scaffold_min_length`
  2. Keep chains within `scaffold_max_deviation` distance of any scaffold
  3. Discard all other chains

## Example: Distant Species Comparison

When comparing Lemur (Lemur_catta) and Marmoset (Callithrix_jacchus):

1. **Initial mapping**: ~37,000 individual 1kb segment matches scattered across genomes
2. **Local chaining**: Most remain as 1kb chains due to scattered matches
3. **Scaffold creation**: ~11 large syntenic blocks identified by connecting distant chains
4. **Final filtering**: Only chains near these 11 scaffolds are kept, reducing output to biologically meaningful regions

## Command Line Usage

```bash
# Default scaffold filtering
wfmash -m query.fa reference.fa  # Uses -S 100k,5k,100k

# Disable scaffold filtering
wfmash -m -S 0 query.fa reference.fa

# Custom scaffold parameters
wfmash -m -S 200k,10k,150k query.fa reference.fa
# scaffold_gap=200kb, min_length=10kb, max_deviation=150kb
```

## When Scaffold Filtering Helps

- **Distant species**: Filters out noise from spurious matches
- **Repetitive regions**: Removes scattered repeat matches
- **Structural variation**: Identifies large-scale syntenic blocks despite local rearrangements

## When to Disable Scaffold Filtering

- **Closely related genomes**: When most mappings form continuous chains
- **Small-scale analysis**: When interested in all local similarities
- **Custom filtering**: When applying your own downstream filtering

## Implementation Details

The scaffold filtering is implemented in `filterByScaffolds()` in `computeMap.hpp`:

1. Creates a copy of input mappings for scaffold construction
2. Applies `mergeMappingsInRange()` with 2× the normal chain gap
3. Filters scaffolds by minimum length
4. Groups mappings by query/reference pair
5. Keeps mappings within scaffold boundaries + deviation

## Debugging Scaffold Behavior

To understand scaffold filtering behavior, uncomment the scaffold output in the code:
- Scaffolds are written to `scaf.paf`
- Shows the large-scale syntenic blocks identified
- Useful for understanding why certain mappings are kept/discarded