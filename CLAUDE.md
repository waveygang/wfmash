# CLAUDE.md - wfmash Development and Testing Guide

## Build Instructions

### Prerequisites
- CMake (version 3.12 or higher)
- C++ compiler with C++17 support
- zlib development libraries

### Building wfmash
```bash
# Create build directory
mkdir -p build
cd build

# Configure with CMake
cmake ..

# Build (use -j for parallel compilation)
make -j4

# The binary will be in build/bin/wfmash
```

### Clean Build
```bash
rm -rf build
mkdir build
cd build
cmake ..
make -j4
```

## Testing Protocol

### Test Data
The primary test file is `cerevisiae.chrV.fa.gz` which should be in the project root.

### Test Commands

#### 1. Raw Mappings (L2 Mapping)
Tests initial mapping generation without merging or filtering:
```bash
wfmash -m -t 8 -M -f none -S 0 cerevisiae.chrV.fa.gz > baseline.raw.paf
```

#### 2. Merged Mappings (Union-Find)
Tests the mapping merge algorithm without filtering:
```bash
wfmash -m -t 8 -f none -S 0 cerevisiae.chrV.fa.gz > baseline.merged.paf
```

#### 3. Full Pipeline
Tests complete mapping pipeline with all processing:
```bash
wfmash -m -t 8 cerevisiae.chrV.fa.gz > baseline.full.paf
```

### CTest Integration
If CTest is configured, you can run tests with:
```bash
cd build
ctest
```

## Comparing Outputs

### PAF Format Comparison
When comparing PAF outputs after changes:
```bash
# Sort and compare outputs (PAF columns can vary in order)
sort baseline.full.paf > baseline.sorted.paf
sort current.full.paf > current.sorted.paf
diff baseline.sorted.paf current.sorted.paf
```

### Key Fields to Check
- Query and target sequence names
- Start and end positions
- Strand orientation
- Mapping quality/identity scores

## Common Issues

### Missing Test Data
If `cerevisiae.chrV.fa.gz` is missing, download from:
```bash
# Example command to download test data
wget [URL_TO_TEST_DATA]/cerevisiae.chrV.fa.gz
```

### Memory Testing
For memory usage analysis:
```bash
# Run with valgrind
valgrind --tool=massif ./build/bin/wfmash -m -t 8 cerevisiae.chrV.fa.gz > /dev/null
ms_print massif.out.[PID]
```

## Development Workflow

### Creating Feature Branches
```bash
# Create new branch from current branch
git checkout -b feature-name

# Make changes and test
# ... edit files ...
make -j4
# ... run tests ...

# Commit changes
git add -p  # Interactive staging
git commit -m "Descriptive commit message"
```

### Testing Changes
After making changes:
1. Rebuild: `cd build && make -j4`
2. Run basic test: `./bin/wfmash -m -t 8 cerevisiae.chrV.fa.gz > test.paf`
3. Compare output: `diff baseline.full.paf test.paf`

## Current Work

### Recent Fixes
- **stdout mapping truncation**: Fixed issue where mapping output was truncated when writing to stdout due to repeated opening/closing of stdout for each reference subset (branch: `fix-stdout-mapping-truncation`)

### Git Workflow Pattern
Always checkpoint work in git commits with distinct changesets:
1. Work on single, focused changes
2. Commit to local non-main branches with descriptive names
3. Write clear commit messages explaining what and why