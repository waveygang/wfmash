# Testing Documentation for wfmash

## Running Tests

### Quick Start
To run all tests after building the project:
```bash
./scripts/run_all_tests.sh
```

### Using CTest Directly
From the build directory:
```bash
cd build
ctest                    # Run all tests
ctest -V                 # Run with verbose output
ctest -R <test-name>     # Run specific test by name
ctest --output-on-failure # Show output only for failed tests
```

## Available Tests

1. **wfmash-time-LPA**: Tests timing and performance with LPA data
2. **wfmash-subset-LPA-to-SAM**: Tests subset functionality with SAM output
3. **wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF**: Tests mapping coverage using yeast genomes
4. **wfmash-pafcheck-yeast**: Validates PAF output format with yeast data
5. **wfmash-maf-validity**: Tests MAF format output validity
6. **wfmash-multi-subset-index**: Tests multi-subset indexing functionality

## Code Refactoring

The `computeMap.hpp` file has been refactored into 5 modular components:

### 1. fragmentManager.hpp
- Fragment data structures (`FragmentData`, `QueryMappingOutput`)
- Fragment management utilities

### 2. mappingCore.hpp
- Core L1/L2 mapping algorithms
- `L1_candidateLocus_t` and `L2_mapLocus_t` structures
- Seed hit computation
- Interval point calculation
- L1 candidate region computation
- L2 mapped region computation

### 3. mappingFilter.hpp
- Mapping filtering algorithms
- Merging logic
- Weak mapping filtering
- Identity filtering
- Sparsification
- Scaffold filtering
- Chain processing

### 4. mappingOutput.hpp
- Result processing
- Boundary sanity checking
- Mapping report generation
- Output formatting

### 5. computeMap.hpp (main)
- Main `Map` class
- Integration of all components
- Query processing workflow
- Subset management

## Building and Testing After Refactoring

1. Clean build:
```bash
rm -rf build
cmake -H. -Bbuild
cmake --build build -- -j 16
```

2. Run tests:
```bash
./scripts/run_all_tests.sh
```

## Verification

The refactoring has been verified to:
- Reduce the main file from 3191 lines to 1001 lines (69% reduction)
- Maintain all functionality
- Pass compilation without errors
- Preserve the existing test suite