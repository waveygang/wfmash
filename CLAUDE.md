# CLAUDE.md - wfmash Development Guide

## What is wfmash?
wfmash is a tool for computing whole-genome alignments using mashmap for mapping and WFA (wavefront alignment) for base-level alignment.

## Repository Structure
```
wfmash/
├── src/
│   ├── map/           # Mapping implementation (mashmap)
│   │   └── include/   # Header files for mapping
│   ├── align/         # Alignment implementation (wflign)
│   ├── common/        # Shared utilities
│   └── interface/     # Main program interface
├── scripts/           # Utility scripts
├── build/            # Build directory (created by cmake)
└── tests/            # Test files
```

## Building wfmash

### Prerequisites
- CMake (version 3.12 or higher)
- C++ compiler with C++17 support
- zlib development libraries

### Build Commands
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

## Testing

### Test Data
The primary test file is `cerevisiae.chrV.fa.gz` which should be in the project root.

### Basic Test Commands

1. **Raw Mappings** (no merging, no filtering):
```bash
./build/bin/wfmash -m -t 8 -M -f none -S 0 cerevisiae.chrV.fa.gz > test.raw.paf
```

2. **Merged Mappings** (with merging, no filtering):
```bash
./build/bin/wfmash -m -t 8 -f none -S 0 cerevisiae.chrV.fa.gz > test.merged.paf
```

3. **Full Pipeline** (all processing):
```bash
./build/bin/wfmash -m -t 8 cerevisiae.chrV.fa.gz > test.full.paf
```

### Running CTest
```bash
cd build
ctest
```

## Key Command Options
- `-m`: Mapping only (no alignment)
- `-t N`: Use N threads
- `-M`: Disable merging
- `-f none`: Disable filtering
- `-S 0`: Disable scaffold filtering