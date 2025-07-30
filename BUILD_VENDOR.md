# Vendor Everything Build

This branch adds a `VENDOR_EVERYTHING` option to build wfmash with all dependencies downloaded and compiled from source.

## Quick Start

```bash
./build-vendor-everything.sh [jobs]
```

This script will:
1. Download and build all dependencies (zlib, bzip2, xz, htslib, gsl, libdeflate)
2. Build wfmash statically linked against these dependencies
3. Copy the final binary to `./wfmash-<version>`

## Manual Build

```bash
rm -rf build
cmake -H. -Bbuild -DVENDOR_EVERYTHING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -- -j 16
```

## What VENDOR_EVERYTHING Does

When `VENDOR_EVERYTHING=ON`, CMake automatically enables:
- `BUILD_STATIC=ON` - Build static binary
- `VENDOR_HTSLIB=ON` - Download and build htslib from source
- `BUILD_DEPS=ON` - Download and build all dependencies from source

Additional dependencies downloaded and built:
- zlib 1.3.1
- bzip2 1.0.8
- xz 5.4.6
- GSL 2.8
- libdeflate 1.20
- htslib 1.20

## Benefits

- No system dependency requirements beyond basic build tools
- Reproducible builds across different systems
- Static linking for easier deployment
- Known working versions of all dependencies

## Build Requirements

- CMake 3.9+
- GCC/Clang with C++17 support
- Basic build tools (make, autotools)
- Internet connection for downloading dependencies