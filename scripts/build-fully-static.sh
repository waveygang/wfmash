#!/bin/bash

# Build script that attempts fully static linking (experimental)
# Usage: ./build-fully-static.sh [jobs]
# Note: This may fail on HPC systems with mandatory dynamic libraries

set -e

JOBS=${1:-16}
BUILD_DIR="build-fully-static"

echo "Building wfmash with FULLY_STATIC=ON (experimental)..."
echo "Using $JOBS parallel jobs"
echo "WARNING: This may fail due to system library dependencies"

# Clean build directory
if [ -d "$BUILD_DIR" ]; then
    echo "Removing existing build directory..."
    rm -rf "$BUILD_DIR"
fi

# Configure with fully static linking
echo "Configuring build with FULLY_STATIC=ON..."
cmake -H. -B"$BUILD_DIR" \
    -DVENDOR_EVERYTHING=ON \
    -DFULLY_STATIC=ON \
    -DCMAKE_BUILD_TYPE=Release

# Build
echo "Building wfmash..."
if cmake --build "$BUILD_DIR" -- -j "$JOBS"; then
    # Get version and copy binary
    if [ -f "$BUILD_DIR/bin/wfmash" ]; then
        VERSION=$("$BUILD_DIR/bin/wfmash" --version 2>&1 | head -n1 | awk '{print $NF}')
        OUTPUT_NAME="wfmash-$VERSION-fully-static"
        
        echo "Build successful! Copying binary to ./$OUTPUT_NAME"
        cp "$BUILD_DIR/bin/wfmash" "./$OUTPUT_NAME"
        
        echo "Done! Binary available at: ./$OUTPUT_NAME"
        echo "File info:"
        ls -lh "./$OUTPUT_NAME"
        
        echo "Dynamic library dependencies (should be minimal):"
        ldd "./$OUTPUT_NAME" || echo "Fully static - no dependencies!"
        
        echo "Testing basic functionality:"
        "./$OUTPUT_NAME" --help | head -5
    else
        echo "Error: wfmash binary not found in $BUILD_DIR/bin/"
        exit 1
    fi
else
    echo ""
    echo "Build failed as expected on this system."
    echo "Common reasons:"
    echo "  - Static versions of system libraries (libc, libm) not available"
    echo "  - HPC monitoring libraries (libxalt) cannot be statically linked"
    echo "  - OpenMP runtime requires dynamic linking"
    echo ""
    echo "The regular BUILD_STATIC option provides a good compromise."
    exit 1
fi