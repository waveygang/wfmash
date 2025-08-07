#!/bin/bash

# Build script that vendors all dependencies and builds wfmash statically
# Usage: ./build-vendor-everything.sh [jobs]

set -e

JOBS=${1:-16}
BUILD_DIR="build"

echo "Building wfmash with all vendored dependencies..."
echo "Using $JOBS parallel jobs"

# Clean build directory
if [ -d "$BUILD_DIR" ]; then
    echo "Removing existing build directory..."
    rm -rf "$BUILD_DIR"
fi

# Configure with vendored everything
echo "Configuring build with VENDOR_EVERYTHING=ON..."
cmake -H. -B"$BUILD_DIR" \
    -DVENDOR_EVERYTHING=ON \
    -DCMAKE_BUILD_TYPE=Release

# Build
echo "Building wfmash..."
cmake --build "$BUILD_DIR" -- -j "$JOBS"

# Get version and copy binary
if [ -f "$BUILD_DIR/bin/wfmash" ]; then
    VERSION=$("$BUILD_DIR/bin/wfmash" --version 2>&1 | head -n1 | awk '{print $NF}')
    OUTPUT_NAME="wfmash-$VERSION"
    
    echo "Build successful! Copying binary to ./$OUTPUT_NAME"
    cp "$BUILD_DIR/bin/wfmash" "./$OUTPUT_NAME"
    
    echo "Done! Binary available at: ./$OUTPUT_NAME"
    echo "File info:"
    ls -lh "./$OUTPUT_NAME"
    
    echo "Testing basic functionality:"
    "./$OUTPUT_NAME" --help | head -5
else
    echo "Error: wfmash binary not found in $BUILD_DIR/bin/"
    exit 1
fi