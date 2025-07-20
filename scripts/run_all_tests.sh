#!/bin/bash
# Test harness for wfmash
# Generated from CMake test configuration

set -e  # Exit on any error

echo "=== WFMASH TEST HARNESS ==="
echo "Running all available tests..."
echo

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Change to build directory
cd build || { echo "Build directory not found. Please run cmake and make first."; exit 1; }

# Function to run a test and report result
run_test() {
    local test_name=$1
    echo -n "Running test: $test_name ... "
    
    if ctest -R "^$test_name$" --output-on-failure > test_output.log 2>&1; then
        echo -e "${GREEN}PASSED${NC}"
        return 0
    else
        echo -e "${RED}FAILED${NC}"
        echo "Error output:"
        cat test_output.log
        return 1
    fi
}

# Counter for test results
total_tests=0
passed_tests=0
failed_tests=0

# List of tests (from ctest -N output)
tests=(
    "wfmash-time-LPA"
    "wfmash-subset-LPA-to-SAM"
    "wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF"
    "wfmash-pafcheck-yeast"
    "wfmash-maf-validity"
    "wfmash-multi-subset-index"
)

# Run each test
for test in "${tests[@]}"; do
    ((total_tests++))
    if run_test "$test"; then
        ((passed_tests++))
    else
        ((failed_tests++))
    fi
    echo
done

# Clean up
rm -f test_output.log

# Summary
echo "=== TEST SUMMARY ==="
echo "Total tests: $total_tests"
echo -e "Passed: ${GREEN}$passed_tests${NC}"
echo -e "Failed: ${RED}$failed_tests${NC}"

if [ $failed_tests -eq 0 ]; then
    echo -e "\n${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "\n${RED}Some tests failed!${NC}"
    exit 1
fi