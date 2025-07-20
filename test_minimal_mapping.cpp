/**
 * @file    test_minimal_mapping.cpp
 * @brief   Unit tests for MinimalMapping conversion functions
 * @author  Test Suite
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <iomanip>
#include "src/map/include/base_types.hpp"

using namespace skch;

// Helper function to compare floating point values
bool nearlyEqual(float a, float b, float epsilon = 0.01f) {
    return std::abs(a - b) < epsilon;
}

// Test basic conversion
void testBasicConversion() {
    std::cout << "Test 1: Basic conversion..." << std::endl;
    
    
    // Create a full MappingResult
    MappingResult original;
    original.querySeqId = 10;
    original.queryLen = 10000;
    original.queryStartPos = 100;
    original.queryEndPos = 600;
    original.refSeqId = 5;
    original.refStartPos = 200;
    original.refEndPos = 700;
    original.nucIdentity = 0.95f;
    original.nucIdentityUpperBound = 0.98f;
    original.blockLength = 500;
    original.strand = strnd::FWD;
    original.discard = 0;
    original.overlapped = false;
    original.kmerComplexity = 0.85;
    original.sketchSize = 100;
    original.conservedSketches = 95;
    
    // Compress
    MinimalMapping compressed = compressMapping(original);
    
    // Check compressed values
    assert(compressed.ref_seqId == 5);
    assert(compressed.ref_pos == 200);
    assert(compressed.query_pos == 100);
    assert(compressed.length == 500);
    assert(compressed.identity == 95); // 0.95 * 100
    assert(compressed.getStrand() == strnd::FWD);
    assert(!compressed.isDiscarded());
    assert(!compressed.isOverlapped());
    
    // Expand back
    MappingResult expanded = expandMinimalMapping(compressed, original.querySeqId, original.queryLen);
    
    // Check expanded values match original (for stored fields)
    assert(expanded.querySeqId == original.querySeqId);
    assert(expanded.queryLen == original.queryLen);
    assert(expanded.queryStartPos == original.queryStartPos);
    assert(expanded.queryEndPos == original.queryStartPos + compressed.length);
    assert(expanded.refSeqId == original.refSeqId);
    assert(expanded.refStartPos == original.refStartPos);
    assert(expanded.refEndPos == original.refStartPos + compressed.length);
    assert(nearlyEqual(expanded.nucIdentity, original.nucIdentity));
    assert(expanded.strand == original.strand);
    assert(expanded.discard == original.discard);
    assert(expanded.overlapped == original.overlapped);
    
    std::cout << "✓ Passed" << std::endl;
}

// Test strand encoding
void testStrandEncoding() {
    std::cout << "Test 2: Strand encoding..." << std::endl;
    
    MappingResult original;
    original.refSeqId = 1;
    original.refStartPos = 0;
    original.refEndPos = 100;
    original.queryStartPos = 0;
    original.queryEndPos = 100;
    original.nucIdentity = 1.0f;
    
    // Test FWD strand
    original.strand = strnd::FWD;
    MinimalMapping compressed = compressMapping(original);
    assert(compressed.getStrand() == strnd::FWD);
    
    // Test REV strand
    original.strand = strnd::REV;
    compressed = compressMapping(original);
    assert(compressed.getStrand() == strnd::REV);
    
    // Test AMBIG strand
    original.strand = strnd::AMBIG;
    compressed = compressMapping(original);
    assert(compressed.getStrand() == strnd::AMBIG);
    
    std::cout << "✓ Passed" << std::endl;
}

// Test flag encoding
void testFlagEncoding() {
    std::cout << "Test 3: Flag encoding..." << std::endl;
    
    MappingResult original;
    original.refSeqId = 1;
    original.refStartPos = 0;
    original.refEndPos = 100;
    original.queryStartPos = 0;
    original.queryEndPos = 100;
    original.nucIdentity = 1.0f;
    original.strand = strnd::FWD;
    
    // Test all flag combinations
    for (int discard = 0; discard <= 1; discard++) {
        for (int overlapped = 0; overlapped <= 1; overlapped++) {
            original.discard = discard;
            original.overlapped = (overlapped == 1);
            
            MinimalMapping compressed = compressMapping(original);
            assert(compressed.isDiscarded() == (discard == 1));
            assert(compressed.isOverlapped() == (overlapped == 1));
            assert(compressed.getStrand() == strnd::FWD); // Ensure strand not affected
        }
    }
    
    std::cout << "✓ Passed" << std::endl;
}

// Test identity precision
void testIdentityPrecision() {
    std::cout << "Test 4: Identity precision..." << std::endl;
    
    MappingResult original;
    original.refSeqId = 1;
    original.refStartPos = 0;
    original.refEndPos = 100;
    original.queryStartPos = 0;
    original.queryEndPos = 100;
    original.strand = strnd::FWD;
    
    // Test various identity values
    std::vector<float> identities = {0.0f, 0.005f, 0.01f, 0.5f, 0.99f, 0.995f, 1.0f};
    
    for (float id : identities) {
        original.nucIdentity = id;
        MinimalMapping compressed = compressMapping(original);
        uint8_t expected = static_cast<uint8_t>(std::round(id * 100));
        assert(compressed.identity == expected);
        
        // Verify precision loss is within 0.5%
        MappingResult expanded = expandMinimalMapping(compressed, 1, 1000);
        assert(nearlyEqual(expanded.nucIdentity, expected / 100.0f, 0.005f));
    }
    
    std::cout << "✓ Passed" << std::endl;
}

// Test memory size
void testMemorySize() {
    std::cout << "Test 5: Memory size..." << std::endl;
    
    std::cout << "  MappingResult size: " << sizeof(MappingResult) << " bytes" << std::endl;
    std::cout << "  MinimalMapping size: " << sizeof(MinimalMapping) << " bytes" << std::endl;
    
    // Verify size reduction
    assert(sizeof(MinimalMapping) <= 16); // Should be 16 bytes with padding
    assert(sizeof(MinimalMapping) < sizeof(MappingResult) / 8); // At least 8x reduction
    
    std::cout << "✓ Passed (>" << (sizeof(MappingResult) / sizeof(MinimalMapping)) << "x reduction)" << std::endl;
}

// Test edge cases
void testEdgeCases() {
    std::cout << "Test 6: Edge cases..." << std::endl;
    
    MappingResult original;
    
    // Test zero length
    original.refSeqId = 0;
    original.refStartPos = 100;
    original.refEndPos = 100;
    original.queryStartPos = 50;
    original.queryEndPos = 50;
    original.nucIdentity = 1.0f;
    original.strand = strnd::FWD;
    
    MinimalMapping compressed = compressMapping(original);
    assert(compressed.length == 0);
    
    // Test maximum values that fit in uint32_t
    original.refSeqId = 0xFFFFFFFF;
    original.refStartPos = 0xFFFFFFFF;
    original.queryStartPos = 0xFFFFFFFF;
    original.refEndPos = 0xFFFFFFFF;
    original.queryEndPos = 0xFFFFFFFF;
    
    compressed = compressMapping(original);
    assert(compressed.ref_seqId == 0xFFFFFFFF);
    assert(compressed.ref_pos == 0xFFFFFFFF);
    assert(compressed.query_pos == 0xFFFFFFFF);
    
    // Test maximum length that fits in uint16_t
    original.refStartPos = 0;
    original.refEndPos = 65535;
    original.queryStartPos = 0;
    original.queryEndPos = 65535;
    
    compressed = compressMapping(original);
    assert(compressed.length == 65535);
    
    std::cout << "✓ Passed" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "=== MinimalMapping Unit Tests ===" << std::endl;
    
    try {
        testBasicConversion();
        testStrandEncoding();
        testFlagEncoding();
        testIdentityPrecision();
        testMemorySize();
        testEdgeCases();
        
        std::cout << "\nAll tests passed! ✓" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << std::endl;
        return 1;
    }
}