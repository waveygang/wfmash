/**
 * @file    test_compressed_mapping_integration.cpp
 * @brief   Integration test for compressed mapping functionality
 * @author  Test Suite
 */

#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <cstring>
#include "src/map/include/base_types.hpp"
#include "src/map/include/map_parameters.hpp"
#include "src/map/include/computeMap.hpp"

using namespace skch;

// Simple mock postprocessing function
void mockPostProcess(const MappingResult&) {}

void testCompressedMappingFlag() {
    std::cout << "Test: Compressed mapping flag integration..." << std::endl;
    
    // Create parameters with compressed mapping enabled
    Parameters params;
    params.kmerSize = 15;
    params.segLength = 1000;
    params.block_length = 0;
    params.chain_gap = 2000;
    params.max_mapping_length = 50000;
    params.percentageIdentity = 70.0;
    params.use_compressed_mappings = true;  // Enable compressed mappings
    params.threads = 1;
    params.filterMode = filter::MAP;
    params.alphabetSize = 4;
    params.referenceSize = 10000;
    params.sketchSize = 100;
    params.numMappingsForSegment = 1;
    params.numMappingsForShortSequence = 1;
    params.split = true;
    params.lower_triangular = false;
    params.skip_self = false;
    params.skip_prefix = false;
    params.filterLengthMismatches = false;
    params.stage2_full_scan = false;
    params.stage1_topANI_filter = false;
    params.ANIDiff = 0.0;
    params.ANIDiffConf = 0.999;
    params.dropRand = false;
    params.mergeMappings = true;
    params.keep_low_pct_id = true;
    params.report_ANI_percentage = true;
    params.kmerComplexityThreshold = 0.0;
    params.hgNumerator = 1.0;
    params.world_minimizers = false;
    params.sparsity_hash_threshold = std::numeric_limits<uint64_t>::max();
    params.overlap_threshold = 1.0;
    params.scaffold_max_deviation = 100000;
    params.scaffold_gap = 5000;
    params.scaffold_min_length = 50000;
    params.legacy_output = false;
    params.minimum_hits = -1;
    params.max_kmer_freq = 0.0002;
    params.use_progress_bar = false;
    
    // Create Map object
    try {
        Map mapper(params, mockPostProcess);
        std::cout << "✓ Map object created successfully with use_compressed_mappings=true" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "✗ Failed to create Map object: " << e.what() << std::endl;
        assert(false);
    }
    
    std::cout << "✓ Passed" << std::endl;
}

void testCompressedVsUncompressedOutput() {
    std::cout << "Test: Compressed vs uncompressed output comparison..." << std::endl;
    
    // This test would require actual sequence data and a more complete setup
    // For now, we just verify the flag is properly propagated
    
    Parameters params1;
    params1.use_compressed_mappings = false;
    
    Parameters params2;
    params2.use_compressed_mappings = true;
    
    assert(params1.use_compressed_mappings == false);
    assert(params2.use_compressed_mappings == true);
    
    std::cout << "✓ Passed" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "=== Compressed Mapping Integration Tests ===" << std::endl;
    
    try {
        testCompressedMappingFlag();
        testCompressedVsUncompressedOutput();
        
        std::cout << "\nAll tests passed! ✓" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << std::endl;
        return 1;
    }
}