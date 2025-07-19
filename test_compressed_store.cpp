/**
 * @file    test_compressed_store.cpp
 * @brief   Unit tests for CompressedMappingStore
 * @author  Test Suite
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include "src/map/include/base_types.hpp"
#include "src/map/include/compressedMapping.hpp"

using namespace skch;

// Helper function to create a test MappingResult
MappingResult createTestMapping(seqno_t querySeqId, seqno_t refSeqId, 
                               offset_t queryStart, offset_t refStart, 
                               offset_t length, float identity) {
    MappingResult m;
    m.querySeqId = querySeqId;
    m.queryLen = 10000;
    m.queryStartPos = queryStart;
    m.queryEndPos = queryStart + length;
    m.refSeqId = refSeqId;
    m.refStartPos = refStart;
    m.refEndPos = refStart + length;
    m.nucIdentity = identity;
    m.nucIdentityUpperBound = identity + 0.02f;
    m.blockLength = length;
    m.blockNucIdentity = identity;
    m.strand = (refStart % 2 == 0) ? strnd::FWD : strnd::REV;
    m.discard = 0;
    m.overlapped = false;
    m.sketchSize = 100;
    m.conservedSketches = static_cast<int>(100 * identity);
    m.kmerComplexity = 0.85;
    m.approxMatches = static_cast<int>(length * identity);
    m.n_merged = 1;
    m.splitMappingId = 0;
    m.selfMapFilter = false;
    m.chainPairScore = std::numeric_limits<double>::max();
    m.chainPairId = std::numeric_limits<int64_t>::min();
    m.chain_id = -1;
    m.chain_length = 1;
    m.chain_pos = 1;
    
    return m;
}

// Test basic store and retrieve
void testBasicStoreRetrieve() {
    std::cout << "Test 1: Basic store and retrieve..." << std::endl;
    
    CompressedMappingStore store(1, 10000);
    
    // Add a single mapping
    MappingResult original = createTestMapping(1, 5, 100, 200, 500, 0.95f);
    store.addMapping(original);
    
    assert(store.size() == 1);
    assert(!store.empty());
    
    // Retrieve and verify
    MappingResult retrieved = store.getMapping(0);
    
    assert(retrieved.querySeqId == original.querySeqId);
    assert(retrieved.queryLen == original.queryLen);
    assert(retrieved.queryStartPos == original.queryStartPos);
    assert(retrieved.queryEndPos == original.queryEndPos);
    assert(retrieved.refSeqId == original.refSeqId);
    assert(retrieved.refStartPos == original.refStartPos);
    assert(retrieved.refEndPos == original.refEndPos);
    assert(std::abs(retrieved.nucIdentity - original.nucIdentity) < 0.01f);
    assert(retrieved.strand == original.strand);
    assert(retrieved.discard == original.discard);
    assert(retrieved.overlapped == original.overlapped);
    
    std::cout << "✓ Passed" << std::endl;
}

// Test multiple mappings
void testMultipleMappings() {
    std::cout << "Test 2: Multiple mappings..." << std::endl;
    
    CompressedMappingStore store;
    std::vector<MappingResult> originals;
    
    // Add 100 mappings
    for (int i = 0; i < 100; i++) {
        MappingResult m = createTestMapping(1, i % 10, i * 100, i * 200, 
                                           100 + i * 10, 0.80f + (i % 20) * 0.01f);
        originals.push_back(m);
        store.addMapping(m);
    }
    
    assert(store.size() == 100);
    
    // Retrieve all and verify
    std::vector<MappingResult> retrieved = store.getAllMappings();
    assert(retrieved.size() == 100);
    
    for (size_t i = 0; i < 100; i++) {
        assert(retrieved[i].queryStartPos == originals[i].queryStartPos);
        assert(retrieved[i].refSeqId == originals[i].refSeqId);
        assert(retrieved[i].refStartPos == originals[i].refStartPos);
        assert(std::abs(retrieved[i].nucIdentity - originals[i].nucIdentity) < 0.01f);
    }
    
    std::cout << "✓ Passed" << std::endl;
}

// Test batch add
void testBatchAdd() {
    std::cout << "Test 3: Batch add mappings..." << std::endl;
    
    CompressedMappingStore store;
    std::vector<MappingResult> batch;
    
    // Create batch
    for (int i = 0; i < 50; i++) {
        batch.push_back(createTestMapping(2, i, i * 50, i * 60, 200, 0.90f));
    }
    
    store.addMappings(batch);
    assert(store.size() == 50);
    
    // Verify first and last
    MappingResult first = store.getMapping(0);
    MappingResult last = store.getMapping(49);
    
    assert(first.queryStartPos == 0);
    assert(last.queryStartPos == 49 * 50);
    
    std::cout << "✓ Passed" << std::endl;
}

// Test discard functionality
void testDiscardFunctionality() {
    std::cout << "Test 4: Discard functionality..." << std::endl;
    
    CompressedMappingStore store;
    
    // Add 10 mappings
    for (int i = 0; i < 10; i++) {
        store.addMapping(createTestMapping(1, i, i * 100, i * 100, 100, 0.95f));
    }
    
    // Mark some as discarded
    store.markDiscard(2);
    store.markDiscard(5);
    store.markDiscard(8);
    
    // Check discarded flags
    assert(!store.getMapping(0).discard);
    assert(store.getMapping(2).discard);
    assert(store.getMapping(5).discard);
    assert(!store.getMapping(7).discard);
    assert(store.getMapping(8).discard);
    
    // Remove discarded
    size_t sizeBefore = store.size();
    store.removeDiscarded();
    assert(store.size() == sizeBefore - 3);
    
    std::cout << "✓ Passed" << std::endl;
}

// Test sorting
void testSorting() {
    std::cout << "Test 5: Sorting functionality..." << std::endl;
    
    CompressedMappingStore store;
    
    // Add mappings in random order
    std::vector<int> positions = {500, 100, 300, 200, 400};
    for (int pos : positions) {
        store.addMapping(createTestMapping(1, 1, pos, pos * 2, 100, 0.95f));
    }
    
    // Sort by query position
    store.sortByQueryPos();
    
    // Verify sorted order
    for (size_t i = 1; i < store.size(); i++) {
        assert(store.getMapping(i-1).queryStartPos <= store.getMapping(i).queryStartPos);
    }
    
    // Sort by reference position
    store.sortByRefPos();
    
    // Verify sorted order
    for (size_t i = 1; i < store.size(); i++) {
        MappingResult prev = store.getMapping(i-1);
        MappingResult curr = store.getMapping(i);
        assert(std::tie(prev.refSeqId, prev.refStartPos) <= 
               std::tie(curr.refSeqId, curr.refStartPos));
    }
    
    std::cout << "✓ Passed" << std::endl;
}

// Test memory efficiency
void testMemoryEfficiency() {
    std::cout << "Test 6: Memory efficiency..." << std::endl;
    
    CompressedMappingStore store;
    
    // Add 1000 mappings
    for (int i = 0; i < 1000; i++) {
        store.addMapping(createTestMapping(1, i % 100, i * 10, i * 20, 500, 0.85f + (i % 15) * 0.01f));
    }
    
    size_t compressed = store.memoryUsage();
    size_t uncompressed = store.uncompressedMemoryUsage();
    float ratio = store.compressionRatio();
    
    std::cout << "  Compressed size: " << compressed << " bytes" << std::endl;
    std::cout << "  Uncompressed size: " << uncompressed << " bytes" << std::endl;
    std::cout << "  Compression ratio: " << ratio << "x" << std::endl;
    
    assert(ratio > 10.0f);  // Should achieve at least 10x compression
    
    std::cout << "✓ Passed" << std::endl;
}

// Test forEach functionality
void testForEach() {
    std::cout << "Test 7: forEach functionality..." << std::endl;
    
    CompressedMappingStore store;
    
    // Add mappings
    for (int i = 0; i < 10; i++) {
        store.addMapping(createTestMapping(1, i, i * 100, i * 100, 100, 0.90f));
    }
    
    // Count mappings with identity > 0.89
    int count = 0;
    store.forEach([&count](const MappingResult& m, size_t idx) {
        if (m.nucIdentity > 0.89f) {
            count++;
        }
    });
    
    assert(count == 10);  // All should have identity 0.90
    
    std::cout << "✓ Passed" << std::endl;
}

// Test thread safety (basic)
void testThreadSafety() {
    std::cout << "Test 8: Thread safety (basic)..." << std::endl;
    
    CompressedMappingStore store;
    const int mappingsPerThread = 100;
    const int numThreads = 4;
    
    std::vector<std::thread> threads;
    
    // Launch threads to add mappings
    for (int t = 0; t < numThreads; t++) {
        threads.emplace_back([&store, t, mappingsPerThread]() {
            for (int i = 0; i < mappingsPerThread; i++) {
                store.addMapping(createTestMapping(t, i, i * 100, i * 100, 100, 0.95f));
            }
        });
    }
    
    // Wait for all threads
    for (auto& thread : threads) {
        thread.join();
    }
    
    assert(store.size() == mappingsPerThread * numThreads);
    
    std::cout << "✓ Passed" << std::endl;
}

// Performance benchmark
void performanceBenchmark() {
    std::cout << "\nPerformance Benchmark:" << std::endl;
    
    const int numMappings = 100000;
    
    // Time uncompressed storage
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<MappingResult> uncompressed;
    uncompressed.reserve(numMappings);
    
    for (int i = 0; i < numMappings; i++) {
        uncompressed.push_back(createTestMapping(1, i % 1000, i * 10, i * 20, 500, 0.85f));
    }
    
    auto uncompressedTime = std::chrono::high_resolution_clock::now() - start;
    
    // Time compressed storage
    start = std::chrono::high_resolution_clock::now();
    CompressedMappingStore compressed;
    compressed.reserve(numMappings);
    
    for (int i = 0; i < numMappings; i++) {
        compressed.addMapping(createTestMapping(1, i % 1000, i * 10, i * 20, 500, 0.85f));
    }
    
    auto compressedTime = std::chrono::high_resolution_clock::now() - start;
    
    // Time retrieval
    start = std::chrono::high_resolution_clock::now();
    std::vector<MappingResult> retrieved = compressed.getAllMappings();
    auto retrievalTime = std::chrono::high_resolution_clock::now() - start;
    
    std::cout << "  Uncompressed storage time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(uncompressedTime).count() 
              << " ms" << std::endl;
    std::cout << "  Compressed storage time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(compressedTime).count() 
              << " ms" << std::endl;
    std::cout << "  Retrieval time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(retrievalTime).count() 
              << " ms" << std::endl;
    std::cout << "  Memory saved: " 
              << (uncompressed.size() * sizeof(MappingResult) - compressed.memoryUsage()) / (1024.0 * 1024.0) 
              << " MB" << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "=== CompressedMappingStore Unit Tests ===" << std::endl;
    
    try {
        testBasicStoreRetrieve();
        testMultipleMappings();
        testBatchAdd();
        testDiscardFunctionality();
        testSorting();
        testMemoryEfficiency();
        testForEach();
        testThreadSafety();
        
        performanceBenchmark();
        
        std::cout << "\nAll tests passed! ✓" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nTest failed with exception: " << e.what() << std::endl;
        return 1;
    }
}