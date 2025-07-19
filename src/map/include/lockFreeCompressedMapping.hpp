/**
 * @file    lockFreeCompressedMapping.hpp
 * @brief   Lock-free compressed mapping storage for maximum performance
 * @author  Memory Optimization
 */

#ifndef LOCK_FREE_COMPRESSED_MAPPING_HPP
#define LOCK_FREE_COMPRESSED_MAPPING_HPP

#include <vector>
#include <algorithm>
#include <atomic>
#include <memory>
#include "map/include/base_types.hpp"

namespace skch
{
  /**
   * @class LockFreeCompressedMappingStore
   * @brief Lock-free storage for mappings using atomic operations
   *
   * This class provides a memory-efficient, lock-free storage for mapping results.
   * It uses a pre-allocated buffer and atomic index for thread-safe appending.
   */
  class LockFreeCompressedMappingStore {
  private:
    // Pre-allocated storage
    std::unique_ptr<MinimalMapping[]> mappings;
    size_t capacity;
    std::atomic<size_t> size{0};
    
    // Query information
    seqno_t querySeqId;
    offset_t queryLen;
    
    // Metadata (set once, read many)
    std::atomic<float> defaultKmerComplexity{0.0};
    std::atomic<int> defaultSketchSize{0};
    
  public:
    // Constructor with pre-allocated capacity
    explicit LockFreeCompressedMappingStore(size_t initialCapacity = 10000, 
                                           seqno_t qSeqId = 0, 
                                           offset_t qLen = 0) 
        : capacity(initialCapacity), querySeqId(qSeqId), queryLen(qLen) {
      mappings = std::make_unique<MinimalMapping[]>(capacity);
    }
    
    // Lock-free add single mapping
    bool addMapping(const MappingResult& m) {
      // Update query info if not set (first mapping sets it)
      if (querySeqId == 0 && m.querySeqId != 0) {
        querySeqId = m.querySeqId;
        queryLen = m.queryLen;
      }
      
      // Atomically claim a slot
      size_t index = size.fetch_add(1, std::memory_order_relaxed);
      
      if (index >= capacity) {
        // Buffer full - in production, could implement wait-free growth
        size.fetch_sub(1, std::memory_order_relaxed);
        return false;
      }
      
      // Write to our exclusive slot (no contention)
      mappings[index] = compressMapping(m);
      
      // Update metadata if needed (race condition acceptable for these)
      if (defaultSketchSize.load(std::memory_order_relaxed) == 0 && m.sketchSize > 0) {
        defaultSketchSize.store(m.sketchSize, std::memory_order_relaxed);
      }
      if (defaultKmerComplexity.load(std::memory_order_relaxed) == 0.0 && m.kmerComplexity > 0.0) {
        defaultKmerComplexity.store(m.kmerComplexity, std::memory_order_relaxed);
      }
      
      return true;
    }
    
    // Batch add - each thread can add its own batch without contention
    size_t addMappingsBatch(const std::vector<MappingResult>& batch) {
      if (batch.empty()) return 0;
      
      // Atomically reserve space for entire batch
      size_t batchSize = batch.size();
      size_t startIndex = size.fetch_add(batchSize, std::memory_order_relaxed);
      
      if (startIndex + batchSize > capacity) {
        // Not enough space - roll back
        size.fetch_sub(batchSize, std::memory_order_relaxed);
        return 0;
      }
      
      // Write entire batch to our exclusive range
      for (size_t i = 0; i < batchSize; ++i) {
        mappings[startIndex + i] = compressMapping(batch[i]);
        
        // Update metadata for first few mappings
        if (i < 10) {
          if (defaultSketchSize.load(std::memory_order_relaxed) == 0 && batch[i].sketchSize > 0) {
            defaultSketchSize.store(batch[i].sketchSize, std::memory_order_relaxed);
          }
          if (defaultKmerComplexity.load(std::memory_order_relaxed) == 0.0 && batch[i].kmerComplexity > 0.0) {
            defaultKmerComplexity.store(batch[i].kmerComplexity, std::memory_order_relaxed);
          }
        }
      }
      
      return batchSize;
    }
    
    // Get all mappings (read-only after parallel phase)
    std::vector<MappingResult> getAllMappings() const {
      size_t currentSize = size.load(std::memory_order_acquire);
      std::vector<MappingResult> results;
      results.reserve(currentSize);
      
      int sketchSize = defaultSketchSize.load(std::memory_order_relaxed);
      float kmerComplexity = defaultKmerComplexity.load(std::memory_order_relaxed);
      
      for (size_t i = 0; i < currentSize; ++i) {
        MappingResult result = expandMinimalMapping(mappings[i], querySeqId, queryLen);
        result.sketchSize = sketchSize;
        result.kmerComplexity = kmerComplexity;
        results.push_back(result);
      }
      
      return results;
    }
    
    // Get compressed mappings directly
    std::vector<MinimalMapping> getCompressedMappings() const {
      size_t currentSize = size.load(std::memory_order_acquire);
      std::vector<MinimalMapping> results;
      results.reserve(currentSize);
      
      for (size_t i = 0; i < currentSize; ++i) {
        results.push_back(mappings[i]);
      }
      
      return results;
    }
    
    // Size and capacity
    size_t getSize() const { 
      return size.load(std::memory_order_acquire);
    }
    
    size_t getCapacity() const { 
      return capacity;
    }
    
    bool empty() const { 
      return getSize() == 0;
    }
    
    void clear() {
      size.store(0, std::memory_order_release);
    }
    
    // Reserve capacity (must be called before parallel phase)
    void reserve(size_t newCapacity) {
      if (newCapacity > capacity && size.load(std::memory_order_acquire) == 0) {
        mappings = std::make_unique<MinimalMapping[]>(newCapacity);
        capacity = newCapacity;
      }
    }
    
    // Memory statistics
    size_t memoryUsage() const {
      return capacity * sizeof(MinimalMapping) + sizeof(*this);
    }
    
    size_t actualMemoryUsage() const {
      return getSize() * sizeof(MinimalMapping) + sizeof(*this);
    }
    
    // Sort operations (only safe after parallel phase completes)
    void sortByQueryPos() {
      size_t currentSize = size.load(std::memory_order_acquire);
      std::sort(mappings.get(), mappings.get() + currentSize,
                [](const MinimalMapping& a, const MinimalMapping& b) {
                  return a.query_pos < b.query_pos;
                });
    }
    
    void sortByRefPos() {
      size_t currentSize = size.load(std::memory_order_acquire);
      std::sort(mappings.get(), mappings.get() + currentSize,
                [](const MinimalMapping& a, const MinimalMapping& b) {
                  return std::tie(a.ref_seqId, a.ref_pos) < 
                         std::tie(b.ref_seqId, b.ref_pos);
                });
    }
  };
}

#endif // LOCK_FREE_COMPRESSED_MAPPING_HPP