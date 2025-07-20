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
#include "base_types.hpp"

namespace skch
{
  /**
   * @class LockFreeCompressedMappingStore
   * @brief Lock-free storage for mappings using atomic operations
   *
   * This class provides a memory-efficient, lock-free storage for mapping results.
   * It uses dynamically growing storage with lock-free operations.
   */
  class LockFreeCompressedMappingStore {
  private:
    // Dynamic storage with growth capability
    struct StorageBlock {
      std::unique_ptr<MinimalMapping[]> data;
      size_t capacity;
      std::atomic<size_t> used{0};
      
      StorageBlock(size_t cap) : capacity(cap) {
        data = std::make_unique<MinimalMapping[]>(capacity);
      }
    };
    
    std::vector<std::unique_ptr<StorageBlock>> blocks;
    std::mutex growth_mutex;  // Only for growing, not for adding
    size_t block_size;
    std::atomic<size_t> total_size{0};
    
    // Query information
    seqno_t querySeqId;
    offset_t queryLen;
    
    // Metadata (set once, read many)
    std::atomic<float> defaultKmerComplexity{0.0};
    std::atomic<int> defaultSketchSize{0};
    
  public:
    // Constructor with initial block size
    explicit LockFreeCompressedMappingStore(size_t initialCapacity = 1000,  // Start small!
                                           seqno_t qSeqId = 0, 
                                           offset_t qLen = 0) 
        : block_size(initialCapacity), querySeqId(qSeqId), queryLen(qLen) {
      // Start with one block
      blocks.emplace_back(std::make_unique<StorageBlock>(block_size));
    }
    
    // Lock-free add single mapping with automatic growth
    bool addMapping(const MappingResult& m) {
      // Update query info if not set (first mapping sets it)
      if (querySeqId == 0 && m.querySeqId != 0) {
        querySeqId = m.querySeqId;
        queryLen = m.queryLen;
      }
      
      // Try to add to existing blocks
      while (true) {
        // Try each block in order
        for (size_t i = 0; i < blocks.size(); ++i) {
          auto& block = blocks[i];
          size_t idx = block->used.fetch_add(1, std::memory_order_relaxed);
          
          if (idx < block->capacity) {
            // Success - we have a slot
            block->data[idx] = compressMapping(m);
            total_size.fetch_add(1, std::memory_order_relaxed);
            
            // Update metadata if needed
            if (defaultSketchSize.load(std::memory_order_relaxed) == 0 && m.sketchSize > 0) {
              defaultSketchSize.store(m.sketchSize, std::memory_order_relaxed);
            }
            if (defaultKmerComplexity.load(std::memory_order_relaxed) == 0.0 && m.kmerComplexity > 0.0) {
              defaultKmerComplexity.store(m.kmerComplexity, std::memory_order_relaxed);
            }
            
            return true;
          } else {
            // This block is full, undo the increment
            block->used.fetch_sub(1, std::memory_order_relaxed);
          }
        }
        
        // All blocks are full, need to grow
        std::lock_guard<std::mutex> lock(growth_mutex);
        
        // Check again in case another thread already grew
        if (blocks.back()->used.load() < blocks.back()->capacity) {
          continue;  // Retry with the new block
        }
        
        // Add a new block with exponential growth
        size_t current_total = 0;
        for (const auto& block : blocks) {
          current_total += block->capacity;
        }
        
        // Growth strategy: double total capacity each time, but cap individual blocks
        size_t new_block_size = std::min(current_total, size_t(10000000)); // Max 10M per block
        blocks.emplace_back(std::make_unique<StorageBlock>(new_block_size));
      }
    }
    
    // Batch add - more efficient than individual adds
    size_t addMappingsBatch(const std::vector<MappingResult>& batch) {
      if (batch.empty()) return 0;
      
      // Just use individual adds for now - the overhead is minimal
      // and it ensures we never lose mappings
      size_t added = 0;
      for (const auto& m : batch) {
        if (addMapping(m)) {
          added++;
        }
      }
      
      return added;
    }
    
    // Get all mappings (read-only after parallel phase)
    std::vector<MappingResult> getAllMappings() const {
      size_t currentSize = total_size.load(std::memory_order_acquire);
      std::vector<MappingResult> results;
      results.reserve(currentSize);
      
      int sketchSize = defaultSketchSize.load(std::memory_order_relaxed);
      float kmerComplexity = defaultKmerComplexity.load(std::memory_order_relaxed);
      
      // Collect from all blocks
      for (const auto& block : blocks) {
        size_t block_used = block->used.load(std::memory_order_acquire);
        for (size_t i = 0; i < block_used && i < block->capacity; ++i) {
          MappingResult result = expandMinimalMapping(block->data[i], querySeqId, queryLen);
          result.sketchSize = sketchSize;
          result.kmerComplexity = kmerComplexity;
          results.push_back(result);
        }
      }
      
      return results;
    }
    
    // Get compressed mappings directly
    std::vector<MinimalMapping> getCompressedMappings() const {
      size_t currentSize = total_size.load(std::memory_order_acquire);
      std::vector<MinimalMapping> results;
      results.reserve(currentSize);
      
      // Collect from all blocks
      for (const auto& block : blocks) {
        size_t block_used = block->used.load(std::memory_order_acquire);
        for (size_t i = 0; i < block_used && i < block->capacity; ++i) {
          results.push_back(block->data[i]);
        }
      }
      
      return results;
    }
    
    // Size and capacity
    size_t getSize() const { 
      return total_size.load(std::memory_order_acquire);
    }
    
    size_t getCapacity() const { 
      size_t total_cap = 0;
      for (const auto& block : blocks) {
        total_cap += block->capacity;
      }
      return total_cap;
    }
    
    bool empty() const { 
      return getSize() == 0;
    }
    
    void clear() {
      // Clear all blocks
      for (auto& block : blocks) {
        block->used.store(0, std::memory_order_release);
      }
      total_size.store(0, std::memory_order_release);
    }
    
    // Reserve is a no-op now since we grow dynamically
    void reserve(size_t newCapacity) {
      // No-op - we grow dynamically
    }
    
    // Memory statistics
    size_t memoryUsage() const {
      return getCapacity() * sizeof(MinimalMapping) + sizeof(*this);
    }
    
    size_t actualMemoryUsage() const {
      return getSize() * sizeof(MinimalMapping) + sizeof(*this);
    }
    
    // Debug information
    void printStats() const {
      std::cerr << "LockFreeCompressedMappingStore stats:\n";
      std::cerr << "  Blocks: " << blocks.size() << "\n";
      std::cerr << "  Total capacity: " << getCapacity() << "\n";
      std::cerr << "  Total used: " << getSize() << "\n";
      std::cerr << "  Memory allocated: " << (getCapacity() * sizeof(MinimalMapping) / 1024 / 1024) << " MB\n";
      std::cerr << "  Memory used: " << (getSize() * sizeof(MinimalMapping) / 1024 / 1024) << " MB\n";
    }
    
    // Sort operations (only safe after parallel phase completes)
    void sortByQueryPos() {
      // First collect all mappings
      auto all_mappings = getCompressedMappings();
      
      // Sort them
      std::sort(all_mappings.begin(), all_mappings.end(),
                [](const MinimalMapping& a, const MinimalMapping& b) {
                  return a.query_pos < b.query_pos;
                });
      
      // Clear and re-add in sorted order
      clear();
      for (const auto& m : all_mappings) {
        // Convert back to MappingResult for adding
        MappingResult mr = expandMinimalMapping(m, querySeqId, queryLen);
        mr.sketchSize = defaultSketchSize.load();
        mr.kmerComplexity = defaultKmerComplexity.load();
        addMapping(mr);
      }
    }
    
    void sortByRefPos() {
      // First collect all mappings
      auto all_mappings = getCompressedMappings();
      
      // Sort them
      std::sort(all_mappings.begin(), all_mappings.end(),
                [](const MinimalMapping& a, const MinimalMapping& b) {
                  return std::tie(a.ref_seqId, a.ref_pos) < 
                         std::tie(b.ref_seqId, b.ref_pos);
                });
      
      // Clear and re-add in sorted order
      clear();
      for (const auto& m : all_mappings) {
        // Convert back to MappingResult for adding
        MappingResult mr = expandMinimalMapping(m, querySeqId, queryLen);
        mr.sketchSize = defaultSketchSize.load();
        mr.kmerComplexity = defaultKmerComplexity.load();
        addMapping(mr);
      }
    }
  };
}

#endif // LOCK_FREE_COMPRESSED_MAPPING_HPP