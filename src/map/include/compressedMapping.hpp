/**
 * @file    compressedMapping.hpp
 * @brief   Compressed mapping storage for memory efficiency
 * @author  Memory Optimization
 */

#ifndef COMPRESSED_MAPPING_HPP
#define COMPRESSED_MAPPING_HPP

#include <vector>
#include <algorithm>
#include <mutex>
#include "map/include/base_types.hpp"

namespace skch
{
  /**
   * @class CompressedMappingStore
   * @brief Stores mappings in compressed format to reduce memory usage
   *
   * This class provides a memory-efficient storage for mapping results
   * by using MinimalMapping internally while providing MappingResult
   * interface externally.
   */
  class CompressedMappingStore {
  private:
    std::vector<MinimalMapping> mappings;
    seqno_t querySeqId;
    offset_t queryLen;
    mutable std::mutex mutex;  // For thread safety
    
    // Additional metadata that's constant for all mappings in the store
    float defaultKmerComplexity = 0.0;
    int defaultSketchSize = 0;
    
  public:
    // Constructor
    CompressedMappingStore(seqno_t qSeqId = 0, offset_t qLen = 0) 
        : querySeqId(qSeqId), queryLen(qLen) {
      mappings.reserve(1000);  // Reserve initial capacity
    }
    
    // Add a mapping to the store
    void addMapping(const MappingResult& m) {
      std::lock_guard<std::mutex> lock(mutex);
      
      // Update metadata if not set
      if (defaultSketchSize == 0 && m.sketchSize > 0) {
        defaultSketchSize = m.sketchSize;
      }
      if (defaultKmerComplexity == 0.0 && m.kmerComplexity > 0.0) {
        defaultKmerComplexity = m.kmerComplexity;
      }
      
      // Compress and store
      mappings.push_back(compressMapping(m));
      
      // Update query info if needed
      if (querySeqId == 0) {
        querySeqId = m.querySeqId;
        queryLen = m.queryLen;
      }
    }
    
    // Add multiple mappings
    void addMappings(const std::vector<MappingResult>& results) {
      std::lock_guard<std::mutex> lock(mutex);
      mappings.reserve(mappings.size() + results.size());
      
      for (const auto& m : results) {
        if (defaultSketchSize == 0 && m.sketchSize > 0) {
          defaultSketchSize = m.sketchSize;
        }
        if (defaultKmerComplexity == 0.0 && m.kmerComplexity > 0.0) {
          defaultKmerComplexity = m.kmerComplexity;
        }
        mappings.push_back(compressMapping(m));
      }
      
      if (querySeqId == 0 && !results.empty()) {
        querySeqId = results[0].querySeqId;
        queryLen = results[0].queryLen;
      }
    }
    
    // Get a mapping by index
    MappingResult getMapping(size_t idx) const {
      std::lock_guard<std::mutex> lock(mutex);
      if (idx >= mappings.size()) {
        throw std::out_of_range("Index out of range in CompressedMappingStore");
      }
      
      MappingResult result = expandMinimalMapping(mappings[idx], querySeqId, queryLen);
      
      // Restore metadata
      result.sketchSize = defaultSketchSize;
      result.kmerComplexity = defaultKmerComplexity;
      
      return result;
    }
    
    // Get all mappings
    std::vector<MappingResult> getAllMappings() const {
      std::lock_guard<std::mutex> lock(mutex);
      std::vector<MappingResult> results;
      results.reserve(mappings.size());
      
      for (const auto& m : mappings) {
        MappingResult result = expandMinimalMapping(m, querySeqId, queryLen);
        result.sketchSize = defaultSketchSize;
        result.kmerComplexity = defaultKmerComplexity;
        results.push_back(result);
      }
      
      return results;
    }
    
    // Direct access to compressed mappings (for advanced operations)
    const std::vector<MinimalMapping>& getCompressedMappings() const {
      return mappings;
    }
    
    // Mark a mapping as discarded
    void markDiscard(size_t idx, bool discard = true) {
      std::lock_guard<std::mutex> lock(mutex);
      if (idx < mappings.size()) {
        mappings[idx].setDiscarded(discard);
      }
    }
    
    // Apply a function to each mapping
    template<typename Func>
    void forEach(Func&& func) {
      std::lock_guard<std::mutex> lock(mutex);
      for (size_t i = 0; i < mappings.size(); ++i) {
        MappingResult result = expandMinimalMapping(mappings[i], querySeqId, queryLen);
        result.sketchSize = defaultSketchSize;
        result.kmerComplexity = defaultKmerComplexity;
        func(result, i);
      }
    }
    
    // Filter mappings (remove discarded)
    void removeDiscarded() {
      std::lock_guard<std::mutex> lock(mutex);
      mappings.erase(
        std::remove_if(mappings.begin(), mappings.end(),
                      [](const MinimalMapping& m) { return m.isDiscarded(); }),
        mappings.end()
      );
    }
    
    // Size and capacity
    size_t size() const { 
      std::lock_guard<std::mutex> lock(mutex);
      return mappings.size(); 
    }
    
    bool empty() const { 
      std::lock_guard<std::mutex> lock(mutex);
      return mappings.empty(); 
    }
    
    void clear() {
      std::lock_guard<std::mutex> lock(mutex);
      mappings.clear();
    }
    
    void reserve(size_t n) {
      std::lock_guard<std::mutex> lock(mutex);
      mappings.reserve(n);
    }
    
    // Memory usage statistics
    size_t memoryUsage() const {
      std::lock_guard<std::mutex> lock(mutex);
      return mappings.size() * sizeof(MinimalMapping) + sizeof(*this);
    }
    
    size_t uncompressedMemoryUsage() const {
      std::lock_guard<std::mutex> lock(mutex);
      return mappings.size() * sizeof(MappingResult) + sizeof(*this);
    }
    
    float compressionRatio() const {
      if (mappings.empty()) return 1.0f;
      return static_cast<float>(uncompressedMemoryUsage()) / memoryUsage();
    }
    
    // Update query information
    void setQueryInfo(seqno_t qSeqId, offset_t qLen) {
      std::lock_guard<std::mutex> lock(mutex);
      querySeqId = qSeqId;
      queryLen = qLen;
    }
    
    // Sorting support
    void sortByQueryPos() {
      std::lock_guard<std::mutex> lock(mutex);
      std::sort(mappings.begin(), mappings.end(),
                [](const MinimalMapping& a, const MinimalMapping& b) {
                  return a.query_pos < b.query_pos;
                });
    }
    
    void sortByRefPos() {
      std::lock_guard<std::mutex> lock(mutex);
      std::sort(mappings.begin(), mappings.end(),
                [](const MinimalMapping& a, const MinimalMapping& b) {
                  return std::tie(a.ref_seqId, a.ref_pos) < 
                         std::tie(b.ref_seqId, b.ref_pos);
                });
    }
  };
}

#endif // COMPRESSED_MAPPING_HPP