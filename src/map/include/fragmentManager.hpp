/**
 * @file    fragmentManager.hpp
 * @brief   Fragment data structures and management for sequence mapping
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef FRAGMENT_MANAGER_HPP
#define FRAGMENT_MANAGER_HPP

#include <memory>
#include <vector>
#include <mutex>
#include <atomic>
#include <string>
#include "map/include/base_types.hpp"
#include "map/include/compressedMapping.hpp"
#include "common/progress.hpp"

namespace skch
{
  // Forward declaration
  struct QueryMappingOutput;

  /**
   * @brief Fragment data structure for processing sequence fragments
   */
  struct FragmentData {
      const char* seq;
      int len;
      int fullLen; 
      seqno_t seqId;
      std::string seqName;
      int refGroup;
      int fragmentIndex;
      std::shared_ptr<QueryMappingOutput> output;
      std::shared_ptr<std::atomic<int>> processedCounter;

      // Constructor
      FragmentData(const char* s, int l, int fl, seqno_t sid, 
                   const std::string& sn, int rg, int idx,
                   std::shared_ptr<QueryMappingOutput> out,
                   std::shared_ptr<std::atomic<int>> counter = nullptr)
          : seq(s), len(l), fullLen(fl), seqId(sid), seqName(sn),
            refGroup(rg), fragmentIndex(idx), output(out), processedCounter(counter) {}
  };

  /**
   * @brief Output structure for query mapping results
   */
  struct QueryMappingOutput {
      std::string queryName;
      CompressedMappingStore results;            // Non-merged mappings (compressed)
      CompressedMappingStore mergedResults;      // Maximally merged mappings (compressed)
      std::mutex mutex;
      progress_meter::ProgressMeter& progress;
      
      QueryMappingOutput(const std::string& name, const std::vector<MappingResult>& r, 
                        const std::vector<MappingResult>& mr, progress_meter::ProgressMeter& p)
          : queryName(name), results(), mergedResults(), progress(p) {
          // Convert existing vectors to compressed storage
          results.addMappings(r);
          mergedResults.addMappings(mr);
      }
  };

  /**
   * @brief Manages fragment lifetime and processing
   */
  class FragmentManager {
  private:
      std::vector<std::shared_ptr<FragmentData>> fragments;
      std::mutex fragments_mutex;

  public:
      void add_fragment(std::shared_ptr<FragmentData> fragment) {
          std::lock_guard<std::mutex> lock(fragments_mutex);
          fragments.push_back(fragment);
      }

      void clear() {
          std::lock_guard<std::mutex> lock(fragments_mutex);
          fragments.clear();
      }

      const std::vector<std::shared_ptr<FragmentData>>& get_fragments() const {
          return fragments;
      }
  };
}

#endif // FRAGMENT_MANAGER_HPP