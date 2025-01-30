/**
 * @file    computeMap.hpp
 * @brief   implements the sequence mapping logic
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_MAP_HPP
#define SKETCH_MAP_HPP

#include <iterator>
#include <limits>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <zlib.h>
#include <cassert>
#include <numeric>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <atomic>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <sstream>
#include "common/atomic_queue/atomic_queue.h"

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/slidingMap.hpp"
#include "map/include/MIIteratorL2.hpp"
#include "map/include/filter.hpp"

//External includes
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "map_stats.hpp"
#include "robin-hood-hashing/robin_hood.h"
// if we ever want to do the union-find chaining in parallel
//#include "common/dset64-gccAtomic.hpp"
// this is for single-threaded use, but is more portable
#include "common/dset64.hpp"
//#include "assert.hpp"
#include "gsl/gsl_randist.h"

namespace skch
{
  struct QueryMappingOutput {
      std::string queryName;
      std::vector<MappingResult> results;        // Non-merged mappings
      std::vector<MappingResult> mergedResults;  // Maximally merged mappings  
      std::mutex mutex;
      progress_meter::ProgressMeter& progress;
      QueryMappingOutput(const std::string& name, const std::vector<MappingResult>& r, 
                        const std::vector<MappingResult>& mr, progress_meter::ProgressMeter& p)
          : queryName(name), results(r), mergedResults(mr), progress(p) {}
  };

  struct FragmentData {
      const char* seq;
      int len;
      int fullLen;
      seqno_t seqId;
      std::string seqName;
      int refGroup;
      int fragmentIndex;
      QueryMappingOutput* output;
      std::atomic<int>* fragments_processed;
  };

  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages
   */
  class Map
  {
    private:
      //Type for Stage L1's predicted candidate location
      struct L1_candidateLocus_t
      {
        seqno_t seqId;                    //sequence id where read is mapped

        /* read could be mapped with its begin location
         * from [rangeStartPos, rangeEndPos]
         */
        offset_t rangeStartPos;
        offset_t rangeEndPos;
        int intersectionSize;
      };

      static constexpr auto L1_locus_intersection_cmp = [](L1_candidateLocus_t& a, L1_candidateLocus_t& b)
      {
        return a.intersectionSize < b.intersectionSize;
      };

      //Cache for commonly used values
      offset_t cached_segment_length;
      int cached_minimum_hits;

      //Type for Stage L2's predicted mapping coordinate within each L1 candidate  
      struct L2_mapLocus_t
      {
        seqno_t seqId;                    //sequence id where read is mapped
        offset_t meanOptimalPos;          //Among multiple consecutive optimal positions, save the avg.
        offset_t optimalStart;            //optimal start mapping position (begin iterator)
        offset_t optimalEnd;              //optimal end mapping position (end iterator)
        int sharedSketchSize;             //count of shared sketch elements
        strand_t strand;
      };

    private:

      //algorithm parameters
      skch::Parameters param;

      //reference sketch
      skch::Sketch* refSketch;

      // Sequence ID manager
      std::unique_ptr<SequenceIdManager> idManager;

      // Vectors to store query and target sequences
      std::vector<std::string> querySequenceNames;
      std::vector<std::string> targetSequenceNames;

      typedef Sketch::MIIter_t MIIter_t;

      //Custom function for post processing the results, by default does nothing
      typedef std::function< void(const MappingResult&) > PostProcessResultsFn_t;
      PostProcessResultsFn_t processMappingResults;

      //Container to store query sequence name and length
      //used only if one-to-one filtering is ON
      std::vector<ContigInfo> qmetadata;

      //Vector for sketch cutoffs. Position [i] indicates the minimum intersection size required
      //for an L1 candidate if the best intersection size is i;
      std::vector<int> sketchCutoffs; 

      // Sequence ID manager
      // Atomic queues for input and output
      typedef atomic_queue::AtomicQueue<InputSeqProgContainer*, 1024, nullptr, true, true, false, false> input_atomic_queue_t;
      typedef atomic_queue::AtomicQueue<QueryMappingOutput*, 1024, nullptr, true, true, false, false> merged_mappings_queue_t;
      typedef atomic_queue::AtomicQueue<std::pair<seqno_t, MappingResultsVector_t*>*, 1024> aggregate_atomic_queue_t;
      typedef atomic_queue::AtomicQueue<std::string*, 1024> writer_atomic_queue_t;
      typedef atomic_queue::AtomicQueue<QueryMappingOutput*, 1024, nullptr, true, true, false, false> query_output_atomic_queue_t;
      typedef atomic_queue::AtomicQueue<FragmentData*, 8192, nullptr, true, true, false, false> fragment_atomic_queue_t;
      
      // Track maximum chain ID seen across all subsets
      std::atomic<offset_t> maxChainIdSeen{0};


    void processFragment(FragmentData* fragment, 
                         std::vector<IntervalPoint>& intervalPoints,
                         std::vector<L1_candidateLocus_t>& l1Mappings,
                         MappingResultsVector_t& l2Mappings,
                         QueryMetaData<MinVec_Type>& Q) {
        intervalPoints.clear();
        l1Mappings.clear();
        l2Mappings.clear();

        Q.seq = const_cast<char*>(fragment->seq);
        Q.len = fragment->len;
        Q.fullLen = fragment->fullLen;
        Q.seqId = fragment->seqId;
        Q.seqName = fragment->seqName;
        Q.refGroup = fragment->refGroup;

        mapSingleQueryFrag(Q, intervalPoints, l1Mappings, l2Mappings);

        std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){
            e.queryLen = fragment->fullLen;
            e.queryStartPos = fragment->fragmentIndex * param.segLength;
            e.queryEndPos = e.queryStartPos + fragment->len;
        });

        {
            std::lock_guard<std::mutex> lock(fragment->output->mutex);
            fragment->output->results.insert(fragment->output->results.end(), l2Mappings.begin(), l2Mappings.end());
            // Initialize mergedResults with same mappings
            if (param.mergeMappings && param.split) {
                fragment->output->mergedResults.insert(fragment->output->mergedResults.end(), l2Mappings.begin(), l2Mappings.end());
            }
        }

        // Update progress after processing the fragment
        fragment->output->progress.increment(1);

        fragment->fragments_processed->fetch_add(1, std::memory_order_relaxed);
        delete fragment;
    }
      
    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] refSketch   reference sketch
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      Map(skch::Parameters p,
          PostProcessResultsFn_t f = nullptr) :
        param(p),
        processMappingResults(f),
        sketchCutoffs(std::min<double>(p.sketchSize, skch::fixed::ss_table_max) + 1, 1),
        idManager(std::make_unique<SequenceIdManager>(
            p.querySequences,
            p.refSequences,
            std::vector<std::string>{p.query_prefix},
            std::vector<std::string>{p.target_prefix},
            std::string(1, p.prefix_delim),
            p.query_list,
            p.target_list)),
        cached_segment_length(p.segLength),
        cached_minimum_hits(p.minimum_hits > 0 ? p.minimum_hits : Stat::estimateMinimumHitsRelaxed(p.sketchSize, p.kmerSize, p.percentageIdentity, skch::fixed::confidence_interval))
          {
              // Initialize sequence names right after creating idManager
              this->querySequenceNames = idManager->getQuerySequenceNames();
              this->targetSequenceNames = idManager->getTargetSequenceNames();

              // Calculate total target length
              uint64_t total_target_length = 0;
              size_t target_seq_count = targetSequenceNames.size();
              std::string target_prefix = param.target_prefix.empty() ? "none" : param.target_prefix;

              for (const auto& seqName : targetSequenceNames) {
                  seqno_t seqId = idManager->getSequenceId(seqName);
                  total_target_length += idManager->getSequenceLength(seqId);
              }

              // Calculate total query length
              uint64_t total_query_length = 0;
              for (const auto& seqName : querySequenceNames) {
                  total_query_length += idManager->getSequenceLength(idManager->getSequenceId(seqName));
              }

              // Count unique groups
              std::unordered_set<int> query_groups, target_groups;
              for (const auto& seqName : querySequenceNames) {
                  query_groups.insert(idManager->getRefGroup(idManager->getSequenceId(seqName)));
              }
              for (const auto& seqName : targetSequenceNames) {
                  target_groups.insert(idManager->getRefGroup(idManager->getSequenceId(seqName)));
              }

              // Calculate average sizes
              double avg_query_size_per_group = query_groups.size() ? (double)total_query_length / query_groups.size() : 0;
              double avg_target_size_per_group = target_groups.size() ? (double)total_target_length / target_groups.size() : 0;

              std::cerr << "[wfmash::mashmap] " 
                        << querySequenceNames.size() << " queries (" << total_query_length << "bp) in "
                        << query_groups.size() << " groups (≈" << std::fixed << std::setprecision(0) << avg_query_size_per_group << "bp/group)" << std::endl
                        << "[wfmash::mashmap] "
                        << target_seq_count << " targets (" << total_target_length << "bp) in "
                        << target_groups.size() << " groups (≈" << std::fixed << std::setprecision(0) << avg_target_size_per_group << "bp/group)" 
                        << std::endl;

              if (p.stage1_topANI_filter) {
                  this->setProbs();
              }
              this->mapQuery();
          }

      // Removed populateIdManager() function

      ~Map() = default;

      private:
      void buildMetadataFromIndex() {
          for (const auto& fileName : param.refSequences) {
              faidx_t* fai = fai_load(fileName.c_str());
              if (fai == nullptr) {
                  std::cerr << "Error: Failed to load FASTA index for file " << fileName << std::endl;
                  exit(1);
              }

              int nseq = faidx_nseq(fai);
              for (int i = 0; i < nseq; ++i) {
                  const char* seq_name = faidx_iseq(fai, i);
                  int seq_len = faidx_seq_len(fai, seq_name);
                  if (seq_len == -1) {
                      std::cerr << "Error: Failed to get length for sequence " << seq_name << std::endl;
                      continue;
                  }
                  // Metadata is now handled by idManager, no need to push_back here
              }

              fai_destroy(fai);
          }
      }

      public:

    private:

      void setProbs()
      {

        float deltaANI = param.ANIDiff;
        float min_p = 1 - param.ANIDiffConf;
        int ss = std::min<double>(param.sketchSize, skch::fixed::ss_table_max);

        // Cache hg pmf results
        std::vector<std::vector<double>> sketchProbs(
            ss + 1,
            std::vector<double>(ss + 1.0)
        );
        for (auto ci = 0; ci <= ss; ci++) 
        {
          for (double y = 0; y <= ci; y++) 
          {
            sketchProbs[ci][y] = gsl_ran_hypergeometric_pdf(y, ss, ss-ci, ci);
          }
        }
        
        // Return true iff Pr(ANI_i >= ANI_max - deltaANI) >= min_p
        const auto distDiff = [this, &sketchProbs, deltaANI, min_p, ss] (int cmax, int ci) {
          double prAboveCutoff = 0;
          for (double ymax = 0; ymax <= cmax; ymax++) {
            // Pr (Ymax = ymax)
            double pymax = sketchProbs[cmax][ymax];

            // yi_cutoff is minimum jaccard numerator required to be within deltaANI of ymax
            double yi_cutoff = deltaANI == 0 ? ymax : (std::floor(skch::Stat::md2j(
                skch::Stat::j2md(ymax / ss, param.kmerSize) + deltaANI, 
                param.kmerSize
            ) * ss));

            // Pr Y_i < yi_cutoff
            //std::cerr << "CMF " << yi_cutoff - 1 << " " << ss << " " << ss-ci << " " << ci << std::endl;
            double pi_acc = (yi_cutoff - 1) >= 0 ? gsl_cdf_hypergeometric_P(yi_cutoff-1, ss, ss-ci, ci) : 0;

            // Pr Y_i >= yi_cutoff
            pi_acc = 1-pi_acc;

            // Pr that mash score from cj leads to an ANI at least deltaJ less than the ANI from cmax
            prAboveCutoff += pymax * pi_acc;
            if (prAboveCutoff > min_p)
            {
              return true;
            }
          }
          return prAboveCutoff > min_p; 
        };

        // Helper vector for binary search
        std::vector<int> ss_range(ss+1);
        std::iota (ss_range.begin(), ss_range.end(), 0);

        for (auto cmax = 1; cmax <= ss; cmax++) 
        {
          // Binary search to find the lowest acceptable ci
          int ci = std::distance(
              ss_range.begin(),
              std::upper_bound(
                ss_range.begin(),
                ss_range.begin() + ss,
                false,
                [&distDiff, cmax] (bool val, int ci) {
                  return distDiff(cmax, ci);
                }
              )
          );
          sketchCutoffs[cmax] = ci;


          // For really high min_p values and some values of cmax, there are no values of
          // ci that satisfy the cutoff, so we just set to the max
          if (sketchCutoffs[cmax] == 0) {
            sketchCutoffs[cmax] = 1;
          }
        }
        //for (auto overlap = 1; overlap <= ss; overlap++) 
        //{
          //DEBUG_ASSERT(sketchCutoffs[overlap] <= overlap);
        //}
      }

      /**
       * @brief   parse over sequences in query file and map each on the reference
       */
      void reader_thread(input_atomic_queue_t& input_queue,
                         std::atomic<bool>& reader_done,
                         progress_meter::ProgressMeter& progress,
                         SequenceIdManager& idManager) {
          // Define allowed_query_names here
          std::unordered_set<std::string> allowed_query_names;
          if (!param.query_list.empty()) {
              std::ifstream filter_list(param.query_list);
              std::string name;
              while (getline(filter_list, name)) {
                  allowed_query_names.insert(name);
              }
          }

          if (!param.querySequences.empty()) {
              const auto& fileName = param.querySequences[0]; // Assume single query input file
              seqiter::for_each_seq_in_file(
                  fileName,
                  querySequenceNames,
                  [&](const std::string& seq_name, const std::string& seq) {
                      seqno_t seqId = idManager.getSequenceId(seq_name);
                      auto input = new InputSeqProgContainer(seq, seq_name, seqId, progress);
                      while (!input_queue.try_push(input)) {
                          std::this_thread::sleep_for(std::chrono::milliseconds(10));
                      }
                  });
          }
          reader_done.store(true);
      }

      void worker_thread(input_atomic_queue_t& input_queue,
                         fragment_atomic_queue_t& fragment_queue,
                         merged_mappings_queue_t& merged_queue,
                         progress_meter::ProgressMeter& progress,
                         std::atomic<bool>& reader_done,
                         std::atomic<bool>& workers_done) {
          while (true) {
              InputSeqProgContainer* input = nullptr;
              if (input_queue.try_pop(input)) {
                  auto output = mapModule(input, fragment_queue);
                  while (!merged_queue.try_push(output)) {
                      std::this_thread::sleep_for(std::chrono::milliseconds(10));
                  }
                  delete input;
              } else if (reader_done.load() && input_queue.was_empty()) {
                  break;
              } else {
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }
      }

      void writer_thread(query_output_atomic_queue_t& output_queue,
                         std::atomic<bool>& workers_done,
                         seqno_t& totalReadsMapped,
                         std::ofstream& outstrm,
                         progress_meter::ProgressMeter& progress,
                         MappingResultsVector_t& allReadMappings) {
          int wait_count = 0;
          while (true) {
              QueryMappingOutput* output = nullptr;
              if (output_queue.try_pop(output)) {
                  wait_count = 0;
                  if(output->results.size() > 0)
                      totalReadsMapped++;
                  if (param.filterMode == filter::ONETOONE) {
                      allReadMappings.insert(allReadMappings.end(), output->results.begin(), output->results.end());
                  } else {
                      reportReadMappings(output->results, output->queryName, outstrm);
                  }
                  delete output;
              } else {
                  if (workers_done.load() && output_queue.was_empty()) {
                      ++wait_count;
                  }
                  if (wait_count > 10) {
                      break;
                  }
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }
      }

      std::vector<std::vector<std::string>> createTargetSubsets(const std::vector<std::string>& targetSequenceNames) {
        std::vector<std::vector<std::string>> target_subsets;
        uint64_t current_subset_size = 0;
        std::vector<std::string> current_subset;

        for (const auto& seqName : targetSequenceNames) {
            seqno_t seqId = idManager->getSequenceId(seqName);
            offset_t seqLen = idManager->getSequenceLength(seqId);
            current_subset.push_back(seqName);
            current_subset_size += seqLen;

            if (current_subset_size >= param.index_by_size || &seqName == &targetSequenceNames.back()) {
                if (!current_subset.empty()) {
                    target_subsets.push_back(current_subset);
                }
                current_subset.clear();
                current_subset_size = 0;
            }
        }
        if (!current_subset.empty()) {
            target_subsets.push_back(current_subset);
        }
        return target_subsets;
      }

      void mapQuery()
      {
        std::cerr << "[wfmash::mashmap] L1 filtering parameters: cached_minimum_hits=" << cached_minimum_hits << std::endl;

        //Count of reads mapped by us
        //Some reads are dropped because of short length
        seqno_t totalReadsPickedForMapping = 0;
        seqno_t totalReadsMapped = 0;

        std::ofstream outstrm(param.outFileName);

        // Get sequence names from ID manager

        // Calculate total target length
        uint64_t total_target_length = 0;
        size_t target_seq_count = targetSequenceNames.size();
        std::string target_prefix = param.target_prefix.empty() ? "none" : param.target_prefix;

        for (const auto& seqName : targetSequenceNames) {
            seqno_t seqId = idManager->getSequenceId(seqName);
            total_target_length += idManager->getSequenceLength(seqId);
        }

        // Calculate total query length
        uint64_t total_query_length = 0;
        for (const auto& seqName : querySequenceNames) {
            total_query_length += idManager->getSequenceLength(idManager->getSequenceId(seqName));
        }


        // Initialize atomic queues and flags
        input_atomic_queue_t input_queue;
        merged_mappings_queue_t merged_queue;
        fragment_atomic_queue_t fragment_queue;
        writer_atomic_queue_t writer_queue;

        this->querySequenceNames = idManager->getQuerySequenceNames();
        this->targetSequenceNames = idManager->getTargetSequenceNames();

        // Count total fragments across all queries
        uint64_t total_fragments = 0;
        for (const auto& seqName : querySequenceNames) {
            auto len = idManager->getSequenceLength(idManager->getSequenceId(seqName));
            // Count full fragments
            int noOverlapFragmentCount = len / param.segLength;
            // Add final fragment if needed
            if (noOverlapFragmentCount >= 1 && len % param.segLength != 0) {
                ++noOverlapFragmentCount;
            }
            total_fragments += noOverlapFragmentCount;
        }

        std::vector<std::vector<std::string>> target_subsets = createTargetSubsets(targetSequenceNames);

        // Calculate and log subset statistics
        uint64_t total_subset_size = 0;
        for (const auto& subset : target_subsets) {
            for (const auto& seqName : subset) {
                seqno_t seqId = idManager->getSequenceId(seqName);
                total_subset_size += idManager->getSequenceLength(seqId);
            }
        }
        double avg_subset_size = target_subsets.size() ? (double)total_subset_size / target_subsets.size() : 0;
        std::cerr << "[wfmash::mashmap] Target subsets: " << target_subsets.size();
        if (param.index_by_size > 0) {
            std::cerr << ", target size: " << param.index_by_size << "bp";
        }
        std::cerr << ", average size: " << std::fixed << std::setprecision(0) << avg_subset_size << "bp" << std::endl;

        typedef std::vector<MappingResult> MappingResultsVector_t;
        std::unordered_map<seqno_t, MappingResultsVector_t> combinedMappings;

        // Build index for the current subset
        // Open the index file once
        std::ifstream indexStream;
        if (!param.indexFilename.empty() && !param.create_index_only) {
            indexStream.open(param.indexFilename.string(), std::ios::binary);
            if (!indexStream) {
                std::cerr << "Error: Unable to open index file: " << param.indexFilename << std::endl;
                exit(1);
            }
        }

        // For each subset of target sequences
        uint64_t subset_count = 0;
        std::cerr << "[wfmash::mashmap] Number of target subsets: " << target_subsets.size() << std::endl;
        for (const auto& target_subset : target_subsets) {
            if (target_subset.empty()) {
                continue;  // Skip empty subsets
            }
            // Calculate total length of sequences in this subset
            uint64_t subset_length = 0;
            for (const auto& seqName : target_subset) {
                seqno_t seqId = idManager->getSequenceId(seqName);
                subset_length += idManager->getSequenceLength(seqId);
            }            

            if (param.create_index_only) {
                // Save the index to a file
                std::cerr << "[wfmash::mashmap] Building and saving index for subset " << subset_count << " with " << target_subset.size() << " sequences" << std::endl;
                refSketch = new skch::Sketch(param, *idManager, target_subset);
                std::string indexFilename = param.indexFilename.string();
                bool append = (subset_count != 0); // Append if not the first subset
                refSketch->writeIndex(target_subset, indexFilename, append);
                std::cerr << "[wfmash::mashmap] Index created for subset " << subset_count 
                          << " and saved to " << indexFilename << std::endl;
            } else {
                if (!param.indexFilename.empty()) {
                    // Load index from file
                    std::cerr << "[wfmash::mashmap] Loading index for subset " << subset_count << " with " << target_subset.size() << " sequences" << std::endl;
                    refSketch = new skch::Sketch(param, *idManager, target_subset, &indexStream);
                } else {
                    std::cerr << "[wfmash::mashmap] Building index for subset " << subset_count << " with " << target_subset.size() 
                             << " sequences (" << subset_length << " bp)" << std::endl;
                    refSketch = new skch::Sketch(param, *idManager, target_subset);
                }
                std::atomic<bool> reader_done(false);
                std::atomic<bool> workers_done(false);
                std::atomic<bool> fragments_done(false);
                processSubset(subset_count, target_subsets.size(), total_fragments, input_queue, merged_queue, 
                              fragment_queue, reader_done, workers_done, fragments_done, combinedMappings);
            }

            // Clean up the current refSketch
            delete refSketch;
            refSketch = nullptr;
            ++subset_count;
        }

        if (indexStream.is_open()) {
            indexStream.close();
        }

        if (param.create_index_only) {
            std::cerr << "[wfmash::mashmap] All indices created successfully. Exiting." << std::endl;
            exit(0);
        }

        // Process combined mappings
        std::atomic<bool> processing_done(false);
        std::atomic<bool> workers_done(false);
        std::atomic<bool> output_done(false);
        aggregate_atomic_queue_t aggregate_queue;

        // Get total count of mappings
        uint64_t totalMappings = 0;
        for (const auto& [querySeqId, mappings] : combinedMappings) {
            totalMappings += mappings.size();
        }

        // Initialize progress logger
        progress_meter::ProgressMeter progress(
            totalMappings * 2,
            "[wfmash::mashmap] merging and filtering");

        // Start worker threads
        std::vector<std::thread> workers;
        for (int i = 0; i < param.threads; ++i) {
            workers.emplace_back(&Map::processCombinedMappingsThread, this, std::ref(aggregate_queue), std::ref(writer_queue), std::ref(processing_done), std::ref(progress));
        }

        // Start output thread
        std::thread output_thread(&Map::outputThread, this, std::ref(outstrm), std::ref(writer_queue), std::ref(processing_done), std::ref(workers_done), std::ref(output_done));

        // Enqueue tasks
        for (auto& [querySeqId, mappings] : combinedMappings) {
            auto* task = new std::pair<seqno_t, MappingResultsVector_t*>(querySeqId, &mappings);
            while (!aggregate_queue.try_push(task)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }

        // Signal that all tasks have been enqueued
        processing_done.store(true);

        // Wait for worker threads to finish
        for (auto& worker : workers) {
            worker.join();
        }

        workers_done.store(true);

        // Wait for output thread to finish
        while (!output_done.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        output_thread.join();

        progress.finish();

      }

      void processSubset(uint64_t subset_count, size_t total_subsets, uint64_t total_fragments,
                         input_atomic_queue_t& input_queue, merged_mappings_queue_t& merged_queue,
                         fragment_atomic_queue_t& fragment_queue, std::atomic<bool>& reader_done,
                         std::atomic<bool>& workers_done, std::atomic<bool>& fragments_done,
                         std::unordered_map<seqno_t, MappingResultsVector_t>& combinedMappings)
      {
          progress_meter::ProgressMeter progress(
              total_fragments,
              "[wfmash::mashmap] mapping ("
              + std::to_string(subset_count + 1) + "/" + std::to_string(total_subsets) + ")");

          // Create temporary storage for this subset's mappings
          std::unordered_map<seqno_t, MappingResultsVector_t> subsetMappings;

          // Launch reader thread
          std::thread reader([&]() {
              reader_thread(input_queue, reader_done, progress, *idManager);
          });

          std::vector<std::thread> fragment_workers;
          for (int i = 0; i < param.threads; ++i) {
              fragment_workers.emplace_back([&]() {
                  fragment_thread(fragment_queue, fragments_done);
              });
          }

          // Launch worker threads
          std::vector<std::thread> workers;
          for (int i = 0; i < param.threads; ++i) {
              workers.emplace_back([&]() {
                  worker_thread(input_queue, fragment_queue, merged_queue, progress, reader_done, workers_done);
              });
          }

          // Launch aggregator thread with subset storage
          std::thread aggregator([&]() {
              aggregator_thread(merged_queue, workers_done, subsetMappings);
          });

          // Wait for all threads to complete
          reader.join();

          for (auto& worker : workers) {
              worker.join();
          }
          workers_done.store(true);
          fragments_done.store(true);

          for (auto& worker : fragment_workers) {
              worker.join();
          }

          aggregator.join();

          // Filter mappings within this subset before merging with previous results
          for (auto& [querySeqId, mappings] : subsetMappings) {
              
              // Merge with existing mappings for this query
              if (combinedMappings.find(querySeqId) == combinedMappings.end()) {
                  combinedMappings[querySeqId] = std::move(mappings);
              } else {
                  combinedMappings[querySeqId].insert(
                      combinedMappings[querySeqId].end(),
                      std::make_move_iterator(mappings.begin()),
                      std::make_move_iterator(mappings.end())
                  );
              }
          }

          // Reset flags for next iteration
          reader_done.store(false);
          workers_done.store(false);
          fragments_done.store(false);

          progress.finish();
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings with fewer than the target number of merged base mappings
       * @param[in]   input   mappings
       * @return              void
       */
      void filterWeakMappings(MappingResultsVector_t &readMappings, int64_t min_count)
      {
          readMappings.erase(
              std::remove_if(readMappings.begin(),
                             readMappings.end(),
                             [&](MappingResult &e){
                                 return e.blockLength < param.block_length //e.queryLen > e.blockLength
                                     || e.n_merged < min_count;
                             }),
              readMappings.end());
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings whose identity and query/ref length don't agree
       * @param[in]   input   mappings
       * @return              void
       */
      void filterFalseHighIdentity(MappingResultsVector_t &readMappings)
      {
          readMappings.erase(
              std::remove_if(readMappings.begin(),
                             readMappings.end(),
                             [&](MappingResult &e){
                                 int64_t q_l = (int64_t)e.queryEndPos - (int64_t)e.queryStartPos;
                                 int64_t r_l = (int64_t)e.refEndPos - (int64_t)e.refStartPos;
                                 uint64_t delta = std::abs(r_l - q_l);
                                 double len_id_bound = (1.0 - (double)delta/(((double)q_l+r_l)/2));
                                 return len_id_bound < std::min(0.7, std::pow(param.percentageIdentity,3));
                             }),
              readMappings.end());
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings whose split ids aren't to be kept
       * @param[in]   input   mappings
       * @param[in]   input
       * @return              void
       */
      void filterFailedSubMappings(MappingResultsVector_t &readMappings,
                                   const robin_hood::unordered_set<offset_t>& kept_chains)
      {
          readMappings.erase(
              std::remove_if(readMappings.begin(),
                             readMappings.end(),
                             [&](MappingResult &e){
                                 return kept_chains.count(e.splitMappingId) == 0;
                             }),
              readMappings.end());
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings by hash value
       * @param[in]   input   mappings
       * @param[in]   input
       * @return              void
       */
      void sparsifyMappings(MappingResultsVector_t &readMappings)
      {
          if (param.sparsity_hash_threshold < std::numeric_limits<uint64_t>::max()) {
              readMappings.erase(
                  std::remove_if(readMappings.begin(),
                                 readMappings.end(),
                                 [&](MappingResult &e){
                                     return e.hash() > param.sparsity_hash_threshold;
                                 }),
                  readMappings.end());
          }
      }

      /**
       * @brief                    helper to main filtering function
       * @details                  filters mappings by group
       * @param[in]   input        unfiltered mappings
       * @param[in]   output       filtered mappings
       * @param[in]   n_mappings   num mappings per segment
       * @param[in]   filter_ref   use Filter::ref instead of Filter::query
       * @return                   void
       */
      void filterByGroup(
          MappingResultsVector_t &unfilteredMappings,
          MappingResultsVector_t &filteredMappings,
          int n_mappings,
          bool filter_ref,
          const SequenceIdManager& idManager,
          progress_meter::ProgressMeter& progress)
      {
        filteredMappings.reserve(unfilteredMappings.size());

        std::sort(unfilteredMappings.begin(), unfilteredMappings.end(), [](const auto& a, const auto& b) 
            { return std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos); });
        auto subrange_begin = unfilteredMappings.begin();
        auto subrange_end = unfilteredMappings.begin();
        if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) 
        {
          std::vector<skch::MappingResult> tmpMappings;
          while (subrange_end != unfilteredMappings.end())
          {
            if (param.skip_prefix)
            {
              int currGroup = idManager.getRefGroup(subrange_begin->refSeqId);
              subrange_end = std::find_if_not(subrange_begin, unfilteredMappings.end(), [this, currGroup, &idManager] (const auto& candidate) {
                  return currGroup == idManager.getRefGroup(candidate.refSeqId);
              });
            }
            else
            {
              subrange_end = unfilteredMappings.end();
            }
            tmpMappings.insert(
                tmpMappings.end(), 
                std::make_move_iterator(subrange_begin), 
                std::make_move_iterator(subrange_end));
            std::sort(tmpMappings.begin(), tmpMappings.end(), [](const auto& a, const auto& b) 
                { return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos) < std::tie(b.queryStartPos, b.refSeqId, b.refStartPos); });
            if (filter_ref)
            {
                skch::Filter::ref::filterMappings(tmpMappings, idManager, n_mappings, param.dropRand, param.overlap_threshold);
            }
            else
            {
                skch::Filter::query::filterMappings(tmpMappings, n_mappings, param.dropRand, param.overlap_threshold, progress);
            }
            filteredMappings.insert(
                filteredMappings.end(), 
                std::make_move_iterator(tmpMappings.begin()), 
                std::make_move_iterator(tmpMappings.end()));
            tmpMappings.clear();
            subrange_begin = subrange_end;
          }
        }
        //Sort the mappings by query (then reference) position
        std::sort(
            filteredMappings.begin(), filteredMappings.end(),
            [](const MappingResult &a, const MappingResult &b) {
                return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos, a.strand) 
                    < std::tie(b.queryStartPos, b.refSeqId, b.refStartPos, b.strand);
            });
      }


      /**
       * @brief               main mapping function given an input read
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @return              output object containing the mappings
       */
      QueryMappingOutput* mapModule(InputSeqProgContainer* input,
                                    fragment_atomic_queue_t& fragment_queue) {

        QueryMappingOutput* output = new QueryMappingOutput{input->name, {}, {}, input->progress};
        std::atomic<int> fragments_processed{0};
        bool split_mapping = true;
        int refGroup = this->idManager->getRefGroup(input->seqId);

        std::vector<FragmentData*> fragments;
        int noOverlapFragmentCount = input->len / param.segLength;

        for (int i = 0; i < noOverlapFragmentCount; i++) {
            auto fragment = new FragmentData{
                &(input->seq)[0u] + i * param.segLength,
                static_cast<int>(param.segLength),
                static_cast<int>(input->len),
                input->seqId,
                input->name,
                refGroup,
                i,
                output,
                &fragments_processed
            };
            fragments.push_back(fragment);
        }

        if (noOverlapFragmentCount >= 1 && input->len % param.segLength != 0) {
            auto fragment = new FragmentData{
                &(input->seq)[0u] + input->len - param.segLength,
                static_cast<int>(param.segLength),
                static_cast<int>(input->len),
                input->seqId,
                input->name,
                refGroup,
                noOverlapFragmentCount,
                output,
                &fragments_processed
            };
            fragments.push_back(fragment);
            noOverlapFragmentCount++;
        }

        for (auto& fragment : fragments) {
            while (!fragment_queue.try_push(fragment)) {
                //std::this_thread::yield(); // too fast
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }

        // Wait for all fragments to be processed
        while (fragments_processed.load(std::memory_order_relaxed) < noOverlapFragmentCount) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        mappingBoundarySanityCheck(input, output->results);
          
        // Filter and get both merged and non-merged mappings
        auto [nonMergedMappings, mergedMappings] = filterSubsetMappings(output->results, input->progress);
        output->results = std::move(nonMergedMappings);
        output->mergedResults = std::move(mergedMappings);

        return output;
      }

      void fragment_thread(fragment_atomic_queue_t& fragment_queue,
                           std::atomic<bool>& fragments_done) {
          std::vector<IntervalPoint> intervalPoints;
          std::vector<L1_candidateLocus_t> l1Mappings;
          MappingResultsVector_t l2Mappings;
          QueryMetaData<MinVec_Type> Q;

          while (!fragments_done.load()) {
              FragmentData* fragment = nullptr;
              if (fragment_queue.try_pop(fragment)) {
                  if (fragment) {
                      processFragment(fragment, intervalPoints, l1Mappings, l2Mappings, Q);
                  }
              } else {
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }
      }

      void processAggregatedMappings(const std::string& queryName, MappingResultsVector_t& mappings, progress_meter::ProgressMeter& progress) {

          // XXX we should fix this combined condition
          if (param.mergeMappings && param.split) {
              auto maximallyMergedMappings = mergeMappingsInRange(mappings, param.chain_gap, progress);
              filterMaximallyMerged(maximallyMergedMappings, param, progress);
              robin_hood::unordered_set<offset_t> kept_chains;
              for (auto &mapping : maximallyMergedMappings) {
                  kept_chains.insert(mapping.splitMappingId);
              }
              mappings.erase(
                  std::remove_if(mappings.begin(), mappings.end(),
                                 [&kept_chains](const MappingResult &mapping) {
                                     return !kept_chains.count(mapping.splitMappingId);
                                 }),
                  mappings.end());
          } else {
              filterNonMergedMappings(mappings, param, progress);
          }

          if (param.filterLengthMismatches) {
              filterFalseHighIdentity(mappings);
          }

          sparsifyMappings(mappings);

          // Apply group filtering aggregated across all targets
          if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
              MappingResultsVector_t filteredMappings;
              filterByGroup(mappings, filteredMappings, param.numMappingsForSegment - 1, param.filterMode == filter::ONETOONE, *idManager, progress);
              mappings = std::move(filteredMappings);
          }

          // Removed reportReadMappings call
      }

      void aggregator_thread(merged_mappings_queue_t& merged_queue,
                             std::atomic<bool>& workers_done,
                             std::unordered_map<seqno_t, MappingResultsVector_t>& combinedMappings) {
          while (true) {
              QueryMappingOutput* output = nullptr;
              if (merged_queue.try_pop(output)) {
                  seqno_t querySeqId = idManager->getSequenceId(output->queryName);
                  auto& mappings = param.mergeMappings && param.split ? output->mergedResults : output->results;
                  // Chain IDs are already compacted in mapModule
                  combinedMappings[querySeqId].insert(
                      combinedMappings[querySeqId].end(),
                      mappings.begin(),
                      mappings.end()
                  );
                  delete output;
              } else if (workers_done.load() && merged_queue.was_empty()) {
                  break;
              } else {
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }
      }

      /**
       * @brief                       routine to handle mapModule's output of mappings
       * @param[in] output            mapping output object
       * @param[in] allReadMappings   vector to store mappings of all reads (optional use depending on filter)
       * @param[in] totalReadsMapped  counter to track count of reads mapped
       * @param[in] outstrm           outstream stream object
       */
      template <typename Vec>
      void mapModuleHandleOutput(MapModuleOutput* output,
                                 Vec &allReadMappings,
                                 seqno_t &totalReadsMapped,
                                 std::ofstream &outstrm) {

          if(output->readMappings.size() > 0)
            totalReadsMapped++;

          if (param.filterMode == filter::ONETOONE)
          {
            //Save for another filtering round
            allReadMappings.insert(allReadMappings.end(), output->readMappings.begin(), output->readMappings.end());
          }
          else
          {
            //Report mapping
            reportReadMappings(output->readMappings, output->qseqName, outstrm);
          }

          delete output;
        }

      /**
       * @brief                       Filter non-merged mappings
       * @param[in/out] readMappings  Mappings computed by Mashmap
       * @param[in]     param         Algorithm parameters
       */
      void filterNonMergedMappings(MappingResultsVector_t &readMappings, const Parameters& param, progress_meter::ProgressMeter& progress)
      {
          if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
              MappingResultsVector_t filteredMappings;
              filterByGroup(readMappings, filteredMappings, param.numMappingsForSegment - 1, false, *idManager, progress);
              readMappings = std::move(filteredMappings);
          }
      }

      /**
       * @brief                   map the parsed query sequence (L1 and L2 mapping)
       * @param[in]   Q           metadata about query sequence
       * @param[in]   outstrm     outstream stream where mappings will be reported
       * @param[out]  l2Mappings  Mapping results in the L2 stage
       */
      template<typename Q_Info, typename IPVec, typename L1Vec, typename VecOut>
        void mapSingleQueryFrag(Q_Info &Q, IPVec& intervalPoints, L1Vec& l1Mappings, VecOut &l2Mappings)
        {
#ifdef ENABLE_TIME_PROFILE_L1_L2
          auto t0 = skch::Time::now();
#endif
          //L1 Mapping
          doL1Mapping(Q, intervalPoints, l1Mappings);
          if (l1Mappings.size() == 0) {
            return;
          }

#ifdef ENABLE_TIME_PROFILE_L1_L2
          std::chrono::duration<double> timeSpentL1 = skch::Time::now() - t0;
          auto t1 = skch::Time::now();
#endif

          auto l1_begin = l1Mappings.begin();
          auto l1_end = l1Mappings.begin();
          while (l1_end != l1Mappings.end())
          {
            if (param.skip_prefix)
            {
              int currGroup = this->idManager->getRefGroup(l1_begin->seqId);
              l1_end = std::find_if_not(l1_begin, l1Mappings.end(), [this, currGroup] (const auto& candidate) {
                  return currGroup == this->idManager->getRefGroup(candidate.seqId);
              });
            }
            else
            {
              l1_end = l1Mappings.end();
            }

            //Sort L1 windows based on intersection size if using hg filter
            if (param.stage1_topANI_filter)
            {
              std::make_heap(l1_begin, l1_end, L1_locus_intersection_cmp);
            }
            doL2Mapping(Q, l1_begin, l1_end, l2Mappings);

            // Set beginning of next range
            l1_begin = l1_end;
          }

          // Sort output mappings
          std::sort(l2Mappings.begin(), l2Mappings.end(), [](const auto& a, const auto& b) 
              { return std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos); });

          // Add chain information
          // All mappings in this batch form a chain
          int32_t chain_id = maxChainIdSeen.fetch_add(1, std::memory_order_relaxed);
          int32_t chain_length = l2Mappings.size();
          int32_t chain_pos = 1;
          for (auto& mapping : l2Mappings) {
              mapping.chain_id = chain_id;
              mapping.chain_length = chain_length;
              mapping.chain_pos = chain_pos++;
          }

#ifdef ENABLE_TIME_PROFILE_L1_L2
          {
            std::chrono::duration<double> timeSpentL2 = skch::Time::now() - t1;
            std::chrono::duration<double> timeSpentMappingFragment = skch::Time::now() - t0;

            std::cerr << Q.seqId << " " << Q.len
              << " " << timeSpentL1.count()
              << " " << timeSpentL2.count()
              << " " << timeSpentMappingFragment.count()
              << "\n";
          }
#endif
        }

      template <typename Q_Info>
        void getSeedHits(Q_Info &Q)
        {
          Q.minmerTableQuery.reserve(param.sketchSize + 1);
          CommonFunc::sketchSequence(Q.minmerTableQuery, Q.seq, Q.len, param.kmerSize, param.alphabetSize, param.sketchSize, Q.seqId);
          if(Q.minmerTableQuery.size() == 0) {
            Q.sketchSize = 0;
            return;
          }

#ifdef DEBUG
          int orig_len = Q.minmerTableQuery.size();
#endif
          const double max_hash_01 = (long double)(Q.minmerTableQuery.back().hash) / std::numeric_limits<hash_t>::max();
          Q.kmerComplexity = (double(Q.minmerTableQuery.size()) / max_hash_01) / ((Q.len - param.kmerSize + 1)*2);

          // Removed frequent kmer filtering

          Q.sketchSize = Q.minmerTableQuery.size();
#ifdef DEBUG
          std::cerr << "INFO, wfmash::mashmap, read id " << Q.seqId << ", minmer count = " << Q.minmerTableQuery.size() << ", bad minmers = " << orig_len - Q.sketchSize << "\n";
#endif
        } 


      /**
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is 
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an 
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details 
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename Vec>
        void getSeedIntervalPoints(Q_Info &Q, Vec& intervalPoints)
        {

#ifdef DEBUG
          std::cerr<< "INFO, wfmash::mashmap, read id " << Q.seqId << ", minmer count = " << Q.minmerTableQuery.size() << " " << Q.len << "\n";
#endif

          //For invalid query (example : just NNNs), we may be left with 0 sketch size
          //Ignore the query in this case
          if(Q.minmerTableQuery.size() == 0)
            return;

          // Priority queue for sorting interval points
          using IP_const_iterator = std::vector<IntervalPoint>::const_iterator;
          std::vector<boundPtr<IP_const_iterator>> pq;
          pq.reserve(Q.sketchSize);
          constexpr auto heap_cmp = [](const auto& a, const auto& b) {return b < a;};

          for(auto it = Q.minmerTableQuery.begin(); it != Q.minmerTableQuery.end(); it++)
          {
            //Check if hash value exists in the reference lookup index
            const auto seedFind = refSketch->minmerPosLookupIndex.find(it->hash);

            if(seedFind != refSketch->minmerPosLookupIndex.end())
            {
              pq.emplace_back(boundPtr<IP_const_iterator> {seedFind->second.cbegin(), seedFind->second.cend()});
            }
          }
          std::make_heap(pq.begin(), pq.end(), heap_cmp);

          while(!pq.empty())
          {
            const IP_const_iterator ip_it = pq.front().it;
            //const auto& ref = this->sketch_metadata[ip_it->seqId];
            const auto& ref_name = this->idManager->getSequenceName(ip_it->seqId);
            //const auto& ref_len = this->idManager.getSeqLen(ip_it->seqId);
            bool skip_mapping = false;
            int queryGroup = idManager->getRefGroup(Q.seqId);
            int targetGroup = idManager->getRefGroup(ip_it->seqId);

            if (param.skip_self && queryGroup == targetGroup) skip_mapping = true;
            if (param.skip_prefix && queryGroup == targetGroup) skip_mapping = true;
            if (param.lower_triangular && Q.seqId <= ip_it->seqId) skip_mapping = true;
    
            if (!skip_mapping) {
              intervalPoints.push_back(*ip_it);
            }
            std::pop_heap(pq.begin(), pq.end(), heap_cmp);
            pq.back().it++;
            if (pq.back().it >= pq.back().end) 
            {
              pq.pop_back();
            }
            else
            {
              std::push_heap(pq.begin(), pq.end(), heap_cmp);
            }
          }

#ifdef DEBUG
          std::cerr << "INFO, wfmash::mashmap, read id " << Q.seqId << ", Count of seed hits in the reference = " << intervalPoints.size() / 2 << "\n";
#endif
        }


      template <typename Q_Info, typename IP_iter, typename Vec2>
        void computeL1CandidateRegions(
            Q_Info &Q, 
            IP_iter ip_begin, 
            IP_iter ip_end, 
            int minimumHits, 
            Vec2 &l1Mappings)
        {
#ifdef DEBUG
          std::cerr << "INFO, skch::Map:computeL1CandidateRegions, read id " << Q.seqId << std::endl;
#endif

          int overlapCount = 0;
          int strandCount = 0;
          int bestIntersectionSize = 0;
          std::vector<L1_candidateLocus_t> localOpts;

          // Keep track of all minmer windows that intersect with [i, i+windowLen]
          int windowLen = std::max<offset_t>(0, Q.len - param.segLength);
          auto trailingIt = ip_begin;
          auto leadingIt = ip_begin;

          // Group together local sketch intersection maximums that are within clusterLen of eachother
          //
          // Since setting up the L2 window [i, j] requires aggregating minmer windows over
          // [i-segLength, i), we might as well group L2 windows together which are closer than
          // segLength  
          int clusterLen = param.segLength;

          // Used to keep track of how many minmer windows for a particular hash are currently "open"
          // Only necessary when windowLen != 0.
          std::unordered_map<hash_t, int> hash_to_freq;

          if (param.stage1_topANI_filter) {
            double maxJaccard = refSketch->hgNumerator / static_cast<double>(Q.sketchSize);
            double cutoff_j = Stat::md2j(1 - param.percentageIdentity + param.ANIDiff, param.kmerSize);
            int minIntersectionSize = std::max(
                static_cast<int>(cutoff_j * Q.sketchSize),
                minimumHits);

            while (leadingIt != ip_end)
            {
              // Catch the trailing iterator up to the leading iterator - windowLen
              while (
                  trailingIt != ip_end 
                  && ((trailingIt->seqId == leadingIt->seqId && trailingIt->pos <= leadingIt->pos - windowLen)
                    || trailingIt->seqId < leadingIt->seqId))
              {
                if (trailingIt->side == side::CLOSE) {
                  if (windowLen != 0)
                    hash_to_freq[trailingIt->hash]--;
                  if (windowLen == 0 || hash_to_freq[trailingIt->hash] == 0) {
                    overlapCount--;
                  }
                }
                trailingIt++;
              }
              auto currentPos = leadingIt->pos;
              while (leadingIt != ip_end && leadingIt->pos == currentPos) {
                if (leadingIt->side == side::OPEN) {
                  if (windowLen == 0 || hash_to_freq[leadingIt->hash] == 0) {
                    overlapCount++;
                  }
                  if (windowLen != 0)
                    hash_to_freq[leadingIt->hash]++;
                }
                leadingIt++;
              }

              //DEBUG_ASSERT(overlapCount >= 0, windowLen, trailingIt->seqId, trailingIt->pos, leadingIt->seqId, leadingIt->pos);
              //DEBUG_ASSERT(windowLen != 0 || overlapCount <= Q.sketchSize, windowLen, trailingIt->seqId, trailingIt->pos, leadingIt->seqId, leadingIt->pos);

              //Is this sliding window the best we have so far?
              bestIntersectionSize = std::max(bestIntersectionSize, overlapCount);
            }

            // Only go back through to find local opts if we know that there are some that are 
            // large enough
            if (bestIntersectionSize < minIntersectionSize) 
            {
              return;
            } else 
            {
              minimumHits = std::max(
                  sketchCutoffs[
                    int(std::min(bestIntersectionSize, Q.sketchSize) 
                      / std::max<double>(1, param.sketchSize / skch::fixed::ss_table_max))
                  ],
                  minIntersectionSize);
            }
          } 
          
          // Clear freq dict, as there will be left open CLOSE points at the end of the last seq
          // that we never got to
          hash_to_freq.clear();

          // Since there can be more than sketchSize windows that overlap w/ [i, i+windowLen]
          // cap the best intersection size 
          bestIntersectionSize = std::min(bestIntersectionSize, Q.sketchSize);

          bool in_candidate = false;
          L1_candidateLocus_t l1_out = {};
          trailingIt = ip_begin;
          leadingIt = ip_begin;

          // Keep track of 3 consecutive points so that we can track local optimums
          overlapCount = 0;
          int prevOverlap = 0;
          int prevPrevOverlap = 0;

          // Need to keep track of two positions, as the previous one will be the local optimum
          SeqCoord prevPos;
          SeqCoord currentPos{leadingIt->seqId, leadingIt->pos};


          while (leadingIt != ip_end)
          {
            prevPrevOverlap = prevOverlap;
            prevOverlap = overlapCount;

            //TODO LEADING it should only hit opens
            // We should only iterate through a new window when we come across an OPEN,
            // right now, this basically happens since every CLOSE should be an OPEN.
            // This doesn't invalidate the logic, just potentially wastes time
            while (
                trailingIt != ip_end 
                && ((trailingIt->seqId == leadingIt->seqId && trailingIt->pos <= leadingIt->pos - windowLen)
                  || trailingIt->seqId < leadingIt->seqId))
            {
              if (trailingIt->side == side::CLOSE) {
                if (windowLen != 0)
                  hash_to_freq[trailingIt->hash]--;
                if (windowLen == 0 || hash_to_freq[trailingIt->hash] == 0) {
                  overlapCount--;
                }
              }
              trailingIt++;
            }
            if (leadingIt->pos != currentPos.pos) {
              prevPos = currentPos;
              currentPos = SeqCoord{leadingIt->seqId, leadingIt->pos};
            }
            while (leadingIt != ip_end && leadingIt->pos == currentPos.pos) 
            {
              if (leadingIt->side == side::OPEN) {
                if (windowLen == 0 || hash_to_freq[leadingIt->hash] == 0) {
                  overlapCount++;
                }
                if (windowLen != 0)
                  hash_to_freq[leadingIt->hash]++;
              }
              leadingIt++;
            }
          if ( prevOverlap >= minimumHits
              //&& prevOverlap > overlapCount && prevOverlap >= prevPrevOverlap)
          ) {
            if (l1_out.seqId != prevPos.seqId && in_candidate) {
              localOpts.push_back(l1_out);
              l1_out = {};
              in_candidate = false;
            }
            if (!in_candidate) {
              l1_out.rangeStartPos = prevPos.pos - windowLen;
              l1_out.rangeEndPos = prevPos.pos - windowLen;
              l1_out.seqId = prevPos.seqId;
              l1_out.intersectionSize = prevOverlap;
              in_candidate = true;
            } else {
              if (param.stage2_full_scan) {
                l1_out.intersectionSize = std::max(l1_out.intersectionSize, prevOverlap);
                l1_out.rangeEndPos = prevPos.pos - windowLen;
              }
              else if (l1_out.intersectionSize < prevOverlap) {
                l1_out.intersectionSize = prevOverlap;
                l1_out.rangeStartPos = prevPos.pos - windowLen;
                l1_out.rangeEndPos = prevPos.pos - windowLen;
              }
            }
          } 
          else {
            if (in_candidate) {
              localOpts.push_back(l1_out);
              l1_out = {};
            }
            in_candidate = false;
          }
        }
        if (in_candidate) {
          localOpts.push_back(l1_out);
        }
        

        // Join together proximal local opts
        for (auto& l1_out : localOpts) 
        {
          if (l1Mappings.empty() 
              || l1_out.seqId != l1Mappings.back().seqId 
              || l1_out.rangeStartPos > l1Mappings.back().rangeEndPos + clusterLen) 
          {
            l1Mappings.push_back(l1_out); 
          } 
          else 
          {
            l1Mappings.back().rangeEndPos = l1_out.rangeEndPos;
            l1Mappings.back().intersectionSize = std::max(l1_out.intersectionSize, l1Mappings.back().intersectionSize);
          }
        }
      }


      /**
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename IPVec, typename L1Vec>
        void doL1Mapping(Q_Info &Q, IPVec& intervalPoints, L1Vec& l1Mappings)
        {
          //1. Compute the minmers
          getSeedHits(Q);

          //Catch all NNNNNN case
          if (Q.sketchSize == 0 || Q.kmerComplexity < param.kmerComplexityThreshold) {
            return;
          }

          //2. Compute windows and sort
          getSeedIntervalPoints(Q, intervalPoints);

          //3. Compute L1 windows
          // Always respect the minimum hits parameter if set
          int minimumHits = param.minimum_hits > 0 ? 
              param.minimum_hits : 
              (Q.len == cached_segment_length ? 
                  cached_minimum_hits : 
                  Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity, skch::fixed::confidence_interval));

          // For each "group"
          auto ip_begin = intervalPoints.begin();
          auto ip_end = intervalPoints.begin();
          while (ip_end != intervalPoints.end())
          {
            if (param.skip_prefix)
            {
              int currGroup = this->idManager->getRefGroup(ip_begin->seqId);
              ip_end = std::find_if_not(ip_begin, intervalPoints.end(), [this, currGroup] (const auto& ip) {
                  return currGroup == this->idManager->getRefGroup(ip.seqId);
              });
            }
            else
            {
              ip_end = intervalPoints.end();
            }
            computeL1CandidateRegions(Q, ip_begin, ip_end, minimumHits, l1Mappings);

            ip_begin = ip_end;
          }
        }


      /**
       * @brief     Build metadata and reference groups for sequences
       * @details   Read FAI files, sort sequences, and assign groups
       */
      void buildRefGroups() {
          std::vector<std::tuple<std::string, size_t, offset_t>> seqInfoWithIndex;
          size_t totalSeqs = 0;

          for (const auto& fileName : param.refSequences) {
              std::string faiName = fileName + ".fai";
              std::ifstream faiFile(faiName);
              
              if (!faiFile.is_open()) {
                  std::cerr << "Error: Unable to open FAI file: " << faiName << std::endl;
                  exit(1);
              }

              std::string line;
              while (std::getline(faiFile, line)) {
                  std::istringstream iss(line);
                  std::string seqName;
                  offset_t seqLength;
                  iss >> seqName >> seqLength;

                  seqInfoWithIndex.emplace_back(seqName, totalSeqs++, seqLength);
              }
          }

          std::sort(seqInfoWithIndex.begin(), seqInfoWithIndex.end());

          std::vector<int> refGroups(totalSeqs);
          // Removed as sketch_metadata is no longer used
          int currentGroup = 0;
          std::string prevPrefix;

          for (const auto& [seqName, originalIndex, seqLength] : seqInfoWithIndex) {
              std::string currPrefix = seqName.substr(0, seqName.find_last_of(param.prefix_delim));
              
              if (currPrefix != prevPrefix) {
                  currentGroup++;
                  prevPrefix = currPrefix;
              }

              refGroups[originalIndex] = currentGroup;
              // Metadata is now handled by idManager, no need to push_back here
          }

          // Removed refIdGroup swap as it's no longer needed

          if (totalSeqs == 0) {
              std::cerr << "[wfmash::mashmap] ERROR: No sequences indexed!" << std::endl;
              exit(1);
          }
      }

      /**
       * @brief                                 Revise L1 candidate regions to more precise locations
       * @param[in]   Q                         query sequence information
       * @param[in]   l1Mappings                candidate regions for query sequence found at L1
       * @param[out]  l2Mappings                Mapping results in the L2 stage
       */
      template <typename Q_Info, typename L1_Iter, typename VecOut>
        void doL2Mapping(Q_Info &Q, L1_Iter l1_begin, L1_Iter l1_end, VecOut &l2Mappings)
        {
          ///2. Walk the read over the candidate regions and compute the jaccard similarity with minimum s sketches
          std::vector<L2_mapLocus_t> l2_vec;
          double bestJaccardNumerator = 0;
          auto loc_iterator = l1_begin;
          while (loc_iterator != l1_end)
          {
            L1_candidateLocus_t& candidateLocus = *loc_iterator;

            if (param.stage1_topANI_filter)
            {
              // Use the global Jaccard numerator here
              double jaccardSimilarity = refSketch->hgNumerator / Q.sketchSize;
              double mash_dist = Stat::j2md(jaccardSimilarity, param.kmerSize);
              double cutoff_ani = std::max(0.0, (1 - mash_dist) - param.ANIDiff);
              double cutoff_j = Stat::md2j(1 - cutoff_ani, param.kmerSize);

              double candidateJaccard = static_cast<double>(candidateLocus.intersectionSize) / Q.sketchSize;

              if (candidateJaccard < cutoff_j) 
              {
                break;
              }
            }

            l2_vec.clear();
            computeL2MappedRegions(Q, candidateLocus, l2_vec);

            for (auto& l2 : l2_vec) 
            {
              //Compute mash distance using calculated jaccard
              float mash_dist = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);

              float nucIdentity = (1 - mash_dist);
              //float nucIdentityUpperBound = getANIUBfromJaccardNum(Q.sketchSize, l2.sharedSketchSize);
              float nucIdentityUpperBound = 1 - Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, skch::fixed::confidence_interval);

              //Report the alignment if it passes our identity threshold and,
              // if we are in all-vs-all mode, it isn't a self-mapping,
              // and if we are self-mapping, the query is shorter than the target
              const auto& ref = this->idManager->getContigInfo(l2.seqId);
              if((param.keep_low_pct_id && nucIdentityUpperBound >= param.percentageIdentity)
                  || nucIdentity >= param.percentageIdentity)
              {
                //Track the best jaccard numerator
                bestJaccardNumerator = std::max<double>(bestJaccardNumerator, l2.sharedSketchSize);

                MappingResult res;

                //Save the output
                {
                  res.queryLen = Q.len;
                  res.refStartPos = l2.meanOptimalPos;
                  res.refEndPos = l2.meanOptimalPos + Q.len;
                  res.queryStartPos = 0;
                  res.queryEndPos = Q.len;
                  res.refSeqId = l2.seqId;
                  res.querySeqId = Q.seqId;
                  res.nucIdentity = nucIdentity;
                  res.nucIdentityUpperBound = nucIdentityUpperBound;
                  res.sketchSize = Q.sketchSize;
                  res.conservedSketches = l2.sharedSketchSize;
                  res.blockLength = std::max(res.refEndPos - res.refStartPos, res.queryEndPos - res.queryStartPos);
                  res.approxMatches = std::round(res.nucIdentity * res.blockLength / 100.0);
                  res.strand = l2.strand; 
                  res.kmerComplexity = Q.kmerComplexity;

                  res.selfMapFilter = ((param.skip_self || param.skip_prefix) && Q.fullLen > ref.len);

                } 
                l2Mappings.push_back(res);
              }
            }

            if (param.stage1_topANI_filter) 
            {
              std::pop_heap(l1_begin, l1_end, L1_locus_intersection_cmp); 
              l1_end--; //"Pop back" 
            }
            else 
            {
              loc_iterator++;
            }
          }
          //std::cerr << "For an segment with " << l1Mappings.size()
            //<< " L1 mappings "
            //<< " there were " << l2Mappings.size() << " L2 mappings\n";
        }

      /**
       * @brief                                 Find optimal mapping within an L1 candidate
       * @param[in]   Q                         query sequence information
       * @param[in]   candidateLocus            L1 candidate location
       * @param[out]  l2_out                    L2 mapping inside L1 candidate
       */
      template <typename Q_Info, typename Vec>
        void computeL2MappedRegions(Q_Info &Q,
            L1_candidateLocus_t &candidateLocus,
            Vec &l2_vec_out)
        {
#ifdef DEBUG
          //std::cerr << "INFO, skch::Map:computeL2MappedRegions, read id " << Q.seqName << "_" << Q.startPos << std::endl; 
#endif
           
          auto& minmerIndex = refSketch->minmerIndex;

          //candidateLocus.rangeStartPos -= param.segLength;
          //candidateLocus.rangeEndPos += param.segLength;
          
          // Get first potential mashimizer
          const MinmerInfo first_minmer = MinmerInfo {0, candidateLocus.rangeStartPos - param.segLength - 1, 0, candidateLocus.seqId, 0};

          //const MinmerInfo first_minmer = MinmerInfo {0, candidateLocus.seqId, -1, 0, 0};
          auto firstOpenIt = std::lower_bound(minmerIndex.begin(), minmerIndex.end(), first_minmer); 

          // Keeps track of the lowest end position
          std::vector<skch::MinmerInfo> slidingWindow;
          slidingWindow.reserve(Q.sketchSize);

          // Used to make a min-heap
          constexpr auto heap_cmp = [](const skch::MinmerInfo& l, const skch::MinmerInfo& r) {return l.wpos_end > r.wpos_end;};

          // windowIt keeps track of the end of window
          auto windowIt = firstOpenIt;

          // Keep track of all minmer windows that intersect with [i, i+windowLen]
          int windowLen = std::max<offset_t>(0, Q.len - param.segLength);

          // Used to keep track of how many minmer windows for a particular hash are currently "open"
          // Only necessary when windowLen != 0.
          std::unordered_map<hash_t, int> hash_to_freq;
          
          // slideMap tracks the S(A or B) and S(A) and S(B)
          SlideMapper<Q_Info> slideMap(Q);

          offset_t beginOptimalPos = 0;
          offset_t lastOptimalPos = 0;
          int bestSketchSize = 1;
          int bestIntersectionSize = 0;
          bool in_candidate = false;
          L2_mapLocus_t l2_out = {};

          // Set up the window
          while (windowIt != minmerIndex.end() && windowIt->seqId == candidateLocus.seqId && windowIt->wpos < candidateLocus.rangeStartPos) 
          {
            if (windowIt->wpos_end > candidateLocus.rangeStartPos) 
            {
              if (windowLen > 0) 
              {
                hash_to_freq[windowIt->hash]++;
              }
              if (windowLen == 0 || hash_to_freq[windowIt->hash] == 1) {
                slidingWindow.push_back(*windowIt);
                std::push_heap(slidingWindow.begin(), slidingWindow.end(), heap_cmp);
                slideMap.insert_minmer(*windowIt);
              }
            }
            windowIt++;
          }

          while (windowIt != minmerIndex.end() && windowIt->seqId == candidateLocus.seqId && windowIt->wpos <= candidateLocus.rangeEndPos + windowLen) 
          {
            int prev_strand_votes = slideMap.strand_votes;
            bool inserted = false;
            while (!slidingWindow.empty() && slidingWindow.front().wpos_end <= windowIt->wpos - windowLen) {

              // Remove minmer from end-ordered heap
              if (windowLen > 0) 
              {
                hash_to_freq[slidingWindow.front().hash]--;
              }
              if (windowLen == 0 || hash_to_freq[slidingWindow.front().hash] == 0) {
                // Remove minmer from  sorted window
                slideMap.delete_minmer(slidingWindow.front());
                std::pop_heap(slidingWindow.begin(), slidingWindow.end(), heap_cmp);
                slidingWindow.pop_back();
              }

            }
            inserted = true;
            if (windowLen > 0) 
            {
              hash_to_freq[windowIt->hash]++;
            }
            if (windowLen == 0 || hash_to_freq[windowIt->hash] == 1) {
              slideMap.insert_minmer(*windowIt);
              slidingWindow.push_back(*windowIt);
              std::push_heap(slidingWindow.begin(), slidingWindow.end(), heap_cmp);
            } else {
              windowIt++;
              continue;
            }

            bestIntersectionSize = std::max(bestIntersectionSize, slideMap.intersectionSize);

            //Is this sliding window the best we have so far?
            if (slideMap.sharedSketchElements > bestSketchSize)
            {
              // Get rid of all candidates seen so far
              l2_vec_out.clear();

              in_candidate = true;
              bestSketchSize = slideMap.sharedSketchElements;
              l2_out.sharedSketchSize = slideMap.sharedSketchElements;

              //Save the position
              l2_out.optimalStart = windowIt->wpos - windowLen;
              l2_out.optimalEnd = windowIt->wpos - windowLen;
            }
            else if(slideMap.sharedSketchElements == bestSketchSize)
            {
              if (!in_candidate) {
                l2_out.sharedSketchSize = slideMap.sharedSketchElements;

                //Save the position
                l2_out.optimalStart = windowIt->wpos - windowLen;
              }

              in_candidate = true;
              //Still save the position
              l2_out.optimalEnd = windowIt->wpos - windowLen;
            } else {
              if (in_candidate) {
                // Save and reset
                l2_out.meanOptimalPos =  (l2_out.optimalStart + l2_out.optimalEnd) / 2;
                l2_out.seqId = windowIt->seqId;
                l2_out.strand = prev_strand_votes >= 0 ? strnd::FWD : strnd::REV;
                if (l2_vec_out.empty() 
                    || l2_vec_out.back().optimalEnd + param.segLength < l2_out.optimalStart)
                {
                  l2_vec_out.push_back(l2_out);
                }
                else 
                {
                  l2_vec_out.back().optimalEnd = l2_out.optimalEnd;
                  l2_vec_out.back().meanOptimalPos = (l2_vec_out.back().optimalStart + l2_vec_out.back().optimalEnd) / 2;
                }
                l2_out = L2_mapLocus_t();
              }
              in_candidate = false;
            }
            if (inserted) {
              windowIt++;
            }
          }
          if (in_candidate) {
            // Save and reset
            l2_out.meanOptimalPos =  (l2_out.optimalStart + l2_out.optimalEnd) / 2;
            l2_out.seqId = std::prev(windowIt)->seqId;
            l2_out.strand = slideMap.strand_votes >= 0 ? strnd::FWD : strnd::REV;
            if (l2_vec_out.empty() 
                || l2_vec_out.back().optimalEnd + param.segLength < l2_out.optimalStart)
            {
              l2_vec_out.push_back(l2_out);
            }
            else 
            {
              l2_vec_out.back().optimalEnd = l2_out.optimalEnd;
              l2_vec_out.back().meanOptimalPos = (l2_vec_out.back().optimalStart + l2_vec_out.back().optimalEnd) / 2;
            }
          }
        }


      /**
       * @brief                       Merge the consecutive fragment mappings reported in each query
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void expandMappings(VecIn &readMappings, int expansion)
        {
            for (auto& m : readMappings) {
                m.refStartPos -= expansion;
                m.refEndPos += expansion;
                m.queryStartPos -= expansion;
                m.queryEndPos += expansion;
            }
        }

      void processMappingFragment(std::vector<MappingResult>::iterator start, std::vector<MappingResult>::iterator end) {
          auto& fragment = *start; // We'll update the first mapping in the fragment

          //std::cerr << "Processing fragment" << std::endl;
          // Compute fragment information
          for (auto it = start; it != end; ++it) {
              //std::cerr << "Merging " << it->queryStartPos << " " << it->queryEndPos << " " 
              //          << it->refStartPos << " " << it->refEndPos << " " << it->nucIdentity
              //          << " " << (it->strand == strnd::FWD ? "+" : "-") << std::endl;
        
              fragment.queryStartPos = std::min(fragment.queryStartPos, it->queryStartPos);
              fragment.refStartPos = std::min(fragment.refStartPos, it->refStartPos);
              fragment.queryEndPos = std::max(fragment.queryEndPos, it->queryEndPos);
              fragment.refEndPos = std::max(fragment.refEndPos, it->refEndPos);
          }

          fragment.blockLength = std::max(fragment.refEndPos - fragment.refStartPos, fragment.queryEndPos - fragment.queryStartPos);
          //fragment.blockLength = fragment.queryEndPos - fragment.queryStartPos;
                                          
          fragment.approxMatches = std::round(fragment.nucIdentity * fragment.blockLength / 100.0);

          fragment.n_merged = std::distance(start, end);

          // Calculate mean nucleotide identity
          fragment.nucIdentity = std::accumulate(start, end, 0.0,
                                                 [](double sum, const MappingResult& e) { return sum + e.nucIdentity; }
              ) / fragment.n_merged;

          // Calculate mean kmer complexity
          fragment.kmerComplexity = std::accumulate(start, end, 0.0,
                                                    [](double sum, const MappingResult& e) { return sum + e.kmerComplexity; }
              ) / fragment.n_merged;

          // Mark other mappings in this fragment for discard
          std::for_each(std::next(start), end, [](MappingResult& e) { e.discard = 1; });
      }

      void adjustConsecutiveMappings(std::vector<MappingResult>::iterator begin_maping,
                                     std::vector<MappingResult>::iterator end_mapping,
                                     const int threshold) {

          if (std::distance(begin_maping, end_mapping) < 2) return;

          //for (size_t i = 1; i < mappings.size(); ++i) {
          for (auto it = std::next(begin_maping); it != end_mapping; ++it) {
              auto& prev = *std::prev(it);
              auto& curr = *it;

              // Check if mappings are on the same reference sequence
              if (prev.refSeqId != curr.refSeqId
                  || prev.strand != curr.strand) continue;

              // Calculate gaps
              int query_gap = curr.queryStartPos - prev.queryEndPos;
              int ref_gap = curr.refStartPos - prev.refEndPos;

              // Check if both gaps are >0 and within the threshold
              if (query_gap > 0 && ref_gap > 0 && query_gap <= threshold && ref_gap <= threshold) {
                  // Calculate midpoints
                  int query_mid = (prev.queryEndPos + curr.queryStartPos) / 2;
                  int ref_mid = (prev.refEndPos + curr.refStartPos) / 2;

                  // Adjust the mappings
                  prev.queryEndPos = query_mid;
                  prev.refEndPos = ref_mid;
                  curr.queryStartPos = query_mid;
                  curr.refStartPos = ref_mid;

                  // Update block lengths
                  prev.blockLength = std::max(prev.refEndPos - prev.refStartPos, 
                                              prev.queryEndPos - prev.queryStartPos);
                  curr.blockLength = std::max(curr.refEndPos - curr.refStartPos, 
                                              curr.queryEndPos - curr.queryStartPos);

                  // Update approximate matches
                  prev.approxMatches = std::round(prev.nucIdentity * prev.blockLength / 100.0);
                  curr.approxMatches = std::round(curr.nucIdentity * curr.blockLength / 100.0);
              }
          }
      }

      double axis_weighted_euclidean_distance(int64_t dx, int64_t dy, double w = 0.5) {
          double euclidean = std::sqrt(dx*dx + dy*dy);
          double axis_factor = 1.0 - (2.0 * std::min(std::abs(dx), std::abs(dy))) / (std::abs(dx) + std::abs(dy));
          return euclidean * (1.0 + w * axis_factor);
      }

      /**
       * @brief                       Filter maximally merged mappings
       * @param[in]   readMappings    Merged mappings computed by Mashmap
       * @param[in]   param           Algorithm parameters
       * @return                      Filtered mappings
       */
      void filterMaximallyMerged(MappingResultsVector_t& readMappings, const Parameters& param, progress_meter::ProgressMeter& progress)
      {
          // Filter weak mappings
          filterWeakMappings(readMappings, std::floor(param.block_length / param.segLength));

          // Apply group filtering if necessary
          if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
              MappingResultsVector_t groupFilteredMappings;
              filterByGroup(readMappings, groupFilteredMappings, param.numMappingsForSegment - 1, false, *idManager, progress);
              readMappings = std::move(groupFilteredMappings);
          }
      }

      /**
       * @brief                       Merge fragment mappings by convolution of a 2D range over the alignment matrix
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       * @param[in]     max_dist      Distance to look in target and query
       */
      template <typename VecIn>
      VecIn mergeMappingsInRange(VecIn &readMappings,
                                 int max_dist,
                                 progress_meter::ProgressMeter& progress) {
          if (!param.split || readMappings.size() < 2) return readMappings;

          //Sort the mappings by query position, then reference sequence id, then reference position
          std::sort(
              readMappings.begin(), readMappings.end(),
              [](const MappingResult &a, const MappingResult &b) {
                  return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos, a.strand)
                      < std::tie(b.queryStartPos, b.refSeqId, b.refStartPos, b.strand);
              });

          //First assign a unique id to each split mapping in the sorted order
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              it->splitMappingId = std::distance(readMappings.begin(), it);
              it->discard = 0;
              it->chainPairScore = std::numeric_limits<double>::max();
              it->chainPairId = std::numeric_limits<int64_t>::min();
          }

          // set up our union find data structure to track merges
          std::vector<dsets::DisjointSets::Aint> ufv(readMappings.size());
          // this initializes everything
          auto disjoint_sets = dsets::DisjointSets(ufv.data(), ufv.size());

          //Start the procedure to identify the chains
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              double best_score = std::numeric_limits<double>::max();
              auto best_it2 = readMappings.end();
              // we we merge only with the best-scored previous mapping in query space
              if (it->chainPairScore != std::numeric_limits<double>::max()) {
                  disjoint_sets.unite(it->splitMappingId, it->chainPairId);
              }
              for (auto it2 = std::next(it); it2 != readMappings.end(); it2++) {
                  // If this mapping is for a different reference sequence, ignore
                  if (it2->refSeqId != it->refSeqId) {
                      continue;
                  }
                  //If this mapping is for exactly the same segment, ignore
                  if (it2->queryStartPos == it->queryStartPos) {
                      continue;
                  }
                  //If this mapping is too far from current mapping being evaluated in the query, stop finding a merge
                  if (it2->queryStartPos > it->queryEndPos + max_dist) {
                      break;
                  }
                  //If the next mapping is within range, we can potentially merge
                  if (it2->strand == it->strand) {
                      // Always calculate query distance the same way, as query always moves forward
                      int64_t query_dist = it2->queryStartPos - it->queryEndPos;

                      // Reference distance calculation depends on strand
                      int64_t ref_dist;
                      if (it->strand == strnd::FWD) {
                          ref_dist = it2->refStartPos - it->refEndPos;
                      } else {
                          // For reverse complement, we need to invert the order
                          ref_dist = it->refStartPos - it2->refEndPos;
                      }

                      // Check if the distance is within acceptable range
                      if (query_dist >= 0 && ref_dist >= -param.segLength/5 && ref_dist <= max_dist) {
                          double dist = std::sqrt(std::pow(query_dist, 2) + std::pow(ref_dist, 2));
                          if (dist < max_dist && best_score > dist && it2->chainPairScore > dist) {
                              best_it2 = it2;
                              best_score = dist;
                          }
                      }
                  }
              }
              if (best_it2 != readMappings.end()) {
                  best_it2->chainPairScore = best_score;
                  best_it2->chainPairId = it->splitMappingId;
              }
          }

          // Assign the merged mapping ids
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              it->splitMappingId = disjoint_sets.find(it->splitMappingId);
          }

          //Sort the mappings by post-merge split mapping id, then by query position, then by target position
          std::sort(
              readMappings.begin(),
              readMappings.end(),
              [](const MappingResult &a, const MappingResult &b) {
                  return std::tie(a.splitMappingId, a.queryStartPos, a.refSeqId, a.refStartPos, a.strand)
                      < std::tie(b.splitMappingId, b.queryStartPos, b.refSeqId, b.refStartPos, b.strand);
              });

          // Create maximallyMergedMappings with length restrictions
          MappingResultsVector_t maximallyMergedMappings;
          for(auto it = readMappings.begin(); it != readMappings.end();) {
              // Find the end of the current chain
              auto it_end = std::find_if(it, readMappings.end(), [&](const MappingResult &e){return e.splitMappingId != it->splitMappingId;} );
              
              // Process the chain into fragments that respect max_mapping_length
              auto fragment_start = it;
              while (fragment_start != it_end) {
                  MappingResult mergedMapping = *fragment_start;
                  auto fragment_end = fragment_start;
                  auto next = std::next(fragment_end);
                  
                  // Always try to merge enough mappings to satisfy minimum block length
                  offset_t current_length = fragment_start->queryEndPos - fragment_start->queryStartPos;
                  
                  while (next != it_end) {
                      offset_t potential_length = next->queryEndPos - fragment_start->queryStartPos;

                      if (potential_length > param.max_mapping_length) {
                          break;
                      }
                      fragment_end = next;
                      current_length = potential_length;
                      next = std::next(next);
                  }

                  // Create merged mapping for this fragment
                  mergedMapping.queryStartPos = fragment_start->queryStartPos;
                  mergedMapping.queryEndPos = fragment_end->queryEndPos;
                  
                  // Handle reference coordinates based on strand
                  if (mergedMapping.strand == strnd::FWD) {
                      // Forward strand - use first mapping's start and last mapping's end
                      mergedMapping.refStartPos = fragment_start->refStartPos;
                      mergedMapping.refEndPos = fragment_end->refEndPos;
                  } else {
                      // Reverse strand - use last mapping's start (highest coordinate) and first mapping's end (lowest coordinate)
                      mergedMapping.refStartPos = fragment_end->refStartPos;
                      mergedMapping.refEndPos = fragment_start->refEndPos;
                  }
                  mergedMapping.blockLength = std::max(
                      mergedMapping.refEndPos - mergedMapping.refStartPos,
                      mergedMapping.queryEndPos - mergedMapping.queryStartPos
                  );
                  
                  // Calculate averages for the fragment
                  double totalNucIdentity = 0.0;
                  double totalKmerComplexity = 0.0;
                  int totalConservedSketches = 0;
                  int totalSketchSize = 0;
                  int fragment_size = 0;
                  
                  for (auto subIt = fragment_start; subIt != std::next(fragment_end); ++subIt) {
                      totalNucIdentity += subIt->nucIdentity;
                      totalKmerComplexity += subIt->kmerComplexity;
                      totalConservedSketches += subIt->conservedSketches;
                      totalSketchSize += subIt->sketchSize;
                      fragment_size++;
                  }
    
                  mergedMapping.n_merged = fragment_size;
                  mergedMapping.nucIdentity = totalNucIdentity / fragment_size;
                  mergedMapping.kmerComplexity = totalKmerComplexity / fragment_size;
                  mergedMapping.conservedSketches = totalConservedSketches;
                  mergedMapping.sketchSize = totalSketchSize;
                  mergedMapping.blockNucIdentity = mergedMapping.nucIdentity;
                  mergedMapping.approxMatches = std::round(mergedMapping.nucIdentity * mergedMapping.blockLength / 100.0);
                  mergedMapping.discard = 0;
                  mergedMapping.overlapped = false;
                  mergedMapping.chainPairScore = std::numeric_limits<double>::max();
                  mergedMapping.chainPairId = std::numeric_limits<int64_t>::min();
                  
                  maximallyMergedMappings.push_back(mergedMapping);
                  
                  // Move to next fragment
                  fragment_start = std::next(fragment_end);
              }
              
              // Process the chain into chunks defined by max_mapping_length
              processChainWithSplits(it, it_end);

              // Move the iterator to the end of the processed chain
              it = it_end;
          }

          // After processing all chains, remove discarded mappings
          readMappings.erase(
              std::remove_if(readMappings.begin(), readMappings.end(), 
                             [](const MappingResult& e) { return e.discard == 1; }),
              readMappings.end()
          );

          return maximallyMergedMappings;
      }

      /**
       * @brief Process a chain of mappings, potentially splitting it into smaller fragments
       * @param begin Iterator to the start of the chain
       * @param end Iterator to the end of the chain
       */
      void processChainWithSplits(std::vector<MappingResult>::iterator begin, std::vector<MappingResult>::iterator end) {
          if (begin == end) return;

          std::vector<bool> is_cuttable(std::distance(begin, end), true);
          
          // Mark positions that are not cuttable (near discontinuities)
          for (auto it = std::next(begin); it != end; ++it) {
              auto prev = std::prev(it);
              if (it->queryStartPos - prev->queryEndPos > param.segLength / 5 ||
                  it->refStartPos - prev->refEndPos > param.segLength / 5) {
                  is_cuttable[std::distance(begin, prev)] = false;
                  is_cuttable[std::distance(begin, it)] = false;
              }
          }

          adjustConsecutiveMappings(begin, end, param.segLength);

          auto fragment_start = begin;
          offset_t accumulate_length = 0;

          for (auto it = begin; it != end; ++it) {
              accumulate_length += it->queryEndPos - it->queryStartPos;
              
              if (accumulate_length >= param.max_mapping_length && is_cuttable[std::distance(begin, it)]) {
                  // Process the fragment up to this point
                  processMappingFragment(fragment_start, std::next(it));
                  
                  // Start a new fragment
                  fragment_start = std::next(it);
                  accumulate_length = 0;
              }
          }

          // Process any remaining fragment
          if (fragment_start != end) {
              processMappingFragment(fragment_start, end);
          }

          // Compute and assign chain statistics
          computeChainStatistics(begin, end);
      }

      /**
       * @brief Compute and assign chain statistics to all mappings in the chain
       * @param begin Iterator to the start of the chain
       * @param end Iterator to the end of the chain
       */
      /**
       * @brief Filter mappings within a subset before aggregation
       * @param mappings Mappings to filter
       * @param param Algorithm parameters
       */
      std::pair<MappingResultsVector_t, MappingResultsVector_t> filterSubsetMappings(MappingResultsVector_t& mappings, progress_meter::ProgressMeter& progress) {
          if (mappings.empty()) return {MappingResultsVector_t(), MappingResultsVector_t()};
          
          // Only merge once and keep both versions
          auto maximallyMergedMappings = mergeMappingsInRange(mappings, param.chain_gap, progress);
          
          // Process both merged and non-merged mappings
          if (param.mergeMappings && param.split) {
              filterMaximallyMerged(maximallyMergedMappings, param, progress);
          } else {
              filterNonMergedMappings(mappings, param, progress);
          }

          // Build dense chain ID mapping
          std::unordered_map<offset_t, offset_t> id_map;
          offset_t next_id = 0;
          
          // First pass - build the mapping from both sets
          for (const auto& mapping : mappings) {
              if (id_map.count(mapping.splitMappingId) == 0) {
                  id_map[mapping.splitMappingId] = next_id++;
              }
          }
          for (const auto& mapping : maximallyMergedMappings) {
              if (id_map.count(mapping.splitMappingId) == 0) {
                  id_map[mapping.splitMappingId] = next_id++;
              }
          }

          // Get atomic offset for this batch of chain IDs
          offset_t base_id = maxChainIdSeen.fetch_add(id_map.size(), std::memory_order_relaxed);
          
          // Apply compacted IDs with offset
          for (auto& mapping : mappings) {
              mapping.splitMappingId = id_map[mapping.splitMappingId] + base_id;
          }
          for (auto& mapping : maximallyMergedMappings) {
              mapping.splitMappingId = id_map[mapping.splitMappingId] + base_id;
          }

          return {std::move(mappings), std::move(maximallyMergedMappings)};
      }

      /**
       * @brief Update chain IDs to prevent conflicts between subsets
       * @param mappings Mappings whose chain IDs need updating
       * @param maxId Current maximum chain ID seen
       */

      void computeChainStatistics(std::vector<MappingResult>::iterator begin, std::vector<MappingResult>::iterator end) {
          offset_t chain_start_query = std::numeric_limits<offset_t>::max();
          offset_t chain_end_query = std::numeric_limits<offset_t>::min();
          offset_t chain_start_ref = std::numeric_limits<offset_t>::max();
          offset_t chain_end_ref = std::numeric_limits<offset_t>::min();
          double accumulate_nuc_identity = 0.0;
          int n_in_full_chain = std::distance(begin, end);

          for (auto it = begin; it != end; ++it) {
              chain_start_query = std::min(chain_start_query, it->queryStartPos);
              chain_end_query = std::max(chain_end_query, it->queryEndPos);
              chain_start_ref = std::min(chain_start_ref, it->refStartPos);
              chain_end_ref = std::max(chain_end_ref, it->refEndPos);
              accumulate_nuc_identity += it->nucIdentity;
          }

          auto chain_nuc_identity = accumulate_nuc_identity / n_in_full_chain;
          auto block_length = std::max(chain_end_query - chain_start_query, chain_end_ref - chain_start_ref);

          for (auto it = begin; it != end; ++it) {
              it->n_merged = n_in_full_chain;
              it->blockLength = block_length;
              it->blockNucIdentity = chain_nuc_identity;
          }
      }

     /**
       * @brief                       This routine is to make sure that all mapping boundaries
       *                              on query and reference are not outside total
       *                              length of sequeunces involved
       * @param[in]     input         input read details
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void mappingBoundarySanityCheck(InputSeqProgContainer* input, VecIn &readMappings)
        {
          for(auto &e : readMappings)
          {
            //reference start pos
            {
              if(e.refStartPos < 0)
                e.refStartPos = 0;
              if(e.refStartPos >= this->idManager->getSequenceLength(e.refSeqId))
                e.refStartPos = this->idManager->getSequenceLength(e.refSeqId) - 1;
            }

            //reference end pos
            {
              if(e.refEndPos < e.refStartPos)
                e.refEndPos = e.refStartPos;
              if(e.refEndPos >= this->idManager->getSequenceLength(e.refSeqId))
                e.refEndPos = this->idManager->getSequenceLength(e.refSeqId) - 1;
            }

            //query start pos
            {
              if(e.queryStartPos < 0)
                e.queryStartPos = 0;
              if(e.queryStartPos >= input->len)
                e.queryStartPos = input->len;
            }

            //query end pos
            {
              if(e.queryEndPos < e.queryStartPos)
                e.queryEndPos = e.queryStartPos;
              if(e.queryEndPos >= input->len)
                e.queryEndPos = input->len;
            }
          }
        }

      /**
       * @brief                         Report the final read mappings to output stream
       * @param[in]   readMappings      mapping results for single or multiple reads
       * @param[in]   queryName         input required if reporting one read at a time
       * @param[in]   outstrm           file output stream object
       */
      void reportReadMappings(MappingResultsVector_t &readMappings, const std::string &queryName,
          std::ostream &outstrm)
      {
        // Sort mappings by chain ID and query position
        std::sort(readMappings.begin(), readMappings.end(),
            [](const MappingResult &a, const MappingResult &b) {
                return std::tie(a.splitMappingId, a.queryStartPos)
                    < std::tie(b.splitMappingId, b.queryStartPos);
            });

        // Assign chain positions within each chain
        int current_chain = -1;
        int chain_pos = 0;
        int chain_length = 0;
        
        // First pass - count chain lengths
        for (size_t i = 0; i < readMappings.size(); ++i) {
            if (readMappings[i].splitMappingId != current_chain) {
                current_chain = readMappings[i].splitMappingId;
                chain_length = 1;
                // Count forward to find chain length
                for (size_t j = i + 1; j < readMappings.size(); ++j) {
                    if (readMappings[j].splitMappingId == current_chain) {
                        chain_length++;
                    } else {
                        break;
                    }
                }
                // Assign length to all mappings in this chain
                readMappings[i].chain_length = chain_length;
                chain_pos = 1;
            } else {
                readMappings[i].chain_length = chain_length;
                chain_pos++;
            }
            readMappings[i].chain_pos = chain_pos;
        }

        //Print the results
        for(auto &e : readMappings)
        {
          float fakeMapQ = e.nucIdentity == 1 ? 255 : std::round(-10.0 * std::log10(1-(e.nucIdentity)));
          std::string sep = param.legacy_output ? " " : "\t";

          outstrm  << (param.filterMode == filter::ONETOONE ? idManager->getSequenceName(e.querySeqId) : queryName)
                   << sep << e.queryLen
                   << sep << e.queryStartPos
                   << sep << e.queryEndPos - (param.legacy_output ? 1 : 0)
                   << sep << (e.strand == strnd::FWD ? "+" : "-")
                   << sep << idManager->getSequenceName(e.refSeqId)
                   << sep << this->idManager->getSequenceLength(e.refSeqId)
                   << sep << e.refStartPos
                   << sep << e.refEndPos - (param.legacy_output ? 1 : 0);

          if (!param.legacy_output) 
          {
            outstrm  << sep << e.conservedSketches
                     << sep << e.blockLength
                     << sep << fakeMapQ
                     << sep << "id:f:" << e.nucIdentity
                     << sep << "kc:f:" << e.kmerComplexity;
            if (!param.mergeMappings) 
            {
              outstrm << sep << "jc:f:" << float(e.conservedSketches) / e.sketchSize;
            } else {
              outstrm << sep << "chain:i:" << e.splitMappingId << "." << e.chain_pos << "." << e.chain_length;
            }
          } else
          {
            outstrm << sep << e.nucIdentity * 100.0;
          }

#ifdef DEBUG
          outstrm << std::endl;
#else
          outstrm << "\n";
#endif

          //User defined processing of the results
          if(processMappingResults != nullptr)
            processMappingResults(e);
        }
      }

    private:
      void processCombinedMappingsThread(aggregate_atomic_queue_t& aggregate_queue, writer_atomic_queue_t& writer_queue, std::atomic<bool>& processing_done, progress_meter::ProgressMeter& progress) {
          int wait_count = 0;
          while (true) {
              std::pair<seqno_t, MappingResultsVector_t*>* task = nullptr;
              if (aggregate_queue.try_pop(task)) {
                  wait_count = 0;
                  auto querySeqId = task->first;
                  auto& mappings = *(task->second);
                  
                  std::string queryName = idManager->getSequenceName(querySeqId);
                  // Final filtering pass on pre-filtered mappings
                  if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                      MappingResultsVector_t filteredMappings;
                      filterByGroup(mappings, filteredMappings, param.numMappingsForSegment - 1, 
                                  param.filterMode == filter::ONETOONE, *idManager, progress);
                      mappings = std::move(filteredMappings);
                  }

                  std::stringstream ss;
                  reportReadMappings(mappings, queryName, ss);
                  
                  auto* output = new std::string(ss.str());
                  while (!writer_queue.try_push(output)) {
                      std::this_thread::sleep_for(std::chrono::milliseconds(10));
                  }
                  delete task;
              } else {
                  if (processing_done.load() && aggregate_queue.was_empty()) {
                      ++wait_count;
                  }
                  if (wait_count > 10) {
                      break;
                  }
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }
      }

      void outputThread(std::ofstream& outstrm, writer_atomic_queue_t& writer_queue, std::atomic<bool>& processing_done,
                        std::atomic<bool>& workers_done, std::atomic<bool>& output_done) {
          int wait_count = 0;
          while (!output_done.load()) {
              std::string* result = nullptr;
              if (writer_queue.try_pop(result)) {
                  wait_count = 0;
                  outstrm << *result;
                  delete result;
              } else {
                  if (processing_done.load() && workers_done.load() && writer_queue.was_empty()) {
                      ++wait_count;
                  }
                  if (wait_count > 10) {
                      output_done.store(true);
                  }
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }
      }

    public:

      /**
       * @brief     An optional utility function to save the
       *            reported results by the L2 stage into a vector
       */
      static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
      {
        v.push_back(reportedL2Result);
      }

  };

}

#endif
