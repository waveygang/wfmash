/**
 * @file    computeMap.hpp
 * @brief   Main Map class that implements the sequence mapping logic with performance optimizations
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
#include <map>
namespace fs = std::filesystem;
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <chrono>
#include <sstream>
#include <semaphore>
#include <memory_resource>
#include "taskflow/taskflow.hpp"

// Include the reentrant FASTA/BGZF index implementation
#define REENTRANT_FAIDX_IMPLEMENTATION
#include "common/faigz.h"

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/map_stats.hpp"

// Include the new separated files
#include "map/include/fragmentManager.hpp"
#include "map/include/mappingCore.hpp"
#include "map/include/mappingFilter.hpp"
#include "map/include/mappingOutput.hpp"

//External includes
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "gsl/gsl_randist.h"

namespace skch
{
  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages with performance optimizations
   */
  class Map
  {
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

      // Track maximum chain ID seen across all subsets
      std::atomic<offset_t> maxChainIdSeen{0};

      //Cache for commonly used values
      offset_t cached_segment_length;
      int cached_minimum_hits;

      // Type aliases for core mapping
      using CoreMapping = MappingCore<Sketch, SequenceIdManager>;
      using FilterUtils = MappingFilterUtils;
      using OutputHandler = MappingOutput;

      // Performance optimization: Configurable thresholds
      static constexpr size_t SMALL_QUERY_THRESHOLD = 100;  // fragments
      static constexpr size_t FRAGMENT_BATCH_SIZE = 16;     // fragments per batch
      static constexpr std::ptrdiff_t MAX_CONCURRENT_QUERIES = 8; // max queries in flight
      static constexpr size_t MEMORY_POOL_SIZE = 64 * 1024 * 1024; // 64MB per thread

      // Structure to hold fragment processing context with pre-allocated memory
      struct FragmentContext {
          std::vector<IntervalPoint> intervalPoints;
          std::vector<L1_candidateLocus_t> l1Mappings;
          MappingResultsVector_t l2Mappings;
          QueryMetaData<MinVec_Type> Q;
          std::vector<MappingResult> results;
          
          // Pre-allocate capacity for better performance
          FragmentContext() {
              intervalPoints.reserve(1000);
              l1Mappings.reserve(100);
              l2Mappings.reserve(100);
              results.reserve(500);
          }
          
          void clear() {
              intervalPoints.clear();
              l1Mappings.clear();
              l2Mappings.clear();
              results.clear();
          }
      };

      // Thread-local memory pool for fragment contexts
      static thread_local std::unique_ptr<FragmentContext> tls_fragmentContext;
      
      // Semaphore for stream processing control
      std::counting_semaphore<MAX_CONCURRENT_QUERIES> queryProcessingSemaphore{MAX_CONCURRENT_QUERIES};

      // Batch of fragments for processing
      struct FragmentBatch {
          std::vector<FragmentData> fragments;
          MappingResultsVector_t results;
          
          void reserve(size_t size) {
              fragments.reserve(size);
              results.reserve(size * 10); // Estimate 10 mappings per fragment
          }
      };

      void processFragment(const FragmentData& fragment, FragmentContext& ctx) {
          ctx.intervalPoints.clear();
          ctx.l1Mappings.clear();
          ctx.l2Mappings.clear();
          ctx.results.clear();

          ctx.Q.seq = const_cast<char*>(fragment.seq);
          ctx.Q.len = fragment.len;
          ctx.Q.fullLen = fragment.fullLen;
          ctx.Q.seqId = fragment.seqId;
          ctx.Q.seqName = fragment.seqName;
          ctx.Q.refGroup = fragment.refGroup;

          mapSingleQueryFrag(ctx.Q, ctx.intervalPoints, ctx.l1Mappings, ctx.l2Mappings);

          std::for_each(ctx.l2Mappings.begin(), ctx.l2Mappings.end(), [&](MappingResult &e){
              // For compact struct: adjust query start position to be relative to full sequence
              e.queryStartPos += fragment.fragmentIndex * param.windowLength;
              // blockLength remains as fragment length (already set during mapping creation)
          });

          if (!ctx.l2Mappings.empty()) {
              ctx.results.insert(ctx.results.end(), 
                                std::make_move_iterator(ctx.l2Mappings.begin()), 
                                std::make_move_iterator(ctx.l2Mappings.end()));
          }
          
          if (fragment.output) {
              fragment.output->progress.increment(fragment.len);
          }
      }

      // Process a batch of fragments
      void processFragmentBatch(const FragmentBatch& batch, FragmentContext& ctx,
                               std::shared_ptr<progress_meter::ProgressMeter> progress) {
          for (const auto& fragment : batch.fragments) {
              ctx.clear();
              
              ctx.Q.seq = const_cast<char*>(fragment.seq);
              ctx.Q.len = fragment.len;
              ctx.Q.fullLen = fragment.fullLen;
              ctx.Q.seqId = fragment.seqId;
              ctx.Q.seqName = fragment.seqName;
              ctx.Q.refGroup = fragment.refGroup;

              mapSingleQueryFrag(ctx.Q, ctx.intervalPoints, ctx.l1Mappings, ctx.l2Mappings);

              std::for_each(ctx.l2Mappings.begin(), ctx.l2Mappings.end(), [&](MappingResult &e){
                  e.queryStartPos += fragment.fragmentIndex * param.windowLength;
              });

              if (!ctx.l2Mappings.empty()) {
                  ctx.results.reserve(ctx.results.size() + ctx.l2Mappings.size());
                  ctx.results.insert(ctx.results.end(), 
                                    std::make_move_iterator(ctx.l2Mappings.begin()), 
                                    std::make_move_iterator(ctx.l2Mappings.end()));
              }
              
              progress->increment(fragment.len);
          }
      }

      // Structure to hold query processing results
      struct QueryResults {
          std::string queryName;
          seqno_t seqId;
          std::string sequence;
          MappingResultsVector_t mappings;
          int refGroup;
          offset_t queryLen;
          size_t fragmentCount;
          
          QueryResults() {
              mappings.reserve(1000); // Pre-allocate for typical case
          }
      };
      
    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      Map(skch::Parameters p, PostProcessResultsFn_t f = nullptr) :
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
        cached_segment_length(p.windowLength),
        cached_minimum_hits(std::max(p.minimum_hits, Stat::estimateMinimumHitsRelaxed(p.sketchSize, p.kmerSize, p.percentageIdentity, skch::fixed::confidence_interval)))
      {
          // Initialize sequence names
          if (!param.query_prefix.empty()) {
              this->querySequenceNames.clear();
              for (const auto& name : idManager->getQuerySequenceNames()) {
                  bool prefix_match = false;
                  for (const auto& prefix : param.query_prefix) {
                      if (name.compare(0, prefix.size(), prefix) == 0) {
                          prefix_match = true;
                          break;
                      }
                  }
                  if (prefix_match) {
                      this->querySequenceNames.push_back(name);
                  }
              }
          } else {
              this->querySequenceNames = idManager->getQuerySequenceNames();
          }
          
          if (!param.target_prefix.empty()) {
              this->targetSequenceNames.clear();
              for (const auto& name : idManager->getTargetSequenceNames()) {
                  if (name.compare(0, param.target_prefix.size(), param.target_prefix) == 0) {
                      this->targetSequenceNames.push_back(name);
                  }
              }
          } else {
              this->targetSequenceNames = idManager->getTargetSequenceNames();
          }

          // Calculate and log statistics
          uint64_t total_target_length = 0;
          size_t target_seq_count = targetSequenceNames.size();
          for (const auto& seqName : targetSequenceNames) {
              seqno_t seqId = idManager->getSequenceId(seqName);
              total_target_length += idManager->getSequenceLength(seqId);
          }

          uint64_t total_query_length = 0;
          for (const auto& seqName : querySequenceNames) {
              total_query_length += idManager->getSequenceLength(idManager->getSequenceId(seqName));
          }

          std::unordered_set<int> query_groups, target_groups;
          for (const auto& seqName : querySequenceNames) {
              query_groups.insert(idManager->getRefGroup(idManager->getSequenceId(seqName)));
          }
          for (const auto& seqName : targetSequenceNames) {
              target_groups.insert(idManager->getRefGroup(idManager->getSequenceId(seqName)));
          }

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

      ~Map() = default;

    private:

      void setProbs()
      {
        float deltaANI = param.ANIDiff;
        float min_p = 1 - param.ANIDiffConf;
        int ss = std::min<double>(param.sketchSize, skch::fixed::ss_table_max);

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
        
        const auto distDiff = [this, &sketchProbs, deltaANI, min_p, ss] (int cmax, int ci) {
          double prAboveCutoff = 0;
          for (double ymax = 0; ymax <= cmax; ymax++) {
            double pymax = sketchProbs[cmax][ymax];
            double yi_cutoff = deltaANI == 0 ? ymax : (std::floor(skch::Stat::md2j(
                skch::Stat::j2md(ymax / ss, param.kmerSize) + deltaANI, 
                param.kmerSize
            ) * ss));
            double pi_acc = (yi_cutoff - 1) >= 0 ? gsl_cdf_hypergeometric_P(yi_cutoff-1, ss, ss-ci, ci) : 0;
            pi_acc = 1-pi_acc;
            prAboveCutoff += pymax * pi_acc;
            if (prAboveCutoff > min_p)
            {
              return true;
            }
          }
          return prAboveCutoff > min_p; 
        };

        std::vector<int> ss_range(ss+1);
        std::iota (ss_range.begin(), ss_range.end(), 0);

        for (auto cmax = 1; cmax <= ss; cmax++) 
        {
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

          if (sketchCutoffs[cmax] == 0) {
            sketchCutoffs[cmax] = 1;
          }
        }
      }

      std::vector<std::vector<std::string>> createTargetSubsets(const std::vector<std::string>& targetSequenceNames) {
        std::vector<std::vector<std::string>> target_subsets;
        uint64_t current_subset_size = 0;
        std::vector<std::string> current_subset;

        int64_t batch_size = param.index_by_size;
        if (batch_size <= 0) {
            batch_size = 5000000;
            if (!param.indexFilename.empty() && !param.create_index_only) {
                std::cerr << "[wfmash::mashmap] Warning: Invalid batch size in index, using default: " 
                          << batch_size << std::endl;
            }
        }

        for (const auto& seqName : targetSequenceNames) {
            seqno_t seqId = idManager->getSequenceId(seqName);
            offset_t seqLen = idManager->getSequenceLength(seqId);
            current_subset.push_back(seqName);
            current_subset_size += seqLen;

            if (current_subset_size >= batch_size || &seqName == &targetSequenceNames.back()) {
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

      void mapQuery() {
          tf::Executor executor(param.threads);
          
          std::unordered_map<seqno_t, MappingResultsVector_t> combinedMappings;
          std::mutex combinedMappings_mutex;

          if (!param.indexFilename.empty() && !param.create_index_only) {
              std::ifstream indexStream(param.indexFilename.string(), std::ios::binary);
              if (!indexStream) {
                  std::cerr << "Error: Unable to open index file for reading: " << param.indexFilename << std::endl;
                  exit(1);
              }
              
              uint64_t magic_number = 0;
              indexStream.read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
              if (magic_number != 0xDEADBEEFCAFEBABE) {
                  std::cerr << "Error: Invalid index file format (wrong magic number)" << std::endl;
                  exit(1);
              }
              
              size_t batch_idx, total_batches;
              indexStream.read(reinterpret_cast<char*>(&batch_idx), sizeof(batch_idx));
              indexStream.read(reinterpret_cast<char*>(&total_batches), sizeof(total_batches));
              
              int64_t batch_size = 0;
              indexStream.read(reinterpret_cast<char*>(&batch_size), sizeof(batch_size));
              if (batch_size > 0) {
                  param.index_by_size = batch_size;
                  std::cerr << "[wfmash::mashmap] Using batch size " << batch_size 
                            << " from index file (" << total_batches << " subsets)" << std::endl;
              }
              indexStream.close();
          }

          auto target_subsets = createTargetSubsets(targetSequenceNames);

          uint64_t total_target_subset_size = 0;
          for (const auto& subset : target_subsets) {
              for (const auto& seqName : subset) {
                  seqno_t seqId = idManager->getSequenceId(seqName);
                  total_target_subset_size += idManager->getSequenceLength(seqId);
              }
          }
          double avg_subset_size = target_subsets.size() ? 
              (double)total_target_subset_size / target_subsets.size() : 0;

          std::cerr << "[wfmash::mashmap] Processing " << target_subsets.size() 
                    << " target subsets (≈" << std::fixed << std::setprecision(0) << avg_subset_size 
                    << "bp/subset)" << std::endl;

          bool exit_after_indices = param.create_index_only;

          // Process each subset serially
          for (size_t subset_idx = 0; subset_idx < target_subsets.size(); ++subset_idx) {
              const auto& target_subset = target_subsets[subset_idx];
              if(target_subset.empty()) continue;
              
              std::cerr << "[wfmash::mashmap] Processing subset " << (subset_idx + 1) 
                        << "/" << target_subsets.size() << " (mapping)" << std::endl;
              
              std::string indexFilename = param.indexFilename.string();
              
              if (param.create_index_only) {
                  std::cerr << "[wfmash::mashmap] Processing subset " << (subset_idx + 1) 
                            << "/" << target_subsets.size() << " (indexing): " << indexFilename << std::endl;
    
                  refSketch = new skch::Sketch(param, *idManager, target_subset);
                  bool append = (subset_idx > 0);
                  refSketch->writeIndex(target_subset, indexFilename, append, subset_idx, target_subsets.size());
                  delete refSketch;
                  refSketch = nullptr;
                  std::cerr << "[wfmash::mashmap] index construction completed" << std::endl;
                  continue;
              }
              
              // Build index for this subset
              buildIndexForSubset(target_subset, subset_idx, target_subsets.size());
              
              uint64_t subset_query_length = 0;
              for (const auto& queryName : querySequenceNames) {
                  subset_query_length += idManager->getSequenceLength(idManager->getSequenceId(queryName));
              }

              auto progress = std::make_shared<progress_meter::ProgressMeter>(
                  subset_query_length,
                  "[wfmash::mashmap] mapping ",
                  param.use_progress_bar
              );

              auto subsetMappings = std::make_shared<std::unordered_map<seqno_t, MappingResultsVector_t>>();
              auto subsetMappings_mutex = std::make_shared<std::mutex>();

              auto outstream = std::make_shared<std::ofstream>();
              
              if (param.filterMode != filter::ONETOONE) {
                  bool append = subset_idx > 0;
                  outstream->open(param.outFileName, append ? std::ios::app : std::ios::out);
                  if (!outstream->is_open()) {
                      std::cerr << "Error: Could not open output file for writing: " << param.outFileName << std::endl;
                      exit(1);
                  }
              }

              // Process all queries for this subset with optimized strategy
              processQueriesForSubsetOptimized(executor, progress, subsetMappings, subsetMappings_mutex, outstream);
              
              // Wait for mapping to complete and notify user
              progress->finish();
              std::cerr << "[wfmash::mashmap] Mapping phase complete, starting post-processing";
              if (param.scaffold_gap > 0) {
                  std::cerr << " (including scaffold filtering)";
              }
              std::cerr << "..." << std::endl;
              
              // Merge results if needed
              if (param.filterMode == filter::ONETOONE) {
                  std::lock_guard<std::mutex> lock(combinedMappings_mutex);
                  for (auto& [querySeqId, mappings] : *subsetMappings) {
                      combinedMappings[querySeqId].insert(
                          combinedMappings[querySeqId].end(),
                          std::make_move_iterator(mappings.begin()),
                          std::make_move_iterator(mappings.end())
                      );
                  }
              }
              
              // Cleanup
              if (param.filterMode != filter::ONETOONE && outstream->is_open()) {
                  outstream->flush();
                  outstream->close();
              }
              
              delete refSketch;
              refSketch = nullptr;
              
              std::cerr << "[wfmash::mashmap] Subset " << (subset_idx + 1) 
                        << "/" << target_subsets.size() << " processing complete" << std::endl;
          }

          if (exit_after_indices) {
              std::cerr << "[wfmash::mashmap] All indices created successfully. Exiting." << std::endl;
              exit(0);
          }

          // Final results processing for ONETOONE mode
          if (param.filterMode == filter::ONETOONE && !exit_after_indices) {
              processFinalOneToOneResults(executor, combinedMappings);
          }
          
          // Explicitly clear combinedMappings to free memory before alignment phase
          combinedMappings.clear();
          std::unordered_map<seqno_t, MappingResultsVector_t>().swap(combinedMappings);
      }

      void buildIndexForSubset(const std::vector<std::string>& target_subset, size_t subset_idx, size_t total_subsets) {
          if (!param.indexFilename.empty()) {
              std::string indexFilename = param.indexFilename.string();
              static std::ifstream indexStream;
              static bool index_opened = false;
              static size_t file_subset_count = 0;
              
              if (!index_opened) {
                  indexStream.open(indexFilename, std::ios::binary);
                  if (!indexStream) {
                      std::cerr << "Error: Unable to open index file for reading: " << indexFilename << std::endl;
                      exit(1);
                  }
                  
                  uint64_t magic_number = 0;
                  indexStream.read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
                  if (magic_number != 0xDEADBEEFCAFEBABE) {
                      std::cerr << "Error: Invalid index file format (wrong magic number)" << std::endl;
                      exit(1);
                  }
                  
                  size_t batch_idx, total_batches;
                  indexStream.read(reinterpret_cast<char*>(&batch_idx), sizeof(batch_idx));
                  indexStream.read(reinterpret_cast<char*>(&total_batches), sizeof(total_batches));
                  
                  file_subset_count = total_batches;
                  std::cerr << "[wfmash::mashmap] Index file contains " << file_subset_count 
                            << " subsets" << std::endl;
                  
                  int64_t batch_size = 0;
                  indexStream.read(reinterpret_cast<char*>(&batch_size), sizeof(batch_size));
                  if (batch_size > 0) {
                      param.index_by_size = batch_size;
                      std::cerr << "[wfmash::mashmap] Using batch size " << batch_size 
                                << " from index" << std::endl;
                  }
                  
                  indexStream.seekg(0, std::ios::beg);
                  index_opened = true;
              }
              
              refSketch = new skch::Sketch(param, *idManager, target_subset, &indexStream);
          } else {
              auto sketch_index_progress = std::make_shared<progress_meter::ProgressMeter>(
                  100,
                  "[wfmash::mashmap] indexing",
                  param.use_progress_bar);
          
              if (!sketch_index_progress->is_finished.load()) {
                  sketch_index_progress->reset_timer();
              }
              
              refSketch = new skch::Sketch(param, *idManager, target_subset, nullptr, sketch_index_progress);
          }
      }

      // Optimized version with hybrid parallel strategy
      void processQueriesForSubsetOptimized(tf::Executor& executor,
                                          std::shared_ptr<progress_meter::ProgressMeter> progress,
                                          std::shared_ptr<std::unordered_map<seqno_t, MappingResultsVector_t>> subsetMappings,
                                          std::shared_ptr<std::mutex> subsetMappings_mutex,
                                          std::shared_ptr<std::ofstream> outstream) {
          
          // Prepare all query tasks first
          std::vector<QueryResults> allQueryResults;
          allQueryResults.reserve(querySequenceNames.size());
          
          // Load all query sequences first (sequential I/O)
          const auto& fileName = param.querySequences[0];
          faidx_meta_t* query_meta = faidx_meta_load(fileName.c_str(), FAI_FASTA, FAI_CREATE);
          if (!query_meta) {
              std::cerr << "Error: Failed to load query FASTA index: " << fileName << std::endl;
              exit(1);
          }
          
          // Create a single reader for sequential access
          faidx_reader_t* reader = faidx_reader_create(query_meta);
          if (!reader) {
              std::cerr << "Error: Failed to create FASTA reader" << std::endl;
              faidx_meta_destroy(query_meta);
              exit(1);
          }
          
          // Read all sequences sequentially
          for (const auto& queryName : querySequenceNames) {
              hts_pos_t seq_len;
              hts_pos_t seq_total_len = faidx_meta_seq_len(query_meta, queryName.c_str());
              if (seq_total_len <= 0) {
                  std::cerr << "Warning: Sequence " << queryName << " not found or empty, skipping" << std::endl;
                  continue;
              }
              
              char* seq_data = faidx_reader_fetch_seq(reader, queryName.c_str(), 0, seq_total_len-1, &seq_len);
              if (!seq_data) {
                  std::cerr << "Warning: Failed to fetch sequence " << queryName << ", skipping" << std::endl;
                  continue;
              }
              
              QueryResults qr;
              qr.queryName = queryName;
              qr.seqId = idManager->getSequenceId(queryName);
              qr.sequence = std::string(seq_data, seq_len);
              qr.refGroup = idManager->getRefGroup(qr.seqId);
              qr.queryLen = seq_len;
              qr.fragmentCount = (qr.queryLen / param.windowLength) + 
                                ((qr.queryLen % param.windowLength != 0) ? 1 : 0);
              allQueryResults.push_back(std::move(qr));
              
              free(seq_data);
          }
          
          faidx_reader_destroy(reader);
          faidx_meta_destroy(query_meta);
          
          // Sort queries by size (largest first) for better load balancing
          std::sort(allQueryResults.begin(), allQueryResults.end(), 
                   [](const QueryResults& a, const QueryResults& b) {
                       return a.fragmentCount > b.fragmentCount;
                   });
          
          // Process queries with stream control
          tf::Taskflow taskflow;
          
          for (auto& queryData : allQueryResults) {
              auto task = taskflow.emplace([this, &queryData, &executor, progress, subsetMappings, subsetMappings_mutex, outstream]() {
                  // Stream control: acquire semaphore
		  this->queryProcessingSemaphore.acquire();
                  
                  // Choose processing strategy based on query size
                  if (queryData.fragmentCount < SMALL_QUERY_THRESHOLD) {
                      // Small query: process serially
                      processQuerySerial(queryData, progress, subsetMappings, subsetMappings_mutex, outstream);
                  } else {
                      // Large query: process fragments in parallel
                      processQueryParallel(queryData, executor, progress, subsetMappings, subsetMappings_mutex, outstream);
                  }
                  
                  // Release semaphore
                  this->queryProcessingSemaphore.release();
              });
          }
          
          // Execute all query tasks
          executor.run(taskflow).wait();
      }

      // Serial processing for small queries
      void processQuerySerial(QueryResults& queryData,
                            std::shared_ptr<progress_meter::ProgressMeter> progress,
                            std::shared_ptr<std::unordered_map<seqno_t, MappingResultsVector_t>> subsetMappings,
                            std::shared_ptr<std::mutex> subsetMappings_mutex,
                            std::shared_ptr<std::ofstream> outstream) {
          
          // Get thread-local fragment context
          if (!tls_fragmentContext) {
              tls_fragmentContext = std::make_unique<FragmentContext>();
          }
          auto& fragmentCtx = *tls_fragmentContext;
          
          MappingResultsVector_t allMappings;
          allMappings.reserve(queryData.fragmentCount * 10); // Estimate
          
          int noOverlapFragmentCount = queryData.queryLen / param.windowLength;
          
          // Process regular fragments
          for(int i = 0; i < noOverlapFragmentCount; i++) {
              FragmentData fragment(
                  queryData.sequence.data() + i * param.windowLength,
                  static_cast<int>(param.windowLength),
                  static_cast<int>(queryData.queryLen),
                  queryData.seqId,
                  queryData.queryName,
                  queryData.refGroup,
                  i,
                  nullptr
              );
              
              processFragment(fragment, fragmentCtx);
              
              if (!fragmentCtx.results.empty()) {
                  allMappings.insert(
                      allMappings.end(),
                      std::make_move_iterator(fragmentCtx.results.begin()),
                      std::make_move_iterator(fragmentCtx.results.end())
                  );
              }
              
              progress->increment(param.windowLength);
          }
          
          // Handle final fragment if needed
          if (noOverlapFragmentCount >= 1 && queryData.queryLen % param.windowLength != 0) {
              FragmentData fragment(
                  queryData.sequence.data() + queryData.queryLen - param.windowLength,
                  static_cast<int>(param.windowLength),
                  static_cast<int>(queryData.queryLen),
                  queryData.seqId,
                  queryData.queryName,
                  queryData.refGroup,
                  noOverlapFragmentCount,
                  nullptr
              );
              
              processFragment(fragment, fragmentCtx);
              
              if (!fragmentCtx.results.empty()) {
                  allMappings.insert(
                      allMappings.end(),
                      std::make_move_iterator(fragmentCtx.results.begin()),
                      std::make_move_iterator(fragmentCtx.results.end())
                  );
              }
              
              progress->increment(queryData.queryLen % param.windowLength);
          }
          
          // Post-process and output results
          finalizeQueryResults(queryData, allMappings, progress, subsetMappings, subsetMappings_mutex, outstream);
      }

      // Parallel processing for large queries
      void processQueryParallel(QueryResults& queryData,
                              tf::Executor& executor,
                              std::shared_ptr<progress_meter::ProgressMeter> progress,
                              std::shared_ptr<std::unordered_map<seqno_t, MappingResultsVector_t>> subsetMappings,
                              std::shared_ptr<std::mutex> subsetMappings_mutex,
                              std::shared_ptr<std::ofstream> outstream) {
          
          // Create fragment batches
          std::vector<FragmentBatch> batches;
          int noOverlapFragmentCount = queryData.queryLen / param.windowLength;
          size_t totalFragments = noOverlapFragmentCount;
          if (noOverlapFragmentCount >= 1 && queryData.queryLen % param.windowLength != 0) {
              totalFragments++;
          }
          
          // Calculate number of batches
          size_t numBatches = (totalFragments + FRAGMENT_BATCH_SIZE - 1) / FRAGMENT_BATCH_SIZE;
          batches.resize(numBatches);
          
          // Distribute fragments to batches
          size_t batchIdx = 0;
          for(int i = 0; i < noOverlapFragmentCount; i++) {
              FragmentData fragment(
                  queryData.sequence.data() + i * param.windowLength,
                  static_cast<int>(param.windowLength),
                  static_cast<int>(queryData.queryLen),
                  queryData.seqId,
                  queryData.queryName,
                  queryData.refGroup,
                  i,
                  nullptr
              );
              
              batches[batchIdx].fragments.push_back(fragment);
              if (batches[batchIdx].fragments.size() >= FRAGMENT_BATCH_SIZE) {
                  batchIdx++;
              }
          }
          
          // Handle final fragment if needed
          if (noOverlapFragmentCount >= 1 && queryData.queryLen % param.windowLength != 0) {
              FragmentData fragment(
                  queryData.sequence.data() + queryData.queryLen - param.windowLength,
                  static_cast<int>(param.windowLength),
                  static_cast<int>(queryData.queryLen),
                  queryData.seqId,
                  queryData.queryName,
                  queryData.refGroup,
                  noOverlapFragmentCount,
                  nullptr
              );
              
              if (batchIdx >= batches.size()) batchIdx = batches.size() - 1;
              batches[batchIdx].fragments.push_back(fragment);
          }
          
          // Process batches in parallel
          tf::Taskflow batchTaskflow;
          std::mutex results_mutex;
          MappingResultsVector_t allMappings;
          
          for (auto& batch : batches) {
              batch.reserve(batch.fragments.size());
              
              batchTaskflow.emplace([this, &batch, progress, &results_mutex, &allMappings]() {
                  // Get thread-local fragment context
                  if (!tls_fragmentContext) {
                      tls_fragmentContext = std::make_unique<FragmentContext>();
                  }
                  auto& ctx = *tls_fragmentContext;
                  
                  // Process all fragments in batch
                  processFragmentBatch(batch, ctx, progress);
                  
                  // Merge results
                  if (!ctx.results.empty()) {
                      std::lock_guard<std::mutex> lock(results_mutex);
                      allMappings.insert(
                          allMappings.end(),
                          std::make_move_iterator(ctx.results.begin()),
                          std::make_move_iterator(ctx.results.end())
                      );
                  }
              });
          }
          
          size_t fragThreads = std::max<size_t>(1, param.threads);
	  tf::Executor fragmentExecutor(fragThreads);
	  fragmentExecutor.run(batchTaskflow).wait();

          // Post-process and output results
          finalizeQueryResults(queryData, allMappings, progress, subsetMappings, subsetMappings_mutex, outstream);
      }

      // Common finalization for query results
      void finalizeQueryResults(QueryResults& queryData,
                              MappingResultsVector_t& allMappings,
                              std::shared_ptr<progress_meter::ProgressMeter> progress,
                              std::shared_ptr<std::unordered_map<seqno_t, MappingResultsVector_t>> subsetMappings,
                              std::shared_ptr<std::mutex> subsetMappings_mutex,
                              std::shared_ptr<std::ofstream> outstream) {
          
          // Boundary sanity check
          InputSeqProgContainer input(queryData.sequence, queryData.queryName, queryData.seqId, *progress);
          OutputHandler::mappingBoundarySanityCheck(&input, allMappings, *idManager);
          
          // Filter results with timing for large sequences
          auto filter_start = std::chrono::high_resolution_clock::now();
          
          // Log filtering info for large sequences
          if (queryData.queryLen > 1000000) {
              std::cerr << "[wfmash::mashmap] Filtering " << allMappings.size() 
                        << " mappings for " << queryData.queryName 
                        << " (" << queryData.queryLen << "bp)";
              if (param.scaffold_gap > 0) {
                  std::cerr << " [scaffold filtering enabled]";
              }
              std::cerr << "..." << std::endl;
          }
          
          auto filteredResult = filterSubsetMappings(allMappings, *progress, queryData.seqId, queryData.queryLen,
                                                   nullptr, nullptr, nullptr);
          
          auto filter_end = std::chrono::high_resolution_clock::now();
          auto filter_duration = std::chrono::duration_cast<std::chrono::milliseconds>(filter_end - filter_start);
          
          // Report long filtering times
          if (filter_duration.count() > 1000) {
              std::cerr << "[wfmash::mashmap] Filtering " << queryData.queryName 
                        << " took " << std::fixed << std::setprecision(1) 
                        << filter_duration.count() / 1000.0 << "s" << std::endl;
          }
          
          auto& mappings = param.mergeMappings && param.split ?
                filteredResult.mergedMappings : filteredResult.nonMergedMappings;
          auto& chainInfo = param.mergeMappings && param.split ?
                filteredResult.mergedChainInfo : filteredResult.nonMergedChainInfo;
          
          // Save or output results
          if (param.filterMode == filter::ONETOONE) {
              std::lock_guard<std::mutex> lock(*subsetMappings_mutex);
              (*subsetMappings)[queryData.seqId].insert(
                  (*subsetMappings)[queryData.seqId].end(),
                  std::make_move_iterator(mappings.begin()),
                  std::make_move_iterator(mappings.end())
              );
          } else {
              static std::mutex output_mutex;
              std::lock_guard<std::mutex> lock(output_mutex);
              
              if (param.mergeMappings && param.split && !chainInfo.empty()) {
                  OutputHandler::reportReadMappings(mappings, chainInfo, queryData.queryName, *outstream, 
                                                  *idManager, param, processMappingResults, queryData.queryLen);
              } else {
                  OutputHandler::reportReadMappings(mappings, queryData.queryName, *outstream, 
                                                  *idManager, param, processMappingResults, queryData.queryLen);
              }
              outstream->flush();
          }
      }

      void processFinalOneToOneResults(tf::Executor& executor,
                                     std::unordered_map<seqno_t, MappingResultsVector_t>& combinedMappings) {
          std::cerr << "[wfmash::mashmap] Processing final one-to-one filtering" << std::endl;
          
          size_t total_mappings = 0;
          for (auto& [querySeqId, mappings] : combinedMappings) {
              total_mappings += mappings.size();
          }
          std::cerr << "[wfmash::mashmap] Processing " << total_mappings 
                    << " mappings from " << combinedMappings.size() 
                    << " queries for one-to-one filtering" << std::endl;
          
          // Structure to track query ID with each mapping for target-based filtering
          struct MappingWithQuery {
              MappingResult mapping;
              seqno_t querySeqId;
          };
          
          std::unordered_map<seqno_t, std::vector<MappingWithQuery>> targetMappings;
          for (auto& [querySeqId, mappings] : combinedMappings) {
              for (auto& mapping : mappings) {
                  targetMappings[mapping.refSeqId].push_back({mapping, querySeqId});
              }
          }
          
          auto filterProgress = std::make_shared<progress_meter::ProgressMeter>(
              targetMappings.size(), 
              "[wfmash::mashmap] One-to-one reference filtering",
              param.use_progress_bar);
          
          // After filtering, group back by query
          std::unordered_map<seqno_t, MappingResultsVector_t> finalMappings;
          
          for (auto& [targetSeqId, mappingsWithQuery] : targetMappings) {
              // Extract just the mappings for filtering
              MappingResultsVector_t mappings;
              mappings.reserve(mappingsWithQuery.size());
              for (const auto& mwq : mappingsWithQuery) {
                  mappings.push_back(mwq.mapping);
              }
              
              MappingResultsVector_t filteredMappings;
              FilterUtils::filterByGroup(mappings, filteredMappings, param.numMappingsForSegment - 1, 
                           true, *idManager, param, *filterProgress);
              
              // Match filtered mappings back to their queries
              for (const auto& filteredMapping : filteredMappings) {
                  // Find which query this mapping came from
                  for (const auto& mwq : mappingsWithQuery) {
                      if (mwq.mapping.refSeqId == filteredMapping.refSeqId &&
                          mwq.mapping.refStartPos == filteredMapping.refStartPos &&
                          mwq.mapping.queryStartPos == filteredMapping.queryStartPos) {
                          finalMappings[mwq.querySeqId].push_back(filteredMapping);
                          break;
                      }
                  }
              }
              
              filterProgress->increment(1);
          }
          filterProgress->finish();
          
          std::ofstream outstrm(param.outFileName);
          size_t final_mapping_count = 0;
          for (auto& [querySeqId, mappings] : finalMappings) {
              std::string queryName = idManager->getSequenceName(querySeqId);
              offset_t queryLen = idManager->getSequenceLength(querySeqId);
              OutputHandler::reportReadMappings(mappings, queryName, outstrm, *idManager, param, processMappingResults, queryLen);
              final_mapping_count += mappings.size();
              if (final_mapping_count % 10000 == 0) {
                  outstrm.flush();
              }
          }
          outstrm.flush();
          outstrm.close();
          std::cerr << "[wfmash::mashmap] Wrote " << final_mapping_count 
                    << " mappings after one-to-one filtering" << std::endl;
      }

      /**
       * @brief Map single query fragment through L1 and L2 stages
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

          if (param.stage1_topANI_filter)
          {
            std::make_heap(l1_begin, l1_end, L1_locus_intersection_cmp);
          }
          doL2Mapping(Q, l1_begin, l1_end, l2Mappings);

          l1_begin = l1_end;
        }

        std::sort(l2Mappings.begin(), l2Mappings.end(), [](const auto& a, const auto& b) 
            { return std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos); });

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

      /**
       * @brief Execute L1 mapping stage
       */
      template <typename Q_Info, typename IPVec, typename L1Vec>
      void doL1Mapping(Q_Info &Q, IPVec& intervalPoints, L1Vec& l1Mappings)
      {
        //1. Compute the minmers
        CoreMapping::getSeedHits(Q, param);

        //Catch all NNNNNN case
        if (Q.sketchSize == 0 || Q.kmerComplexity < param.kmerComplexityThreshold) {
          return;
        }

        //2. Compute windows and sort
        CoreMapping::getSeedIntervalPoints(Q, intervalPoints, refSketch, *idManager, param);

        //3. Compute L1 windows
        int minimumHits = (Q.len == cached_segment_length) ?
            cached_minimum_hits :
            std::max(param.minimum_hits, Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity, skch::fixed::confidence_interval));

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
          CoreMapping::computeL1CandidateRegions(Q, ip_begin, ip_end, minimumHits, param, sketchCutoffs, l1Mappings);

          ip_begin = ip_end;
        }
      }

      /**
       * @brief Execute L2 mapping stage
       */
      template <typename Q_Info, typename L1_Iter, typename VecOut>
      void doL2Mapping(Q_Info &Q, L1_Iter l1_begin, L1_Iter l1_end, VecOut &l2Mappings)
      {
        std::vector<L2_mapLocus_t> l2_vec;
        double bestJaccardNumerator = 0;
        auto loc_iterator = l1_begin;
        
        while (loc_iterator != l1_end)
        {
          L1_candidateLocus_t& candidateLocus = *loc_iterator;

          if (param.stage1_topANI_filter)
          {
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
          CoreMapping::computeL2MappedRegions(Q, candidateLocus, l2_vec, refSketch, param);

          for (auto& l2 : l2_vec) 
          {
            float mash_dist = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);
            float nucIdentity = (1 - mash_dist);
            float nucIdentityUpperBound = 1 - Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, skch::fixed::confidence_interval);

            const auto& ref = this->idManager->getContigInfo(l2.seqId);
            if((param.keep_low_pct_id && nucIdentityUpperBound >= param.percentageIdentity)
                || nucIdentity >= param.percentageIdentity)
            {
              bestJaccardNumerator = std::max<double>(bestJaccardNumerator, l2.sharedSketchSize);

              MappingResult res;
              {
                // Fields that exist in compact struct
                res.refSeqId = l2.seqId;
                res.refStartPos = l2.meanOptimalPos;
                res.queryStartPos = 0;
                res.blockLength = Q.len; // mapping covers full fragment length
                res.conservedSketches = l2.sharedSketchSize;
                res.n_merged = 1; // Single mapping, not merged
                res.setNucIdentity(nucIdentity);
                res.setKmerComplexity(Q.kmerComplexity);
                res.setStrand(l2.strand);
                res.setDiscard(false);
                res.setOverlapped(false);
              } 
              l2Mappings.push_back(res);
            }
          }

          if (param.stage1_topANI_filter) 
          {
            std::pop_heap(l1_begin, l1_end, L1_locus_intersection_cmp); 
            l1_end--;
          }
          else 
          {
            loc_iterator++;
          }
        }
      }

      /**
       * @brief Result structure for filtered mappings with chain info
       */
      struct FilteredMappingsResult {
          MappingResultsVector_t nonMergedMappings;
          MappingResultsVector_t mergedMappings;
          ChainInfoVector_t nonMergedChainInfo;
          ChainInfoVector_t mergedChainInfo;
      };

      /**
       * @brief Filter mappings for a subset before aggregation
       */
      FilteredMappingsResult filterSubsetMappings(
          MappingResultsVector_t& mappings, 
          progress_meter::ProgressMeter& progress,
          seqno_t querySeqId,
          offset_t queryLen,
          std::shared_ptr<progress_meter::ProgressMeter> scaffold_progress,
          std::shared_ptr<std::atomic<size_t>> scaffold_total_work,
          std::shared_ptr<std::atomic<size_t>> scaffold_completed_work)
      {
          FilteredMappingsResult result;
          
          if (mappings.empty()) return result;

          MappingResultsVector_t rawMappings = mappings;
          
          // Pass context to the merge function - now returns mappings with chain info
          auto mappingsWithChains = FilterUtils::mergeMappingsInRangeWithChains(mappings, param.chain_gap, param, progress, querySeqId, queryLen);
          auto& maximallyMergedMappings = mappingsWithChains.mappings;
          auto& chainInfo = mappingsWithChains.chainInfo;

          if (param.mergeMappings && param.split) {
              // Pass context (queryLen) to weak mapping filter
              FilterUtils::filterWeakMappings(maximallyMergedMappings, 
                  std::floor(param.block_length / param.windowLength), param, *idManager, queryLen);
              
              if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                  MappingResultsVector_t groupFilteredMappings;
                  FilterUtils::filterByGroup(maximallyMergedMappings, groupFilteredMappings, 
                      param.numMappingsForSegment - 1, false, *idManager, param, progress);
                  maximallyMergedMappings = std::move(groupFilteredMappings);
              }

              if (param.filterLengthMismatches) {
                  FilterUtils::filterFalseHighIdentity(maximallyMergedMappings, param);
              }
              
              FilterUtils::sparsifyMappings(maximallyMergedMappings, param);
              // Pass context to scaffold filter
              FilterUtils::filterByScaffolds(maximallyMergedMappings, rawMappings, param, *idManager, progress, querySeqId, queryLen,
                                           scaffold_progress, scaffold_total_work, scaffold_completed_work);
          } else {
              // Same for non-merged path
              if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                  MappingResultsVector_t filteredMappings;
                  FilterUtils::filterByGroup(mappings, filteredMappings, 
                      param.numMappingsForSegment - 1, false, *idManager, param, progress);
                  mappings = std::move(filteredMappings);
              }
              FilterUtils::filterByScaffolds(mappings, rawMappings, param, *idManager, progress, querySeqId, queryLen,
                                           scaffold_progress, scaffold_total_work, scaffold_completed_work);
          }

          // Return results with chain info
          result.nonMergedMappings = std::move(mappings);
          result.mergedMappings = std::move(maximallyMergedMappings);
          
          // For non-merged mappings, each is its own chain
          result.nonMergedChainInfo.resize(result.nonMergedMappings.size());
          for (size_t i = 0; i < result.nonMergedMappings.size(); ++i) {
              result.nonMergedChainInfo[i] = {static_cast<uint32_t>(i), 1, 1};
          }
          result.mergedChainInfo = std::move(chainInfo);
          
          return result;
      }

    public:
      /**
       * @brief Utility function to save L2 results to vector
       */
      static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
      {
        OutputHandler::insertL2ResultsToVec(v, reportedL2Result);
      }
  };

  // Initialize thread-local storage
  thread_local std::unique_ptr<Map::FragmentContext> Map::tls_fragmentContext;
}

#endif // SKETCH_MAP_HPP
