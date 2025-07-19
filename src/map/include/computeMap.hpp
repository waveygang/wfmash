/**
 * @file    computeMap.hpp
 * @brief   Main Map class that implements the sequence mapping logic
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
#include <mutex>
#include <thread>
#include <sstream>
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
#include "map/include/compressedMapping.hpp"

//External includes
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "gsl/gsl_randist.h"

namespace skch
{
  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages
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

    void processFragment(const FragmentData& fragment, 
                         std::vector<IntervalPoint>& intervalPoints,
                         std::vector<L1_candidateLocus_t>& l1Mappings,
                         MappingResultsVector_t& l2Mappings,
                         QueryMetaData<MinVec_Type>& Q,
                         std::vector<MappingResult>& thread_local_results) {
        
        intervalPoints.clear();
        l1Mappings.clear();
        thread_local_results.clear();

        Q.seq = const_cast<char*>(fragment.seq);
        Q.len = fragment.len;
        Q.fullLen = fragment.fullLen;
        Q.seqId = fragment.seqId;
        Q.seqName = fragment.seqName;
        Q.refGroup = fragment.refGroup;

        // Always use compressed storage internally for memory efficiency
        CompressedMappingStore compressedStore(fragment.seqId, fragment.fullLen);
        mapSingleQueryFrag(Q, intervalPoints, l1Mappings, compressedStore);
        
        // Retrieve and update mappings with query positions
        auto mappings = compressedStore.getAllMappings();
        for (auto& e : mappings) {
            e.queryLen = fragment.fullLen;
            e.queryStartPos = fragment.fragmentIndex * param.segLength;
            e.queryEndPos = e.queryStartPos + fragment.len;
        }
        
        thread_local_results = std::move(mappings);
        
        // Note: l2Mappings parameter is now unused but kept for interface compatibility
        l2Mappings.clear();
        
        if (fragment.output) {
            fragment.output->progress.increment(fragment.len);
        }
    }
      
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
        cached_segment_length(p.segLength),
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
          tf::Taskflow taskflow;

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
              
              auto subset_flow = std::make_shared<tf::Taskflow>();

              uint64_t subset_query_length = 0;
              for (const auto& queryName : querySequenceNames) {
                  subset_query_length += idManager->getSequenceLength(idManager->getSequenceId(queryName));
              }

              auto progress = std::make_shared<progress_meter::ProgressMeter>(
                  subset_query_length,
                  "[wfmash::mashmap] mapping",
                  param.use_progress_bar
                  );

              auto buildIndex_task = subset_flow->emplace([this, target_subset=target_subset, subset_idx, total_subsets=target_subsets.size()]() {
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
                      size_t seq_count = refSketch->getSequenceCount();
                  } else {
                      auto sketch_index_progress = std::make_shared<progress_meter::ProgressMeter>(
                          100,
                          "[wfmash::mashmap] indexing",
                          param.use_progress_bar);
                  
                      if (!sketch_index_progress->is_finished.load()) {
                          sketch_index_progress->reset_timer();
                      }
                      
                      refSketch = new skch::Sketch(param, *idManager, target_subset, nullptr, sketch_index_progress);
                      std::cerr << "[wfmash::mashmap] building index data structures..." << std::endl;
                  }
              }).name("build_index_" + std::to_string(subset_idx));

              auto subsetMappings = std::make_shared<std::unordered_map<seqno_t, MappingResultsVector_t>>();
              auto subsetMappings_mutex = std::make_shared<std::mutex>();

              auto outstream = std::make_shared<std::ofstream>();
              auto outstream_mutex = std::make_shared<std::mutex>();
              
              if (param.filterMode != filter::ONETOONE) {
                  bool append = subset_idx > 0;
                  outstream->open(param.outFileName, append ? std::ios::app : std::ios::out);
                  if (!outstream->is_open()) {
                      std::cerr << "Error: Could not open output file for writing: " << param.outFileName << std::endl;
                      exit(1);
                  }
              }

              auto processQueries_task = subset_flow->emplace([this, progress, subsetMappings, subsetMappings_mutex, 
                                                           outstream, outstream_mutex, subset_flow](tf::Subflow& sf) {
                  const auto& fileName = param.querySequences[0];
    
                  faidx_meta_t* query_meta = faidx_meta_load(fileName.c_str(), FAI_FASTA, FAI_CREATE);
                  if (!query_meta) {
                      std::cerr << "Error: Failed to load query FASTA index: " << fileName << std::endl;
                      exit(1);
                  }
                  
                  auto getThreadLocalReader = [](faidx_meta_t* meta) -> faidx_reader_t* {
                      thread_local faidx_reader_t* reader = nullptr;
                      if (reader == nullptr) {
                          reader = faidx_reader_create(meta);
                          if (!reader) {
                              throw std::runtime_error("Failed to create thread-local reader");
                          }
                      }
                      return reader;
                  };
                  
                  for (const auto& queryName : querySequenceNames) {
                      faidx_reader_t* reader = getThreadLocalReader(query_meta);
                      
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
                      
                      std::string sequence(seq_data, seq_len);
                      free(seq_data);
            
                      auto query_task = sf.emplace([&, queryName, sequence](tf::Subflow& query_sf) {
                          seqno_t seqId = idManager->getSequenceId(queryName);
                          auto input = std::make_shared<InputSeqProgContainer>(
                              sequence, queryName, seqId, *progress);
                          auto output = std::make_shared<QueryMappingOutput>(
                              queryName, MappingResultsVector_t{}, MappingResultsVector_t{}, *progress);
            
                          int refGroup = idManager->getRefGroup(seqId);
                          int noOverlapFragmentCount = input->len / param.segLength;

                          std::mutex results_mutex;

                          // Regular fragments
                          for(int i = 0; i < noOverlapFragmentCount; i++) {
                              query_sf.emplace([&, i]() {
                                  std::vector<MappingResult> all_fragment_results;
                                  auto fragment = std::make_shared<FragmentData>(
                                      &(sequence)[0u] + i * param.segLength,
                                      static_cast<int>(param.segLength),
                                      static_cast<int>(input->len),
                                      seqId,
                                      queryName,
                                      refGroup,
                                      i,
                                      output
                                  );

                                  static std::once_flag first_fragment;
                                  std::call_once(first_fragment, [&output]() {
                                      output->progress.reset_timer();
                                  });
                                  
                                  std::vector<IntervalPoint> intervalPoints;
                                  std::vector<L1_candidateLocus_t> l1Mappings;
                                  MappingResultsVector_t l2Mappings;
                                  QueryMetaData<MinVec_Type> Q;
                                  processFragment(*fragment, intervalPoints, l1Mappings, l2Mappings, Q, all_fragment_results);
                                  
                                  if (!all_fragment_results.empty()) {
                                      std::lock_guard<std::mutex> lock(results_mutex);
                                      output->results.addMappings(all_fragment_results);
                                  }
                              });
                          }

                          // Handle final fragment if needed
                          if (noOverlapFragmentCount >= 1 && input->len % param.segLength != 0) {
                              query_sf.emplace([&]() {
                                  std::vector<MappingResult> all_fragment_results;
                                  auto fragment = std::make_shared<FragmentData>(
                                      &(sequence)[0u] + input->len - param.segLength,
                                      static_cast<int>(param.segLength),
                                      static_cast<int>(input->len),
                                      seqId,
                                      queryName,
                                      refGroup,
                                      noOverlapFragmentCount,
                                      output
                                  );

                                  std::vector<IntervalPoint> intervalPoints;
                                  std::vector<L1_candidateLocus_t> l1Mappings;
                                  MappingResultsVector_t l2Mappings;
                                  QueryMetaData<MinVec_Type> Q;
                                  processFragment(*fragment, intervalPoints, l1Mappings, l2Mappings, Q, all_fragment_results);
                                  
                                  if (!all_fragment_results.empty()) {
                                      std::lock_guard<std::mutex> lock(results_mutex);
                                      output->results.addMappings(all_fragment_results);
                                  }
                              });
                          }

                          query_sf.join();

                          // Get uncompressed results for processing
                          auto uncompressed_results = output->results.getAllMappings();
                          OutputHandler::mappingBoundarySanityCheck(input.get(), uncompressed_results, *idManager);
                          auto [nonMergedMappings, mergedMappings] = 
                              filterSubsetMappings(uncompressed_results, output->progress);

                          auto& mappings = param.mergeMappings && param.split ?
                                mergedMappings : nonMergedMappings;

                          if (param.filterMode == filter::ONETOONE) {
                              std::lock_guard<std::mutex> lock(*subsetMappings_mutex);
                              (*subsetMappings)[seqId].insert(
                                  (*subsetMappings)[seqId].end(),
                                  std::make_move_iterator(mappings.begin()),
                                  std::make_move_iterator(mappings.end())
                              );
                          } else {
                              std::lock_guard<std::mutex> lock(*outstream_mutex);
                              OutputHandler::reportReadMappings(mappings, queryName, *outstream, *idManager, param, processMappingResults);
                          }
                      }).name("query_" + queryName);
                  }
                  
                  faidx_meta_destroy(query_meta);
              }).name("process_queries");

              auto merge_task = subset_flow->emplace([this,
                                                    subsetMappings,
                                                    &combinedMappings,
                                                    &combinedMappings_mutex]() {
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
              }).name("merge_results");

              auto cleanup_task = subset_flow->emplace([this, outstream]() {
                  if (param.filterMode != filter::ONETOONE && outstream->is_open()) {
                      outstream->close();
                  }
                  
                  delete refSketch;
                  refSketch = nullptr;
              }).name("cleanup");

              buildIndex_task.precede(processQueries_task);
              processQueries_task.precede(merge_task);
              merge_task.precede(cleanup_task);

              executor.run(*subset_flow).wait();
        
              progress->finish();
          }

          if (exit_after_indices) {
              std::cerr << "[wfmash::mashmap] All indices created successfully. Exiting." << std::endl;
              exit(0);
          }

          // Final results processing for ONETOONE mode
          if (param.filterMode == filter::ONETOONE && !exit_after_indices) {
              tf::Taskflow final_flow;
              std::cerr << "[wfmash::mashmap] Processing final one-to-one filtering" << std::endl;
              
              auto final_task = final_flow.emplace([&]() {
                  size_t total_mappings = 0;
                  for (auto& [querySeqId, mappings] : combinedMappings) {
                      total_mappings += mappings.size();
                  }
                  std::cerr << "[wfmash::mashmap] Processing " << total_mappings 
                            << " mappings from " << combinedMappings.size() 
                            << " queries for one-to-one filtering" << std::endl;
                  
                  std::unordered_map<seqno_t, MappingResultsVector_t> targetMappings;
                  for (auto& [querySeqId, mappings] : combinedMappings) {
                      for (auto& mapping : mappings) {
                          targetMappings[mapping.refSeqId].push_back(mapping);
                      }
                  }
                  
                  auto filterProgress = std::make_shared<progress_meter::ProgressMeter>(
                      targetMappings.size(), 
                      "[wfmash::mashmap] One-to-one reference filtering",
                      param.use_progress_bar);
                  
                  std::unordered_map<seqno_t, MappingResultsVector_t> finalMappings;
                  for (auto& [targetSeqId, mappings] : targetMappings) {
                      MappingResultsVector_t filteredMappings;
                      FilterUtils::filterByGroup(mappings, filteredMappings, param.numMappingsForSegment - 1, 
                                   true, *idManager, param, *filterProgress);
                      
                      for (auto& mapping : filteredMappings) {
                          finalMappings[mapping.querySeqId].push_back(mapping);
                      }
                      filterProgress->increment(1);
                  }
                  filterProgress->finish();
                  
                  std::ofstream outstrm(param.outFileName);
                  size_t final_mapping_count = 0;
                  for (auto& [querySeqId, mappings] : finalMappings) {
                      std::string queryName = idManager->getSequenceName(querySeqId);
                      OutputHandler::reportReadMappings(mappings, queryName, outstrm, *idManager, param, processMappingResults);
                      final_mapping_count += mappings.size();
                  }
                  
                  std::cerr << "[wfmash::mashmap] Wrote " << final_mapping_count 
                            << " mappings after one-to-one filtering" << std::endl;
              });
              
              executor.run(final_flow).wait();
          }
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

      /**
       * @brief Map single query fragment through L1 and L2 stages with compressed storage
       * @note This is the memory-efficient path using CompressedMappingStore
       */
      template<typename Q_Info, typename IPVec, typename L1Vec>
      void mapSingleQueryFrag(Q_Info &Q, IPVec& intervalPoints, L1Vec& l1Mappings, CompressedMappingStore &compressedStore)
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

        // Temporary vector for L2 mappings
        MappingResultsVector_t tempL2Mappings;
        
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
          
          tempL2Mappings.clear();
          doL2Mapping(Q, l1_begin, l1_end, tempL2Mappings);
          
          // Add mappings to compressed store immediately
          compressedStore.addMappings(tempL2Mappings);

          l1_begin = l1_end;
        }
        
        // Sort all mappings in the compressed store
        compressedStore.sortByRefPos();
        
        // Assign chain IDs to compressed mappings
        int32_t chain_id = maxChainIdSeen.fetch_add(1, std::memory_order_relaxed);
        int32_t chain_length = compressedStore.size();
        int32_t chain_pos = 1;
        
        // We need to update chain info, but compressed store doesn't support direct modification
        // So we'll retrieve, update, and re-store
        auto allMappings = compressedStore.getAllMappings();
        compressedStore.clear();
        
        for (auto& mapping : allMappings) {
            mapping.chain_id = chain_id;
            mapping.chain_length = chain_length;
            mapping.chain_pos = chain_pos++;
            compressedStore.addMapping(mapping);
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
            l1_end--;
          }
          else 
          {
            loc_iterator++;
          }
        }
      }

      /**
       * @brief Filter mappings for a subset before aggregation
       */
      std::pair<MappingResultsVector_t, MappingResultsVector_t> filterSubsetMappings(
          MappingResultsVector_t& mappings, 
          progress_meter::ProgressMeter& progress) 
      {
          if (mappings.empty()) return {MappingResultsVector_t(), MappingResultsVector_t()};

          MappingResultsVector_t rawMappings = mappings;
          
          auto maximallyMergedMappings = FilterUtils::mergeMappingsInRange(mappings, param.chain_gap, param, progress);

          if (param.mergeMappings && param.split) {
              FilterUtils::filterWeakMappings(maximallyMergedMappings, 
                  std::floor(param.block_length / param.segLength), param, *idManager);
              
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
              FilterUtils::filterByScaffolds(maximallyMergedMappings, rawMappings, param, *idManager, progress);
          } else {
              // Filter non-merged mappings
              if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                  MappingResultsVector_t filteredMappings;
                  FilterUtils::filterByGroup(mappings, filteredMappings, 
                      param.numMappingsForSegment - 1, false, *idManager, param, progress);
                  mappings = std::move(filteredMappings);
              }
              FilterUtils::filterByScaffolds(mappings, rawMappings, param, *idManager, progress);
          }

          // Build dense chain ID mapping
          std::unordered_map<offset_t, offset_t> id_map;
          offset_t next_id = 0;
          
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

          offset_t base_id = maxChainIdSeen.fetch_add(id_map.size(), std::memory_order_relaxed);
          
          for (auto& mapping : mappings) {
              mapping.splitMappingId = id_map[mapping.splitMappingId] + base_id;
          }
          for (auto& mapping : maximallyMergedMappings) {
              mapping.splitMappingId = id_map[mapping.splitMappingId] + base_id;
          }

          return {std::move(mappings), std::move(maximallyMergedMappings)};
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
}

#endif // SKETCH_MAP_HPP