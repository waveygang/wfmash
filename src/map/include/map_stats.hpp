/**
 * @file    map_stats.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef MAP_STATS_HPP 
#define MAP_STATS_HPP

#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>
#include <unordered_set>
#include <iomanip>
#include <mutex>
#include <thread>
#include <atomic>

#ifdef USE_BOOST
    #include <boost/math/distributions/binomial.hpp>
    using namespace::boost::math;
#else
    #include <gsl/gsl_cdf.h>
#endif

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/sequenceIds.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/streamingMinHash.hpp"
#include "common/taskflow/taskflow.hpp"

//External includes
#include "common/murmur3.h"
#include "common/kseq.h"
#include "common/prettyprint.hpp"
#include "common/seqiter.hpp"
#include "common/faigz.h"

namespace skch
{
  /**
   * @namespace skch::Stat
   * @brief     Implements utility functions that involve statistical computation
   */
  namespace Stat
  {

    /**
     * @brief         jaccard estimate to mash distance
     * @param[in] j   jaccard estimate
     * @param[in] k   kmer size 
     * @return        mash distance [0.0 - 1.0]
     */
    inline float j2md(float j, int k)
    {
      if(j == 0)
        return 1.0; //jaccard estimate 0 -> 1.0 mash distance

      if(j == 1)
        return 0.0; //jaccard estimate 1 -> 0.0 mash distance

      float mash_dist = 1 - std::pow(2*j / (1+j), 1.0/k);
      return mash_dist;
    }

    /**
     * @brief         mash distance to jaccard estimate
     * @param[in] d   mash distance [0.0 - 1.0]
     * @param[in] k   kmer size 
     * @return        jaccard estimate 
     */
    inline float md2j(float d, int k)
    {
      float sim = 1 - d;
      float jaccard = std::pow(sim, k) / (2 - std::pow(sim, k));
      return jaccard;
    }

    /**
     * @brief               Given a distance d, compute the lower bound on d within required confidence interval 
     * @details             If a given match has distance d in the L2 stage, we compare its lower distance bound 
     *                      against the assumed cutoff to decide its significance. This makes the mapping algorithm
     *                      more sensitive to true alignments
     * @param[in]   d       calculated mash distance
     * @param[in]   s       sketch size
     * @param[in]   k       kmer size
     * @param[in]   ci      confidence interval [0-1], example 0.9 implies 90% confidence interval
     * @return              computed lower bound on d within 'ci' confidence interval
     */
    inline float md_lower_bound(float d, int s, int k, float ci)
    {
      //One side interval probability
      float q2 = (1.0 - ci)/2;

      //Computing count of sketches using confidence interval
#ifdef USE_BOOST
      
      //Inverse binomial 
      int x = quantile(complement(binomial(s, md2j(d,k)), q2));

#else   
      //GSL 
      int x = std::max( int(ceil(s * md2j(d,k))), 1 );    //Begin search from jaccard * s
      while(x <= s)
      {
        //probability of having x or more shared sketches
        double cdf_complement = gsl_cdf_binomial_Q(x-1, md2j(d,k), s);

        if (cdf_complement < q2)
        {
          x--;  //Last guess was right
          break;
        }

        x++;
      }
#endif

      float jaccard = float(x) / s;
      float low_d = j2md(jaccard, k);
      return low_d; 
    }

    /**
     * @brief                 Estimate minimum number of shared sketches to achieve the desired identity
     * @param[in] s           sketch size
     * @param[in] k           kmer size
     * @param[in] identity    percentage identity [0-1]
     * @return                minimum count of hits
     */
    inline int estimateMinimumHits(int s, int k, float perc_identity)
    {
      //Compute the estimate
      float mash_dist = 1.0 - perc_identity;
      float jaccard = md2j(mash_dist, k);

      //function to convert jaccard to min hits
      //Atleast these many minimizers should match for achieving the required jaccard identity
      int minimumSharedMinimizers = ceil (1.0 * s * jaccard); 

      return minimumSharedMinimizers;
    }

    /**
     * @brief                 Estimate minimum number of shared sketches 
     *                        s.t. upper bound identity is >= desired identity
     *                        Upper bound is computed using the 90% confidence interval
     * @param[in] s           sketch size
     * @param[in] k           kmer size
     * @param[in] identity    percentage identity [0-1]
     * @return                count of min. shared minimizers
     */
    inline int estimateMinimumHitsRelaxed(int s, int k, float perc_identity, float confidence_interval)
    {
      // The desired value has be between [0, min  s.t. identity >= perc_identity]
      auto searchRange = std::pair<int, int>( estimateMinimumHits(s, k, perc_identity) , 0);

      int minimumSharedMinimizers_relaxed = searchRange.first;

      for(int i = searchRange.first ; i >= searchRange.second; i--)
      {
        float jaccard = 1.0 * i/s;
        float d = j2md(jaccard, k);

        float d_lower = md_lower_bound(d, s, k, confidence_interval);

        //Upper bound identity
        float id_upper = 1.0 - d_lower;

        //Check if it satisfies the criteria
        if(id_upper >= perc_identity)
          minimumSharedMinimizers_relaxed = i;
        else
          break;    //Stop the search
      }

      return minimumSharedMinimizers_relaxed;
    }

    /**
     * @brief                     calculate p-value for a given alignment identity, sketch size..
     * @param[in] s               sketch size
     * @param[in] k               kmer size
     * @param[in] alphabetSize    alphabet size
     * @param[in] identity        mapping identity cut-off
     * @param[in] lengthQuery     query length
     * @param[in] lengthReference reference length
     * @return                    p-value
     */
    inline double estimate_pvalue (
        int s, 
        int k, 
        int alphabetSize, 
        float identity,
        int64_t lengthQuery, 
        uint64_t lengthReference, 
        float confidence_interval)
    {
      //total space size of k-mers
      double kmerSpace = pow(alphabetSize, k);

      //probability of a kmer match by random in |query| sized sequence 
      double pX, pY; 
      pX = pY = 1. / (1. + kmerSpace / lengthQuery);

      //Jaccard similarity of two random given sequences
      double r = pX * pY / (pX + pY - pX * pY);

      int x = estimateMinimumHitsRelaxed(s, k, identity, confidence_interval);

      //P (x or more minimizers match)
      double cdf_complement;
      if(x == 0)
      {
        cdf_complement = 1.0;
      }
      else
      {
#ifdef USE_BOOST
      cdf_complement = cdf(complement(binomial(s, r), x-1));
#else
      cdf_complement =  gsl_cdf_binomial_Q(x-1, r, s);
#endif
      }

      double pVal = lengthReference * cdf_complement;

      return pVal;
    }

    /**
     * @brief                       calculate minimum window size for sketching that satisfies
     *                              the given p-value threshold
     * @param[in] pValue_cutoff     cut off p-value threshold
     * @param[in] confidence_interval confidence interval to relax jaccard cutoff for mapping
     * @param[in] k                 kmer size
     * @param[in] alphabetSize      alphabet size
     * @param[in] identity          mapping identity cut-off
     * @param[in] segmentLength     mashmap's internal minimum query sequence length
     * @param[in] lengthReference   reference length
     * @return                      optimal window size for sketching
     */
    inline int64_t recommendedSketchSize(
        double pValue_cutoff, 
        float confidence_interval,
        int k, 
        int alphabetSize,
        float identity,
        int64_t segmentLength, 
        uint64_t lengthReference)
    {
      int64_t lengthQuery = segmentLength - k;

      int optimalSketchSize;
      for (optimalSketchSize = 10; optimalSketchSize < lengthQuery; optimalSketchSize += 10) {
        //Compute pvalue
        double pVal = estimate_pvalue(optimalSketchSize, k, alphabetSize, identity, lengthQuery, lengthReference, confidence_interval);

        //Check if pvalue is <= cutoff
        if(pVal <= pValue_cutoff)
        {
          break;
        }
      }

      return optimalSketchSize;
    }

    /**
     * @brief Helper function to update a MinHash sketch with k-mers from a sequence
     * @param[in] seq The sequence string
     * @param[in,out] sketch Vector maintaining the MinHash sketch (smallest hashes)
     * @param[in] k K-mer size
     * @param[in] sketch_size Number of minimizers to retain
     */
    inline void update_minhash_sketch(const std::string& seq, std::vector<hash_t>& sketch, int k, int sketch_size) {
      if ((int)seq.length() < k) return;
      
      // Use streaming MinHash for efficiency
      StreamingMinHash minhash(sketch_size, 0); // windowSize not used for pure MinHash
      
      // Add existing sketch values
      for (hash_t h : sketch) {
        minhash.add(h);
      }
      
      char* seq_data = const_cast<char*>(seq.c_str());
      
      // Process all k-mers
      for (offset_t i = 0; i <= (offset_t)seq.length() - k; ++i) {
        hash_t hash = CommonFunc::getHash(seq_data + i, k);
        minhash.add(hash);
      }
      
      // Get updated sketch
      sketch = minhash.getSketch();
    }

    /**
     * @brief Finalize MinHash sketch by sorting to get the smallest hashes
     * @param[in,out] sketch Vector containing the MinHash sketch
     * @param[in] sketch_size Number of minimizers to retain
     */
    inline void finalize_minhash_sketch(std::vector<hash_t>& sketch, int sketch_size) {
      // Convert from max heap to sorted vector of smallest hashes
      std::sort_heap(sketch.begin(), sketch.end());
      
      // Keep only the smallest sketch_size elements
      if ((int)sketch.size() > sketch_size) {
        sketch.resize(sketch_size);
      }
      
      // Sort in ascending order for efficient comparison
      std::sort(sketch.begin(), sketch.end());
    }

    /**
     * @brief Estimate optimal identity threshold based on sequence groups
     * @param[in] params Mapping parameters
     * @param[in] idManager Sequence ID manager with grouping information
     * @return Estimated identity threshold [0,1]
     */
    inline double estimate_identity_for_groups(const skch::Parameters& params, skch::SequenceIdManager& idManager) {
      // Use k=21 to match Mash defaults for genome comparison
      // This provides good sensitivity while avoiding too many random matches
      const int estimation_k = 21;
      // Use configurable sketch size for ANI estimation
      const int estimation_sketch_size = params.ani_sketch_size;

      std::cerr << "[wfmash::auto-identity] Starting identity estimation with k=" << estimation_k 
                << ", sketch_size=" << estimation_sketch_size << std::endl;

      // Maps to maintain MinHash sketches per group
      std::map<int, std::vector<hash_t>> query_group_sketches;
      std::map<int, std::vector<hash_t>> target_group_sketches;
      std::mutex sketch_mutex;
      std::atomic<int> query_seq_count(0);
      std::atomic<int> target_seq_count(0);

      // Determine number of threads to use
      int num_threads = params.threads > 0 ? params.threads : std::thread::hardware_concurrency();
      if (num_threads == 0) num_threads = 1;
      
      std::cerr << "[wfmash::auto-identity] Using " << num_threads << " threads for k-mer extraction" << std::endl;
      
      // First, count total sequences to process
      size_t total_sequences = 0;
      for (const auto& [name, id] : idManager.getSequenceNameToIdMap()) {
        total_sequences++;
      }
      std::cerr << "[wfmash::auto-identity] Total sequences in idManager: " << total_sequences << std::endl;

      // Process sequences in parallel using GroupedStreamingMinHash
      GroupedStreamingMinHash query_grouped_sketch(estimation_sketch_size, 0);
      GroupedStreamingMinHash target_grouped_sketch(estimation_sketch_size, 0);
      
      // Collect sequence names to process (NOT the sequences themselves)
      struct SeqToProcess {
        std::string name;
        std::string file;
        bool is_query;
      };
      std::vector<SeqToProcess> sequences_to_process;
      
      // Get sequence names from idManager
      const auto& seqMap = idManager.getSequenceNameToIdMap();
      std::unordered_set<std::string> valid_names;
      for (const auto& [name, id] : seqMap) {
        valid_names.insert(name);
      }
      
      // Collect query sequence names
      for (const auto& file : params.querySequences) {
        faidx_meta_t* meta = faidx_meta_load(file.c_str(), FAI_FASTA, FAI_CREATE);
        if (!meta) {
          std::cerr << "[wfmash::auto-identity] Warning: Failed to load FASTA index: " << file << std::endl;
          continue;
        }
        
        int nseq = faidx_meta_nseq(meta);
        for (int i = 0; i < nseq; i++) {
          const char* seqname = faidx_meta_iseq(meta, i);
          if (seqname && valid_names.count(seqname) > 0) {
            sequences_to_process.push_back({seqname, file, true});
          }
        }
        faidx_meta_destroy(meta);
      }
      
      // Collect target sequence names
      for (const auto& file : params.refSequences) {
        faidx_meta_t* meta = faidx_meta_load(file.c_str(), FAI_FASTA, FAI_CREATE);
        if (!meta) {
          std::cerr << "[wfmash::auto-identity] Warning: Failed to load FASTA index: " << file << std::endl;
          continue;
        }
        
        int nseq = faidx_meta_nseq(meta);
        for (int i = 0; i < nseq; i++) {
          const char* seqname = faidx_meta_iseq(meta, i);
          if (seqname && valid_names.count(seqname) > 0) {
            sequences_to_process.push_back({seqname, file, false});
          }
        }
        faidx_meta_destroy(meta);
      }
      
      std::cerr << "[wfmash::auto-identity] Sequences to process: " << sequences_to_process.size() 
                << " (query: " << std::count_if(sequences_to_process.begin(), sequences_to_process.end(), 
                                                [](const SeqToProcess& s) { return s.is_query; })
                << ", target: " << std::count_if(sequences_to_process.begin(), sequences_to_process.end(),
                                                  [](const SeqToProcess& s) { return !s.is_query; })
                << ")" << std::endl;
      
      // Create progress meter for ANI sketching
      progress_meter::ProgressMeter ani_progress(
          sequences_to_process.size(),
          "[wfmash::auto-identity] sketching",
          params.use_progress_bar);
      
      // Use taskflow for parallel processing
      tf::Executor executor(num_threads);
      tf::Taskflow taskflow;
      
      // Process all files in parallel
      taskflow.emplace([&](tf::Subflow& sf) {
        // Load meta for each file once
        std::map<std::string, faidx_meta_t*> file_metas;
        for (const auto& file : params.querySequences) {
          if (file_metas.find(file) == file_metas.end()) {
            faidx_meta_t* meta = faidx_meta_load(file.c_str(), FAI_FASTA, FAI_CREATE);
            if (!meta) {
              std::cerr << "[wfmash::auto-identity] Error: Failed to load FASTA index: " << file << std::endl;
              continue;
            }
            file_metas[file] = meta;
          }
        }
        for (const auto& file : params.refSequences) {
          if (file_metas.find(file) == file_metas.end()) {
            faidx_meta_t* meta = faidx_meta_load(file.c_str(), FAI_FASTA, FAI_CREATE);
            if (!meta) {
              std::cerr << "[wfmash::auto-identity] Error: Failed to load FASTA index: " << file << std::endl;
              continue;
            }
            file_metas[file] = meta;
          }
        }
        
        // Function to get thread-local reader (exact pattern from computeMap)
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
        
        // Process each sequence in parallel
        for (const auto& seq_info : sequences_to_process) {
          sf.emplace([&, seq_info]() {
            try {
              auto meta_it = file_metas.find(seq_info.file);
              if (meta_it == file_metas.end()) return;
              
              faidx_meta_t* meta = meta_it->second;
              faidx_reader_t* reader = getThreadLocalReader(meta);
              
              // Get sequence info from idManager
              auto it = seqMap.find(seq_info.name);
              if (it == seqMap.end()) return;
              
              seqno_t seqId = it->second;
              int groupId = idManager.getRefGroup(seqId);
              
              hts_pos_t seq_len;
              hts_pos_t seq_total_len = faidx_meta_seq_len(meta, seq_info.name.c_str());
              if (seq_total_len <= 0) {
                std::cerr << "[wfmash::auto-identity] Warning: Sequence " << seq_info.name 
                          << " not found or empty, skipping" << std::endl;
                return;
              }
              
              char* seq_data = faidx_reader_fetch_seq(reader, seq_info.name.c_str(), 0, seq_total_len-1, &seq_len);
              if (!seq_data) {
                std::cerr << "[wfmash::auto-identity] Warning: Failed to fetch sequence " 
                          << seq_info.name << ", skipping" << std::endl;
                return;
              }
              
              // Process the sequence using GroupedStreamingMinHash
              if (seq_info.is_query) {
                query_grouped_sketch.processSequence(
                    seq_data,
                    seq_len,
                    seqId,
                    groupId,
                    estimation_k,
                    4,  // alphabetSize for DNA
                    nullptr);
                query_seq_count++;
              } else {
                target_grouped_sketch.processSequence(
                    seq_data,
                    seq_len,
                    seqId,
                    groupId,
                    estimation_k,
                    4,  // alphabetSize for DNA
                    nullptr);
                target_seq_count++;
              }
              
              free(seq_data);
              ani_progress.increment(1);
              
            } catch (const std::exception& e) {
              std::cerr << "[wfmash::auto-identity] Warning: Failed to process sequence '" 
                        << seq_info.name << "': " << e.what() << std::endl;
            }
          });
        }
        
        sf.join();
        
        // Clean up metas
        for (auto& [file, meta] : file_metas) {
          if (meta) faidx_meta_destroy(meta);
        }
      });
      
      // Execute the taskflow
      executor.run(taskflow).wait();
      
      ani_progress.finish();
      
      // Get the final sketches from GroupedStreamingMinHash
      auto query_sketches = query_grouped_sketch.getAllGroupSketches();
      auto target_sketches = target_grouped_sketch.getAllGroupSketches();
      
      // Convert to the expected format
      for (const auto& [groupId, sketch] : query_sketches) {
        query_group_sketches[groupId] = sketch;
      }
      for (const auto& [groupId, sketch] : target_sketches) {
        target_group_sketches[groupId] = sketch;
      }

      std::cerr << "[wfmash::auto-identity] Processed " << query_seq_count 
                << " query sequences into " << query_group_sketches.size() << " groups" << std::endl;
      std::cerr << "[wfmash::auto-identity] Processed " << target_seq_count 
                << " target sequences into " << target_group_sketches.size() << " groups" << std::endl;

      if (query_seq_count == 0 || target_seq_count == 0) {
        std::cerr << "[wfmash::auto-identity] Warning: No sequences found in FAI files. " 
                  << "Make sure .fai index files exist (run 'samtools faidx' on your FASTA files)" << std::endl;
        return skch::fixed::percentage_identity;
      }

      // Finalize sketches for each group
      std::cerr << "[wfmash::auto-identity] Finalizing MinHash sketches..." << std::endl;
      
      for (auto& [groupId, sketch] : query_group_sketches) {
        finalize_minhash_sketch(sketch, estimation_sketch_size);
      }
      
      for (auto& [groupId, sketch] : target_group_sketches) {
        finalize_minhash_sketch(sketch, estimation_sketch_size);
      }

      std::vector<double> all_anis;
      int comparison_count = 0;

      for (const auto& q_pair : query_group_sketches) {
        for (const auto& t_pair : target_group_sketches) {
          // Self-comparison logic for all-vs-all mode
          bool is_self_mode = (&params.querySequences == &params.refSequences);
          if (is_self_mode && q_pair.first > t_pair.first) {
            continue; // Avoid redundant (B vs A) check if (A vs B) is done
          }
          
          const auto& q_sketch = q_pair.second;
          const auto& t_sketch = t_pair.second;

          if (q_sketch.empty() || t_sketch.empty()) continue;

          // For MinHash, both sketches are already sorted (smallest hashes first)
          // Count intersection of the two sketches
          size_t intersection_count = 0;
          size_t i = 0, j = 0;
          while (i < q_sketch.size() && j < t_sketch.size()) {
            if (q_sketch[i] == t_sketch[j]) {
              intersection_count++;
              i++;
              j++;
            } else if (q_sketch[i] < t_sketch[j]) {
              i++;
            } else {
              j++;
            }
          }

          if (intersection_count == 0) continue;

          // For MinHash, Jaccard = intersection / sketch_size (when both sketches have same size)
          double jaccard = static_cast<double>(intersection_count) / std::min(q_sketch.size(), t_sketch.size());
          
          double mash_dist = j2md(jaccard, estimation_k);
          double ani = 1.0 - mash_dist;
          all_anis.push_back(ani);
          
          comparison_count++;
          std::cerr << "[wfmash::auto-identity] Group " << q_pair.first << " vs " << t_pair.first 
                    << ": " << intersection_count << "/" << std::min(q_sketch.size(), t_sketch.size()) 
                    << " sketches overlap, Jaccard=" << std::fixed << std::setprecision(4) << jaccard 
                    << ", ANI=" << std::fixed << std::setprecision(2) << ani * 100 << "%" << std::endl;
        }
      }

      std::cerr << "[wfmash::auto-identity] Performed " << comparison_count << " group comparisons" << std::endl;

      if (all_anis.empty()) {
        std::cerr << "[wfmash::auto-identity] Warning: No k-mer overlap found between any query and target groups. Falling back to default identity." << std::endl;
        return skch::fixed::percentage_identity;
      }

      // Sort ANI values to compute percentiles
      std::sort(all_anis.begin(), all_anis.end());
      
      // Calculate the requested percentile
      size_t percentile_idx = (params.ani_percentile * all_anis.size()) / 100;
      if (percentile_idx >= all_anis.size()) {
        percentile_idx = all_anis.size() - 1;
      }
      
      // Apply adjustment
      double selected_ani = all_anis[percentile_idx];
      double adjusted_ani = selected_ani + (params.ani_adjustment / 100.0);
      
      // Clamp to valid range [0, 1]
      if (adjusted_ani < 0.0) adjusted_ani = 0.0;
      if (adjusted_ani > 1.0) adjusted_ani = 1.0;
      
      // Log the distribution of ANI values
      std::cerr << "[wfmash::auto-identity] ANI distribution: min=" 
                << std::fixed << std::setprecision(2) << all_anis.front() * 100 << "%, "
                << "25th percentile=" << all_anis[all_anis.size() / 4] * 100 << "%, "
                << "median=" << all_anis[all_anis.size() / 2] * 100 << "%, "
                << "75th percentile=" << all_anis[(3 * all_anis.size()) / 4] * 100 << "%, "
                << "max=" << all_anis.back() * 100 << "%" << std::endl;
      
      std::cerr << "[wfmash::auto-identity] Selected ani" << params.ani_percentile 
                << " (" << params.ani_percentile << "th percentile) = " 
                << std::fixed << std::setprecision(2) << selected_ani * 100 << "%";
      
      if (params.ani_adjustment != 0) {
        std::cerr << ", adjusted by " << std::showpos << params.ani_adjustment << std::noshowpos 
                  << "% to " << std::fixed << std::setprecision(2) << adjusted_ani * 100 << "%";
      }
      std::cerr << std::endl;
      
      return adjusted_ani;
    }
  }
}

#endif
