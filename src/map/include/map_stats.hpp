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
      // Use 4096 elements for accurate ANI estimation at lower identities
      const int estimation_sketch_size = 4096;

      std::cerr << "[wfmash::auto-identity] Starting identity estimation with k=" << estimation_k 
                << ", sketch_size=" << estimation_sketch_size << std::endl;

      // Maps to maintain MinHash sketches per group with per-group mutexes
      struct GroupSketch {
        StreamingMinHash sketch;
        std::mutex mutex;
        GroupSketch(size_t sketch_size) : sketch(sketch_size, 0) {}
      };
      std::map<int, std::unique_ptr<GroupSketch>> query_group_sketches;
      std::map<int, std::unique_ptr<GroupSketch>> target_group_sketches;
      std::mutex sketch_map_mutex;  // Only for creating new groups
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
      
      // Collect sequence names to process (NOT the sequences themselves)
      struct SeqToProcess {
        std::string name;
        std::string file;
        bool is_query;
        bool is_target;
      };
      
      // Use a map to track unique sequences and their roles
      std::unordered_map<std::string, SeqToProcess> unique_sequences;
      
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
            auto& seq_info = unique_sequences[seqname];
            seq_info.name = seqname;
            seq_info.file = file;
            seq_info.is_query = true;
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
            auto& seq_info = unique_sequences[seqname];
            seq_info.name = seqname;
            if (seq_info.file.empty()) seq_info.file = file;  // Use first file found
            seq_info.is_target = true;
          }
        }
        faidx_meta_destroy(meta);
      }
      
      // Convert map to vector for processing
      std::vector<SeqToProcess> sequences_to_process;
      sequences_to_process.reserve(unique_sequences.size());
      for (const auto& [name, info] : unique_sequences) {
        sequences_to_process.push_back(info);
      }
      
      // Count different types of sequences
      size_t query_only = 0, target_only = 0, both = 0;
      for (const auto& seq : sequences_to_process) {
        if (seq.is_query && seq.is_target) both++;
        else if (seq.is_query) query_only++;
        else if (seq.is_target) target_only++;
      }
      
      std::cerr << "[wfmash::auto-identity] Unique sequences to sketch: " << sequences_to_process.size() 
                << " (query-only: " << query_only
                << ", target-only: " << target_only
                << ", both: " << both << ")" << std::endl;
      
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
        
        // Function to get thread-local reader (per-file, per-thread)
        auto getThreadLocalReader = [&](faidx_meta_t* meta, const std::string& filename) -> faidx_reader_t* {
          thread_local std::map<std::string, faidx_reader_t*> readers;
          
          auto it = readers.find(filename);
          if (it == readers.end() || it->second == nullptr) {
            faidx_reader_t* reader = faidx_reader_create(meta);
            if (!reader) {
              throw std::runtime_error("Failed to create thread-local reader for " + filename);
            }
            readers[filename] = reader;
            return reader;
          }
          return it->second;
        };
        
        // Process each sequence in parallel with thread-local sketches
        for (const auto& seq_info : sequences_to_process) {
          sf.emplace([&, seq_info]() {
            try {
              auto meta_it = file_metas.find(seq_info.file);
              if (meta_it == file_metas.end()) return;
              
              faidx_meta_t* meta = meta_it->second;
              faidx_reader_t* reader = getThreadLocalReader(meta, seq_info.file);
              
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
              
              // Use lock-free streaming MinHash (no windowSize needed for pure MinHash)
              StreamingMinHash local_sketch(estimation_sketch_size, 0);
              
              // Process all k-mers in the sequence
              std::vector<char> kmerBuf(estimation_k);
              std::vector<char> seqRev(estimation_k);
              int ambig_kmer_count = 0;
              
              // Check initial k-mer for N's
              for (int j = 0; j < estimation_k && j < seq_len; j++) {
                char c = std::toupper(seq_data[j]);
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                  ambig_kmer_count = estimation_k;
                  break;
                }
              }
              
              // Process all k-mers
              for (offset_t i = 0; i <= (offset_t)seq_len - estimation_k; i++) {
                // Check for N at the end of current k-mer
                char endChar = std::toupper(seq_data[i + estimation_k - 1]);
                if (endChar != 'A' && endChar != 'C' && endChar != 'G' && endChar != 'T') {
                  ambig_kmer_count = estimation_k;
                }
                
                if (ambig_kmer_count == 0) {
                  // Copy and uppercase k-mer
                  for (int j = 0; j < estimation_k; j++) {
                    kmerBuf[j] = std::toupper(seq_data[i + j]);
                  }
                  
                  // Hash forward k-mer
                  hash_t hashFwd = CommonFunc::getHash(kmerBuf.data(), estimation_k);
                  hash_t hashBwd = std::numeric_limits<hash_t>::max();
                  
                  // Hash reverse complement for DNA
                  CommonFunc::reverseComplement(kmerBuf.data(), seqRev.data(), estimation_k);
                  hashBwd = CommonFunc::getHash(seqRev.data(), estimation_k);
                  
                  // Use canonical k-mer (minimum of forward and reverse)
                  if (hashFwd != hashBwd) {
                    hash_t canonicalHash = std::min(hashFwd, hashBwd);
                    local_sketch.add_unsafe(canonicalHash);  // No lock needed for local sketch
                  }
                }
                
                if (ambig_kmer_count > 0) {
                  ambig_kmer_count--;
                }
              }
              
              free(seq_data);
              
              // Get the final sketch
              auto sketch = local_sketch.getSketch();
              
              // Merge into group sketches with per-group locking
              if (seq_info.is_query) {
                // Ensure group exists
                {
                  std::lock_guard<std::mutex> lock(sketch_map_mutex);
                  if (query_group_sketches.find(groupId) == query_group_sketches.end()) {
                    query_group_sketches[groupId] = std::make_unique<GroupSketch>(estimation_sketch_size);
                  }
                }
                
                // Lock only this group's sketch
                std::lock_guard<std::mutex> lock(query_group_sketches[groupId]->mutex);
                for (hash_t h : sketch) {
                  query_group_sketches[groupId]->sketch.add_unsafe(h);
                }
                query_seq_count++;
              }
              
              if (seq_info.is_target) {
                // Ensure group exists
                {
                  std::lock_guard<std::mutex> lock(sketch_map_mutex);
                  if (target_group_sketches.find(groupId) == target_group_sketches.end()) {
                    target_group_sketches[groupId] = std::make_unique<GroupSketch>(estimation_sketch_size);
                  }
                }
                
                // Lock only this group's sketch
                std::lock_guard<std::mutex> lock(target_group_sketches[groupId]->mutex);
                for (hash_t h : sketch) {
                  target_group_sketches[groupId]->sketch.add_unsafe(h);
                }
                target_seq_count++;
              }
              
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

      if (query_seq_count == 0 || target_seq_count == 0) {
        std::cerr << "[wfmash::auto-identity] Warning: No sequences found in FAI files. " 
                  << "Make sure .fai index files exist (run 'samtools faidx' on your FASTA files)" << std::endl;
        return skch::fixed::percentage_identity;
      }

      // Convert sketches to vectors for comparison
      std::cerr << "[wfmash::auto-identity] Finalizing MinHash sketches..." << std::endl;
      
      std::map<int, std::vector<hash_t>> query_group_sketch_vectors;
      std::map<int, std::vector<hash_t>> target_group_sketch_vectors;
      
      for (auto& [groupId, groupSketch] : query_group_sketches) {
        query_group_sketch_vectors[groupId] = groupSketch->sketch.getSketch();
      }
      
      for (auto& [groupId, groupSketch] : target_group_sketches) {
        target_group_sketch_vectors[groupId] = groupSketch->sketch.getSketch();
      }
      
      std::cerr << "[wfmash::auto-identity] Processed " << query_seq_count 
                << " query sequences into " << query_group_sketch_vectors.size() << " groups" << std::endl;
      std::cerr << "[wfmash::auto-identity] Processed " << target_seq_count 
                << " target sequences into " << target_group_sketch_vectors.size() << " groups" << std::endl;

      // Compute total length for each group
      std::unordered_map<int, uint64_t> group_lengths;
      for (const auto& [name, id] : idManager.getSequenceNameToIdMap()) {
        int groupId = idManager.getRefGroup(id);
        uint64_t seqLen = idManager.getSequenceLength(id);
        group_lengths[groupId] += seqLen;
      }

      // Store ANI values with their weights (product of group lengths)
      std::vector<std::pair<double, uint64_t>> weighted_anis; // (ANI, weight)
      int comparison_count = 0;

      for (const auto& q_pair : query_group_sketch_vectors) {
        for (const auto& t_pair : target_group_sketch_vectors) {
          // Skip self-comparisons (same group ID)
          if (q_pair.first == t_pair.first) {
            continue; // Skip comparing a group to itself
          }
          
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
          
          // Weight is the product of both group lengths
          uint64_t weight = group_lengths[q_pair.first] * group_lengths[t_pair.first];
          weighted_anis.push_back({ani, weight});
          
          comparison_count++;
          std::cerr << "[wfmash::auto-identity] Group " << q_pair.first << " vs " << t_pair.first 
                    << ": " << intersection_count << "/" << std::min(q_sketch.size(), t_sketch.size()) 
                    << " sketches overlap, Jaccard=" << std::fixed << std::setprecision(4) << jaccard 
                    << ", ANI=" << std::fixed << std::setprecision(2) << ani * 100 << "%"
                    << " (weight=" << weight << ")"
                    << " [" << idManager.getGroupPrefix(q_pair.first) << " vs " << idManager.getGroupPrefix(t_pair.first) << "]"
                    << std::endl;
        }
      }

      // Count skipped self-comparisons
      int self_comparisons = 0;
      for (const auto& q_pair : query_group_sketch_vectors) {
        if (target_group_sketch_vectors.count(q_pair.first) > 0) {
          self_comparisons++;
        }
      }
      
      std::cerr << "[wfmash::auto-identity] Performed " << comparison_count << " group comparisons";
      if (self_comparisons > 0) {
        std::cerr << " (skipped " << self_comparisons << " self-comparisons)";
      }
      std::cerr << std::endl;

      if (weighted_anis.empty()) {
        std::cerr << "[wfmash::auto-identity] Warning: No k-mer overlap found between any query and target groups. Falling back to default identity." << std::endl;
        return skch::fixed::percentage_identity;
      }

      // Sort by ANI value for weighted percentile calculation
      std::sort(weighted_anis.begin(), weighted_anis.end(), 
                [](const auto& a, const auto& b) { return a.first < b.first; });
      
      // Calculate total weight
      uint64_t total_weight = 0;
      for (const auto& [ani, weight] : weighted_anis) {
        total_weight += weight;
      }
      
      // Calculate weighted percentile
      uint64_t target_weight = (params.ani_percentile * total_weight) / 100;
      uint64_t cumulative_weight = 0;
      double selected_ani = weighted_anis.back().first; // default to max if not found
      
      for (const auto& [ani, weight] : weighted_anis) {
        cumulative_weight += weight;
        if (cumulative_weight >= target_weight) {
          selected_ani = ani;
          break;
        }
      }
      
      // Apply adjustment
      double adjusted_ani = selected_ani + (params.ani_adjustment / 100.0);
      
      // Clamp to valid range [0, 1]
      if (adjusted_ani < 0.0) adjusted_ani = 0.0;
      if (adjusted_ani > 1.0) adjusted_ani = 1.0;
      
      // Extract just ANI values for simple statistics
      std::vector<double> all_anis;
      for (const auto& [ani, weight] : weighted_anis) {
        all_anis.push_back(ani);
      }
      
      // Log the distribution of ANI values (unweighted for reference)
      std::cerr << "[wfmash::auto-identity] ANI distribution: min=" 
                << std::fixed << std::setprecision(2) << all_anis.front() * 100 << "%, "
                << "25th percentile=" << all_anis[all_anis.size() / 4] * 100 << "%, "
                << "median=" << all_anis[all_anis.size() / 2] * 100 << "%, "
                << "75th percentile=" << all_anis[(3 * all_anis.size()) / 4] * 100 << "%, "
                << "max=" << all_anis.back() * 100 << "%" << std::endl;
      
      std::cerr << "[wfmash::auto-identity] Selected ani" << params.ani_percentile 
                << " (" << params.ani_percentile << "th length-weighted percentile) = " 
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
