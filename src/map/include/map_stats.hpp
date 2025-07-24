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

//External includes
#include "common/murmur3.h"
#include "common/kseq.h"
#include "common/prettyprint.hpp"
#include "common/seqiter.hpp"

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
     * @brief Helper function to add k-mers from a sequence to a hash multiset
     * @param[in] seq The sequence string
     * @param[out] all_hashes Multiset to store all k-mer hashes
     * @param[in] k K-mer size
     */
    inline void add_sequence_kmers(const std::string& seq, std::multiset<hash_t>& all_hashes, int k) {
      if ((int)seq.length() < k) return;

      char* seq_data = const_cast<char*>(seq.c_str());
      for (offset_t i = 0; i <= (offset_t)seq.length() - k; ++i) {
        all_hashes.insert(CommonFunc::getHash(seq_data + i, k));
      }
    }

    /**
     * @brief Compute MinHash sketch from a multiset of all k-mer hashes
     * @param[in] all_hashes All k-mer hashes from the group
     * @param[out] sketch Vector to store the MinHash sketch
     * @param[in] sketch_size Number of minimizers to retain
     */
    inline void compute_group_sketch(const std::multiset<hash_t>& all_hashes, 
                                     std::vector<hash_t>& sketch, int sketch_size) {
      sketch.clear();
      sketch.reserve(sketch_size);
      
      // MinHash: take the smallest sketch_size hash values
      auto it = all_hashes.begin();
      for (int i = 0; i < sketch_size && it != all_hashes.end(); ++i, ++it) {
        sketch.push_back(*it);
      }
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
      // For whole genome comparisons, use larger sketch size for better accuracy
      // 10000 provides excellent accuracy even for divergent genomes
      const int estimation_sketch_size = 10000;

      std::cerr << "[wfmash::auto-identity] Starting identity estimation with k=" << estimation_k 
                << ", sketch_size=" << estimation_sketch_size << std::endl;

      // Maps to collect all k-mers per group
      std::map<int, std::multiset<hash_t>> query_group_kmers;
      std::map<int, std::multiset<hash_t>> target_group_kmers;
      std::mutex kmer_mutex;
      std::atomic<int> query_seq_count(0);
      std::atomic<int> target_seq_count(0);

      // Determine number of threads to use
      int num_threads = params.threads > 0 ? params.threads : std::thread::hardware_concurrency();
      if (num_threads == 0) num_threads = 1;
      
      std::cerr << "[wfmash::auto-identity] Using " << num_threads << " threads for k-mer extraction" << std::endl;

      // First, collect all sequences to process
      struct SeqInfo {
        std::string name;
        std::string seq;
        bool is_query;
      };
      std::vector<SeqInfo> all_sequences;
      
      // Collect query sequences
      for (const auto& file : params.querySequences) {
        std::unordered_set<std::string> keep_seq; // Empty set means keep all
        std::string keep_prefix = ""; // Empty prefix means no prefix filtering
        seqiter::for_each_seq_in_file(file, keep_seq, keep_prefix, [&](const std::string& name, const std::string& seq) {
          // Only process sequences that are already in the SequenceIdManager
          const auto& seqMap = idManager.getSequenceNameToIdMap();
          if (seqMap.find(name) != seqMap.end()) {
            all_sequences.push_back({name, seq, true});
          }
        });
      }

      // Collect target sequences
      for (const auto& file : params.refSequences) {
        std::unordered_set<std::string> keep_seq; // Empty set means keep all
        std::string keep_prefix = ""; // Empty prefix means no prefix filtering
        seqiter::for_each_seq_in_file(file, keep_seq, keep_prefix, [&](const std::string& name, const std::string& seq) {
          // Only process sequences that are already in the SequenceIdManager
          const auto& seqMap = idManager.getSequenceNameToIdMap();
          if (seqMap.find(name) != seqMap.end()) {
            all_sequences.push_back({name, seq, false});
          }
        });
      }

      // Process sequences in parallel
      std::vector<std::thread> threads;
      size_t chunk_size = (all_sequences.size() + num_threads - 1) / num_threads;
      
      for (int t = 0; t < num_threads; t++) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, all_sequences.size());
        
        if (start >= all_sequences.size()) break;
        
        threads.emplace_back([&, start, end]() {
          // Thread-local temporary k-mer collections
          std::map<int, std::multiset<hash_t>> local_query_kmers;
          std::map<int, std::multiset<hash_t>> local_target_kmers;
          
          for (size_t i = start; i < end; i++) {
            const auto& seq_info = all_sequences[i];
            try {
              const auto& seqMap = idManager.getSequenceNameToIdMap();
              auto it = seqMap.find(seq_info.name);
              if (it == seqMap.end()) continue;
              
              seqno_t seqId = it->second;
              int groupId = idManager.getRefGroup(seqId);
              
              if (seq_info.is_query) {
                add_sequence_kmers(seq_info.seq, local_query_kmers[groupId], estimation_k);
                query_seq_count++;
              } else {
                add_sequence_kmers(seq_info.seq, local_target_kmers[groupId], estimation_k);
                target_seq_count++;
              }
            } catch (const std::exception& e) {
              std::cerr << "[wfmash::auto-identity] Warning: Failed to process sequence '" 
                        << seq_info.name << "': " << e.what() << std::endl;
            }
          }
          
          // Merge local k-mers into global maps
          std::lock_guard<std::mutex> lock(kmer_mutex);
          for (const auto& [groupId, kmers] : local_query_kmers) {
            query_group_kmers[groupId].insert(kmers.begin(), kmers.end());
          }
          for (const auto& [groupId, kmers] : local_target_kmers) {
            target_group_kmers[groupId].insert(kmers.begin(), kmers.end());
          }
        });
      }
      
      // Wait for all threads to complete
      for (auto& t : threads) {
        t.join();
      }

      std::cerr << "[wfmash::auto-identity] Collected k-mers from " << query_seq_count 
                << " query sequences into " << query_group_kmers.size() << " groups" << std::endl;
      std::cerr << "[wfmash::auto-identity] Collected k-mers from " << target_seq_count 
                << " target sequences into " << target_group_kmers.size() << " groups" << std::endl;

      if (query_seq_count == 0 || target_seq_count == 0) {
        std::cerr << "[wfmash::auto-identity] Warning: No sequences found in FAI files. " 
                  << "Make sure .fai index files exist (run 'samtools faidx' on your FASTA files)" << std::endl;
        return skch::fixed::percentage_identity;
      }

      // Now compute MinHash sketches for each group
      std::cerr << "[wfmash::auto-identity] Computing MinHash sketches for each group..." << std::endl;
      std::map<int, std::vector<hash_t>> query_group_sketches;
      std::map<int, std::vector<hash_t>> target_group_sketches;
      
      for (const auto& [groupId, kmers] : query_group_kmers) {
        compute_group_sketch(kmers, query_group_sketches[groupId], estimation_sketch_size);
      }
      
      for (const auto& [groupId, kmers] : target_group_kmers) {
        compute_group_sketch(kmers, target_group_sketches[groupId], estimation_sketch_size);
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
