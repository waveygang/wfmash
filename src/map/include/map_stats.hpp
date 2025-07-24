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
     * @brief Helper function to compute minihash sketch for a single sequence
     * @param[in] seq The sequence string
     * @param[out] sketch_set Set to store the computed sketch
     * @param[in] k K-mer size
     * @param[in] sketch_size Number of minimizers to retain
     */
    inline void get_minihash_sketch(const std::string& seq, std::set<hash_t>& sketch_set, int k, int sketch_size) {
      if ((int)seq.length() < k) return;

      std::vector<hash_t> hashes;
      char* seq_data = const_cast<char*>(seq.c_str());

      for (offset_t i = 0; i <= (offset_t)seq.length() - k; ++i) {
        hashes.push_back(CommonFunc::getHash(seq_data + i, k));
      }

      // Use nth_element for efficiency if sketch_size is much smaller than total hashes
      if ((int)hashes.size() > sketch_size) {
        std::nth_element(hashes.begin(), hashes.begin() + sketch_size, hashes.end());
        hashes.resize(sketch_size);
      }
      
      sketch_set.insert(hashes.begin(), hashes.end());
    }

    /**
     * @brief Estimate optimal identity threshold based on sequence groups
     * @param[in] params Mapping parameters
     * @param[in] idManager Sequence ID manager with grouping information
     * @return Estimated identity threshold [0,1]
     */
    inline double estimate_identity_for_groups(const skch::Parameters& params, const skch::SequenceIdManager& idManager) {
      // Use k=21 to match Mash defaults for genome comparison
      // This provides good sensitivity while avoiding too many random matches
      const int estimation_k = 21;
      // Use larger sketch size for better accuracy with whole genomes
      // Mash uses 1000 by default, but for large/divergent genomes we need more
      // 10000 provides good accuracy even for divergent comparisons
      const int estimation_sketch_size = 10000;

      std::cerr << "[wfmash::auto-identity] Starting identity estimation with k=" << estimation_k 
                << ", sketch_size=" << estimation_sketch_size << std::endl;

      std::map<int, std::set<hash_t>> query_group_sketches;
      std::map<int, std::set<hash_t>> target_group_sketches;

      // Sketch query sequences and group them using the provided idManager
      int query_seq_count = 0;
      for (const auto& file : params.querySequences) {
        std::unordered_set<std::string> keep_seq; // Empty set means keep all
        std::string keep_prefix = ""; // Empty prefix means no prefix filtering
        seqiter::for_each_seq_in_file(file, keep_seq, keep_prefix, [&](const std::string& name, const std::string& seq) {
          seqno_t seqId = idManager.getSequenceId(name);
          int groupId = idManager.getRefGroup(seqId);
          get_minihash_sketch(seq, query_group_sketches[groupId], estimation_k, estimation_sketch_size);
          query_seq_count++;
        });
      }

      // Sketch target sequences
      int target_seq_count = 0;
      for (const auto& file : params.refSequences) {
        std::unordered_set<std::string> keep_seq; // Empty set means keep all
        std::string keep_prefix = ""; // Empty prefix means no prefix filtering
        seqiter::for_each_seq_in_file(file, keep_seq, keep_prefix, [&](const std::string& name, const std::string& seq) {
          seqno_t seqId = idManager.getSequenceId(name);
          int groupId = idManager.getRefGroup(seqId);
          get_minihash_sketch(seq, target_group_sketches[groupId], estimation_k, estimation_sketch_size);
          target_seq_count++;
        });
      }

      std::cerr << "[wfmash::auto-identity] Processed " << query_seq_count << " query sequences into " 
                << query_group_sketches.size() << " groups" << std::endl;
      std::cerr << "[wfmash::auto-identity] Processed " << target_seq_count << " target sequences into " 
                << target_group_sketches.size() << " groups" << std::endl;

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

          std::vector<hash_t> intersection;
          std::set_intersection(q_sketch.begin(), q_sketch.end(),
                                t_sketch.begin(), t_sketch.end(),
                                std::back_inserter(intersection));

          if (intersection.empty()) continue;

          size_t union_size = q_sketch.size() + t_sketch.size() - intersection.size();
          double jaccard = static_cast<double>(intersection.size()) / union_size;
          
          double mash_dist = j2md(jaccard, estimation_k);
          double ani = 1.0 - mash_dist;
          all_anis.push_back(ani);
          
          comparison_count++;
          std::cerr << "[wfmash::auto-identity] Group " << q_pair.first << " vs " << t_pair.first 
                    << ": " << intersection.size() << "/" << std::min(q_sketch.size(), t_sketch.size()) 
                    << " sketches overlap, Jaccard=" << std::fixed << std::setprecision(4) << jaccard 
                    << ", ANI=" << std::fixed << std::setprecision(2) << ani * 100 << "%" << std::endl;
        }
      }

      std::cerr << "[wfmash::auto-identity] Performed " << comparison_count << " group comparisons" << std::endl;

      if (all_anis.empty()) {
        std::cerr << "[wfmash::auto-identity] Warning: No k-mer overlap found between any query and target groups. Falling back to default identity." << std::endl;
        return skch::fixed::percentage_identity;
      }

      // Use the 25th percentile (a quartile) as a robust and conservative estimate
      std::sort(all_anis.begin(), all_anis.end());
      size_t percentile_idx = all_anis.size() / 4;
      
      // Log the distribution of ANI values
      std::cerr << "[wfmash::auto-identity] ANI distribution: min=" 
                << std::fixed << std::setprecision(2) << all_anis.front() * 100 << "%, "
                << "25th percentile=" << all_anis[percentile_idx] * 100 << "%, "
                << "median=" << all_anis[all_anis.size() / 2] * 100 << "%, "
                << "max=" << all_anis.back() * 100 << "%" << std::endl;
      
      std::cerr << "[wfmash::auto-identity] Selected 25th percentile as threshold: " 
                << std::fixed << std::setprecision(2) << all_anis[percentile_idx] * 100 << "%" << std::endl;
      
      return all_anis[percentile_idx];
    }
  }
}

#endif
