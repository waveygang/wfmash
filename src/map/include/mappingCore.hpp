/**
 * @file    mappingCore.hpp
 * @brief   Core L1/L2 mapping algorithms
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef MAPPING_CORE_HPP
#define MAPPING_CORE_HPP

#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/slidingMap.hpp"
#include "map/include/MIIteratorL2.hpp"
#include "map/include/map_stats.hpp"

namespace skch
{
  //Type for Stage L1's predicted candidate location
  struct L1_candidateLocus_t
  {
    seqno_t seqId;                    //sequence id where read is mapped
    offset_t rangeStartPos;
    offset_t rangeEndPos;
    int intersectionSize;
  };

  //Type for Stage L2's predicted mapping coordinate
  struct L2_mapLocus_t
  {
    seqno_t seqId;                    //sequence id where read is mapped
    offset_t meanOptimalPos;          //Among multiple consecutive optimal positions, save the avg.
    offset_t optimalStart;            //optimal start mapping position
    offset_t optimalEnd;              //optimal end mapping position
    int sharedSketchSize;             //count of shared sketch elements
    strand_t strand;
  };

  //Comparator for L1 candidate locations
  static constexpr auto L1_locus_intersection_cmp = [](L1_candidateLocus_t& a, L1_candidateLocus_t& b)
  {
    return a.intersectionSize < b.intersectionSize;
  };

  /**
   * @brief Template class containing core mapping algorithms
   */
  template<typename SketchType, typename IdManagerType>
  class MappingCore {
  public:
    using MI_Type = typename SketchType::MI_Type;
    using MIIter_t = typename SketchType::MIIter_t;

    /**
     * @brief Compute seed hits for a query sequence
     */
    template <typename Q_Info>
    static void getSeedHits(Q_Info &Q, const Parameters& param)
    {
      Q.minmerTableQuery.reserve(param.sketchSize + 1);
      CommonFunc::sketchSequence(Q.minmerTableQuery, Q.seq, Q.len, param.kmerSize, 
                                param.alphabetSize, param.sketchSize, Q.seqId);
      if(Q.minmerTableQuery.size() == 0) {
        Q.sketchSize = 0;
        return;
      }

      const double max_hash_01 = (long double)(Q.minmerTableQuery.back().hash) / std::numeric_limits<hash_t>::max();
      Q.kmerComplexity = (double(Q.minmerTableQuery.size()) / max_hash_01) / 
                         ((Q.len - param.kmerSize + 1)*2);
      Q.sketchSize = Q.minmerTableQuery.size();
    }

    /**
     * @brief Find seed interval points using reference sketch
     */
    template <typename Q_Info, typename Vec>
    static void getSeedIntervalPoints(Q_Info &Q, Vec& intervalPoints, 
                                     const SketchType* refSketch,
                                     const IdManagerType& idManager,
                                     const Parameters& param)
    {
      if(Q.minmerTableQuery.size() == 0)
        return;

      // Priority queue for sorting interval points
      using IP_const_iterator = typename std::vector<IntervalPoint>::const_iterator;
      std::vector<boundPtr<IP_const_iterator>> pq;
      pq.reserve(Q.sketchSize);
      constexpr auto heap_cmp = [](const auto& a, const auto& b) {return b < a;};

      for(auto it = Q.minmerTableQuery.begin(); it != Q.minmerTableQuery.end(); it++)
      {
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
        bool skip_mapping = false;
        int queryGroup = idManager.getRefGroup(Q.seqId);
        int targetGroup = idManager.getRefGroup(ip_it->seqId);

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
    }

    /**
     * @brief Compute L1 candidate regions from interval points
     */
    template <typename Q_Info, typename IP_iter, typename Vec2>
    static void computeL1CandidateRegions(
        Q_Info &Q, 
        IP_iter ip_begin, 
        IP_iter ip_end, 
        int minimumHits,
        const Parameters& param,
        const std::vector<int>& sketchCutoffs,
        Vec2 &l1Mappings)
    {
      int overlapCount = 0;
      int strandCount = 0;
      int bestIntersectionSize = 0;
      std::vector<L1_candidateLocus_t> localOpts;

      int windowLen = std::max<offset_t>(0, Q.len - param.segLength);
      auto trailingIt = ip_begin;
      auto leadingIt = ip_begin;
      int clusterLen = param.segLength;

      std::unordered_map<hash_t, int> hash_to_freq;

      // First pass to find best intersection size
      if (param.stage1_topANI_filter) {
        while (leadingIt != ip_end)
        {
          while (trailingIt != ip_end 
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
          bestIntersectionSize = std::max(bestIntersectionSize, overlapCount);
        }

        // Adjust minimum hits based on best intersection
        if (bestIntersectionSize < minimumHits) {
          return;
        } else {
          minimumHits = std::max(
              sketchCutoffs[
                int(std::min(bestIntersectionSize, Q.sketchSize) 
                  / std::max<double>(1, param.sketchSize / skch::fixed::ss_table_max))
              ],
              minimumHits);
        }
      }

      // Second pass to find candidate regions
      hash_to_freq.clear();
      bestIntersectionSize = std::min(bestIntersectionSize, Q.sketchSize);

      bool in_candidate = false;
      L1_candidateLocus_t l1_out = {};
      trailingIt = ip_begin;
      leadingIt = ip_begin;
      overlapCount = 0;
      int prevOverlap = 0;
      int prevPrevOverlap = 0;
      SeqCoord prevPos;
      SeqCoord currentPos{leadingIt->seqId, leadingIt->pos};

      while (leadingIt != ip_end)
      {
        prevPrevOverlap = prevOverlap;
        prevOverlap = overlapCount;

        while (trailingIt != ip_end 
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
        if (prevOverlap >= minimumHits)
        {
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

      // Join proximal local opts
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
     * @brief Compute L2 mapped regions within an L1 candidate
     */
    template <typename Q_Info, typename Vec>
    static void computeL2MappedRegions(Q_Info &Q,
        L1_candidateLocus_t &candidateLocus,
        Vec &l2_vec_out,
        const SketchType* refSketch,
        const Parameters& param)
    {
      auto& minmerIndex = refSketch->minmerIndex;
      const MinmerInfo first_minmer = MinmerInfo {0, candidateLocus.rangeStartPos - param.segLength - 1, 0, candidateLocus.seqId, 0};
      auto firstOpenIt = std::lower_bound(minmerIndex.begin(), minmerIndex.end(), first_minmer); 

      std::vector<skch::MinmerInfo> slidingWindow;
      slidingWindow.reserve(Q.sketchSize);
      constexpr auto heap_cmp = [](const skch::MinmerInfo& l, const skch::MinmerInfo& r) {return l.wpos_end > r.wpos_end;};

      auto windowIt = firstOpenIt;
      int windowLen = std::max<offset_t>(0, Q.len - param.segLength);
      std::unordered_map<hash_t, int> hash_to_freq;
      
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

      // Process the window
      while (windowIt != minmerIndex.end() && windowIt->seqId == candidateLocus.seqId && windowIt->wpos <= candidateLocus.rangeEndPos + windowLen) 
      {
        int prev_strand_votes = slideMap.strand_votes;
        bool inserted = false;
        
        while (!slidingWindow.empty() && slidingWindow.front().wpos_end <= windowIt->wpos - windowLen) {
          if (windowLen > 0) 
          {
            hash_to_freq[slidingWindow.front().hash]--;
          }
          if (windowLen == 0 || hash_to_freq[slidingWindow.front().hash] == 0) {
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

        if (slideMap.sharedSketchElements > bestSketchSize)
        {
          l2_vec_out.clear();
          in_candidate = true;
          bestSketchSize = slideMap.sharedSketchElements;
          l2_out.sharedSketchSize = slideMap.sharedSketchElements;
          l2_out.optimalStart = windowIt->wpos - windowLen;
          l2_out.optimalEnd = windowIt->wpos - windowLen;
        }
        else if(slideMap.sharedSketchElements == bestSketchSize)
        {
          if (!in_candidate) {
            l2_out.sharedSketchSize = slideMap.sharedSketchElements;
            l2_out.optimalStart = windowIt->wpos - windowLen;
          }
          in_candidate = true;
          l2_out.optimalEnd = windowIt->wpos - windowLen;
        } else {
          if (in_candidate) {
            l2_out.meanOptimalPos = (l2_out.optimalStart + l2_out.optimalEnd) / 2;
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
        l2_out.meanOptimalPos = (l2_out.optimalStart + l2_out.optimalEnd) / 2;
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
  };
}

#endif // MAPPING_CORE_HPP