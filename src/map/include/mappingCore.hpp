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

      // Gather matching points, then radix sort
      const std::size_t ip_start = intervalPoints.size();
      const int queryGroup = idManager.getRefGroup(Q.seqId);
      const bool checkGroup = param.skip_self || param.skip_prefix;
      for (const auto& m : Q.minmerTableQuery)
      {
        const auto seedFind = refSketch->minmerPosLookupIndex.find(m.hash);
        if (seedFind == refSketch->minmerPosLookupIndex.end()) continue;
        for (const auto& ip : seedFind->second)
        {
          if (checkGroup && queryGroup == idManager.getRefGroup(ip.seqId)) continue;
          if (param.lower_triangular && Q.seqId <= ip.seqId) continue;
          intervalPoints.push_back(ip);
        }
      }
      radixSortIntervalPoints(intervalPoints, ip_start);
    }

    /**
     * @brief LSD radix sort of [start, end) by (seqId, pos, side), the
     *        IntervalPoint::operator< order. Ties (same key, other hash) are
     *        order-insensitive downstream. Falls back to std::sort for small
     *        ranges or unpackable keys.
     */
    template <typename Vec>
    static void radixSortIntervalPoints(Vec& ip, std::size_t start)
    {
      const std::size_t n = ip.size() - start;
      if (n < 128) {
        std::sort(ip.begin() + start, ip.end());
        return;
      }
      // Key: seqId bits 33..63, pos 1..32, OPEN bit 0
      thread_local std::vector<uint64_t> keys;
      keys.resize(n);
      for (std::size_t i = 0; i < n; ++i) {
        const IntervalPoint& p = ip[start + i];
        if (p.seqId < 0 || p.pos < 0 || (uint64_t)p.pos >= (UINT64_C(1) << 32)) {
          std::sort(ip.begin() + start, ip.end());
          return;
        }
        keys[i] = ((uint64_t)(uint32_t)p.seqId << 33) | ((uint64_t)p.pos << 1) | (p.side == side::OPEN ? 1u : 0u);
      }
      // Six 11-bit passes; constant digits are skipped
      constexpr int RB = 11, RN = 1 << RB, NP = (64 + RB - 1) / RB;
      constexpr uint64_t RM = RN - 1;
      thread_local std::vector<std::size_t> hist;
      hist.assign((std::size_t)NP * RN, 0);
      for (std::size_t i = 0; i < n; ++i)
        for (int d = 0; d < NP; ++d) hist[(std::size_t)d * RN + ((keys[i] >> (d * RB)) & RM)]++;
      thread_local std::vector<uint32_t> ordA, ordB;
      ordA.resize(n); ordB.resize(n);
      for (uint32_t i = 0; i < (uint32_t)n; ++i) ordA[i] = i;
      uint32_t* src = ordA.data();
      uint32_t* dst = ordB.data();
      for (int d = 0; d < NP; ++d) {
        std::size_t* h = hist.data() + (std::size_t)d * RN;
        if (h[(keys[src[0]] >> (d * RB)) & RM] == n) continue;
        std::size_t sum = 0;
        for (int c = 0; c < RN; ++c) { const std::size_t t = h[c]; h[c] = sum; sum += t; }
        for (std::size_t i = 0; i < n; ++i) {
          const uint32_t idx = src[i];
          dst[h[(keys[idx] >> (d * RB)) & RM]++] = idx;
        }
        std::swap(src, dst);
      }
      thread_local std::vector<typename Vec::value_type> tmp;
      tmp.resize(n);
      for (std::size_t i = 0; i < n; ++i) tmp[i] = ip[start + src[i]];
      std::copy(tmp.begin(), tmp.end(), ip.begin() + start);
    }

    // One pos-group of the fused L1 sweep
    struct SweepStep {
      seqno_t seqId;
      offset_t pos;
      int overlapAfter;
    };

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
      if (ip_begin == ip_end) return;

      int overlapCount = 0;
      int bestIntersectionSize = 0;
      std::vector<L1_candidateLocus_t> localOpts;

      int windowLen = std::max<offset_t>(0, Q.len - param.windowLength);
      auto trailingIt = ip_begin;
      auto leadingIt = ip_begin;
      int clusterLen = param.windowLength;

      std::unordered_map<hash_t, int> hash_to_freq;

      bool in_candidate = false;
      L1_candidateLocus_t l1_out = {};

      // Emission state machine, fed the previous group
      auto emit = [&](int prevOverlap, const SeqCoord& prevPos) {
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
      };

      if (param.stage1_topANI_filter) {
        // One counting sweep records every pos-group; emission replays it
        // after the minimumHits raise. A seqId run with too few opens to
        // reach minimumHits collapses to one zero-overlap step.
        thread_local std::vector<SweepStep> steps;
        steps.clear();
        auto runStart = ip_begin;
        bool cleanStart = true;   // previous group did not spill into this run
        while (runStart != ip_end)
        {
          const seqno_t runSeqId = runStart->seqId;
          const auto runEnd = std::partition_point(runStart, ip_end,
              [runSeqId](const auto& p) { return p.seqId <= runSeqId; });
          if (windowLen == 0 && cleanStart
              && (int)(std::distance(runStart, runEnd) / 2) < minimumHits
              && !(runEnd != ip_end && runEnd->pos == std::prev(runEnd)->pos))
          {
            steps.push_back(SweepStep{runSeqId, runStart->pos, 0});
            trailingIt = leadingIt = runStart = runEnd;
            continue;
          }
          while (leadingIt < runEnd)
          {
            if (windowLen == 0) {
              // A CLOSE at pos C drains at group C: running sum
              const seqno_t groupSeqId = leadingIt->seqId;
              const offset_t groupPos = leadingIt->pos;
              while (leadingIt != ip_end && leadingIt->pos == groupPos) {
                overlapCount += (leadingIt->side == side::OPEN) ? 1 : -1;
                leadingIt++;
              }
              trailingIt = leadingIt;
              bestIntersectionSize = std::max(bestIntersectionSize, overlapCount);
              steps.push_back(SweepStep{groupSeqId, groupPos, overlapCount});
              continue;
            }
            while (trailingIt != ip_end 
                && ((trailingIt->seqId == leadingIt->seqId && trailingIt->pos <= leadingIt->pos - windowLen)
                  || trailingIt->seqId < leadingIt->seqId))
            {
              if (trailingIt->side == side::CLOSE) {
                hash_to_freq[trailingIt->hash]--;
                if (hash_to_freq[trailingIt->hash] == 0) {
                  overlapCount--;
                }
              }
              trailingIt++;
            }
            const seqno_t groupSeqId = leadingIt->seqId;
            const offset_t groupPos = leadingIt->pos;
            while (leadingIt != ip_end && leadingIt->pos == groupPos) {
              if (leadingIt->side == side::OPEN) {
                if (hash_to_freq[leadingIt->hash] == 0) {
                  overlapCount++;
                }
                hash_to_freq[leadingIt->hash]++;
              }
              leadingIt++;
            }
            bestIntersectionSize = std::max(bestIntersectionSize, overlapCount);
            steps.push_back(SweepStep{groupSeqId, groupPos, overlapCount});
          }
          cleanStart = (leadingIt == runEnd);
          runStart = leadingIt;
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

        // Step t emits with step t-1's overlap and coordinate
        int prevOverlap = 0;
        SeqCoord prevPos{};
        for (const auto& st : steps) {
          emit(prevOverlap, prevPos);
          prevOverlap = st.overlapAfter;
          prevPos = SeqCoord{st.seqId, st.pos};
        }
      }
      else
      {
        overlapCount = 0;
        int prevOverlap = 0;
        SeqCoord prevPos;
        SeqCoord currentPos{leadingIt->seqId, leadingIt->pos};

        while (leadingIt != ip_end)
        {
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
          emit(prevOverlap, prevPos);
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
      const MinmerInfo first_minmer = MinmerInfo {0, candidateLocus.rangeStartPos - param.windowLength - 1, 0, candidateLocus.seqId, 0};
      auto firstOpenIt = std::lower_bound(minmerIndex.begin(), minmerIndex.end(), first_minmer); 

      std::vector<skch::MinmerInfo> slidingWindow;
      slidingWindow.reserve(Q.sketchSize);
      constexpr auto heap_cmp = [](const skch::MinmerInfo& l, const skch::MinmerInfo& r) {return l.wpos_end > r.wpos_end;};

      auto windowIt = firstOpenIt;
      int windowLen = std::max<offset_t>(0, Q.len - param.windowLength);
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
                || l2_vec_out.back().optimalEnd + param.windowLength < l2_out.optimalStart)
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
            || l2_vec_out.back().optimalEnd + param.windowLength < l2_out.optimalStart)
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