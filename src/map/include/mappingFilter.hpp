/**
 * @file    mappingFilter.hpp
 * @brief   Filtering and merging logic for mappings
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef MAPPING_FILTER_HPP
#define MAPPING_FILTER_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <numeric>
#include <set>
#include "map/include/base_types.hpp"
#include "map/include/filter.hpp"
#include "map/include/sequenceIds.hpp"
#include "common/dset64.hpp"
#include "common/progress.hpp"

namespace skch
{
  /**
   * @brief Class containing filtering and merging algorithms for mappings
   */
  class MappingFilterUtils {
  public:
    /**
     * @brief Filter mappings with fewer than target merged base mappings
     */
    static void filterWeakMappings(MappingResultsVector_t &readMappings, 
                                  int64_t min_count,
                                  const Parameters& param,
                                  const SequenceIdManager& idManager)
    {
      readMappings.erase(
          std::remove_if(readMappings.begin(),
                       readMappings.end(),
                       [&](MappingResult &e){
                           bool is_boundary_mapping = 
                               e.queryStartPos < param.segLength ||
                               e.queryEndPos > e.queryLen - param.segLength ||
                               e.refStartPos < param.segLength ||
                               e.refEndPos > idManager.getSequenceLength(e.refSeqId) - param.segLength;
                           
                           if (is_boundary_mapping) {
                               return e.blockLength < param.block_length / 2 || 
                                      e.n_merged < min_count / 2;
                           } else {
                               return e.blockLength < param.block_length || 
                                      e.n_merged < min_count;
                           }
                       }),
          readMappings.end());
    }

    /**
     * @brief Filter mappings whose identity and query/ref length don't agree
     */
    static void filterFalseHighIdentity(MappingResultsVector_t &readMappings,
                                       const Parameters& param)
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
     * @brief Filter mappings by hash value for sparsification
     */
    static void sparsifyMappings(MappingResultsVector_t &readMappings,
                                const Parameters& param)
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
     * @brief Filter mappings by group
     */
    static void filterByGroup(
        MappingResultsVector_t &unfilteredMappings,
        MappingResultsVector_t &filteredMappings,
        int n_mappings,
        bool filter_ref,
        const SequenceIdManager& idManager,
        const Parameters& param,
        progress_meter::ProgressMeter& progress)
    {
      filteredMappings.reserve(unfilteredMappings.size());

      std::sort(unfilteredMappings.begin(), unfilteredMappings.end(), 
          [](const auto& a, const auto& b) { 
              return std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos); 
          });
          
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
            subrange_end = std::find_if_not(subrange_begin, unfilteredMappings.end(), 
                [&currGroup, &idManager] (const auto& candidate) {
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
              
          std::sort(tmpMappings.begin(), tmpMappings.end(), 
              [](const auto& a, const auto& b) { 
                  return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos) < 
                         std::tie(b.queryStartPos, b.refSeqId, b.refStartPos); 
              });
              
          if (filter_ref)
          {
              skch::Filter::ref::filterMappings(tmpMappings, idManager, n_mappings, 
                                               param.dropRand, param.overlap_threshold);
          }
          else
          {
              skch::Filter::query::filterMappings(tmpMappings, n_mappings, 
                                                 param.dropRand, param.overlap_threshold, progress);
          }
          
          filteredMappings.insert(
              filteredMappings.end(), 
              std::make_move_iterator(tmpMappings.begin()), 
              std::make_move_iterator(tmpMappings.end()));
          tmpMappings.clear();
          subrange_begin = subrange_end;
        }
      }
      
      std::sort(
          filteredMappings.begin(), filteredMappings.end(),
          [](const MappingResult &a, const MappingResult &b) {
              return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos, a.strand) 
                  < std::tie(b.queryStartPos, b.refSeqId, b.refStartPos, b.strand);
          });
    }

    /**
     * @brief Adjust consecutive mappings to close small gaps
     */
    static void adjustConsecutiveMappings(std::vector<MappingResult>::iterator begin_maping,
                                         std::vector<MappingResult>::iterator end_mapping,
                                         const int threshold)
    {
      if (std::distance(begin_maping, end_mapping) < 2) return;

      for (auto it = std::next(begin_maping); it != end_mapping; ++it) {
          auto& prev = *std::prev(it);
          auto& curr = *it;

          if (prev.refSeqId != curr.refSeqId || prev.strand != curr.strand) continue;

          int query_gap = curr.queryStartPos - prev.queryEndPos;
          int ref_gap = curr.refStartPos - prev.refEndPos;

          if (query_gap > 0 && ref_gap > 0 && query_gap <= threshold && ref_gap <= threshold) {
              int query_mid = (prev.queryEndPos + curr.queryStartPos) / 2;
              int ref_mid = (prev.refEndPos + curr.refStartPos) / 2;

              prev.queryEndPos = query_mid;
              prev.refEndPos = ref_mid;
              curr.queryStartPos = query_mid;
              curr.refStartPos = ref_mid;

              prev.blockLength = std::max(prev.refEndPos - prev.refStartPos, 
                                          prev.queryEndPos - prev.queryStartPos);
              curr.blockLength = std::max(curr.refEndPos - curr.refStartPos, 
                                          curr.queryEndPos - curr.queryStartPos);

              prev.approxMatches = std::round(prev.nucIdentity * prev.blockLength / 100.0);
              curr.approxMatches = std::round(curr.nucIdentity * curr.blockLength / 100.0);
          }
      }
    }

    /**
     * @brief Process a fragment of mappings
     */
    static void processMappingFragment(std::vector<MappingResult>::iterator start, 
                                      std::vector<MappingResult>::iterator end)
    {
      auto& fragment = *start;

      for (auto it = start; it != end; ++it) {
          fragment.queryStartPos = std::min(fragment.queryStartPos, it->queryStartPos);
          fragment.refStartPos = std::min(fragment.refStartPos, it->refStartPos);
          fragment.queryEndPos = std::max(fragment.queryEndPos, it->queryEndPos);
          fragment.refEndPos = std::max(fragment.refEndPos, it->refEndPos);
      }

      fragment.blockLength = std::max(fragment.refEndPos - fragment.refStartPos, 
                                     fragment.queryEndPos - fragment.queryStartPos);
      fragment.approxMatches = std::round(fragment.nucIdentity * fragment.blockLength / 100.0);
      fragment.n_merged = std::distance(start, end);

      fragment.nucIdentity = std::accumulate(start, end, 0.0,
          [](double sum, const MappingResult& e) { return sum + e.nucIdentity; }
      ) / fragment.n_merged;

      fragment.kmerComplexity = std::accumulate(start, end, 0.0,
          [](double sum, const MappingResult& e) { return sum + e.kmerComplexity; }
      ) / fragment.n_merged;

      std::for_each(std::next(start), end, [](MappingResult& e) { e.discard = 1; });
    }

    /**
     * @brief Merge mappings within specified range
     */
    template <typename VecIn>
    static VecIn mergeMappingsInRange(VecIn &readMappings,
                                     int max_dist,
                                     const Parameters& param,
                                     progress_meter::ProgressMeter& progress)
    {
      if (!param.split || readMappings.size() < 2) return readMappings;

      // Sort and assign unique IDs
      std::sort(readMappings.begin(), readMappings.end(),
          [](const MappingResult &a, const MappingResult &b) {
              return std::tie(a.refSeqId, a.strand, a.queryStartPos, a.refStartPos)
                   < std::tie(b.refSeqId, b.strand, b.queryStartPos, b.refStartPos);
          });

      for (auto it = readMappings.begin(); it != readMappings.end(); ++it) {
          it->splitMappingId = std::distance(readMappings.begin(), it);
          it->discard = 0;
          it->chainPairScore = std::numeric_limits<double>::max();
          it->chainPairId = std::numeric_limits<int64_t>::min();
      }

      // Set up union-find for merging
      std::vector<dsets::DisjointSets::Aint> ufv(readMappings.size());
      auto disjoint_sets = dsets::DisjointSets(ufv.data(), ufv.size());

      // Process each group
      auto group_begin = readMappings.begin();
      while (group_begin != readMappings.end()) {
          auto group_end = std::find_if_not(group_begin, readMappings.end(),
              [refSeqId = group_begin->refSeqId, strand = group_begin->strand](const MappingResult &m) {
                  return m.refSeqId == refSeqId && m.strand == strand;
              });

          // Find best pairs within group
          for (auto it = group_begin; it != group_end; ++it) {
              if (it->chainPairScore != std::numeric_limits<double>::max()) {
                  disjoint_sets.unite(it->splitMappingId, it->chainPairId);
              }
              double best_score = std::numeric_limits<double>::max();
              auto best_it2 = group_end;

              auto end_it2 = std::upper_bound(it + 1, group_end, it->queryEndPos + max_dist,
                  [](offset_t val, const MappingResult &m) {
                      return val < m.queryStartPos;
                  });

              for (auto it2 = it + 1; it2 != end_it2; ++it2) {
                  if (it2->queryStartPos == it->queryStartPos) continue;
                  int64_t query_dist = it2->queryStartPos - it->queryEndPos;
                  int64_t ref_dist = (it->strand == strnd::FWD) ? 
                                     it2->refStartPos - it->refEndPos : 
                                     it->refStartPos - it2->refEndPos;
                                     
                  if (query_dist <= max_dist && ref_dist >= -param.segLength/5 && ref_dist <= max_dist) {
                      double dist_sq = static_cast<double>(query_dist) * query_dist + 
                                       static_cast<double>(ref_dist) * ref_dist;
                      double max_dist_sq = static_cast<double>(max_dist) * max_dist;
                      if (dist_sq < max_dist_sq && dist_sq < best_score && dist_sq < it2->chainPairScore) {
                          best_it2 = it2;
                          best_score = dist_sq;
                      }
                  }
              }
              if (best_it2 != group_end) {
                  best_it2->chainPairScore = best_score;
                  best_it2->chainPairId = it->splitMappingId;
              }
          }
          group_begin = group_end;
      }

      // Assign merged IDs
      for (auto &mapping : readMappings) {
          mapping.splitMappingId = disjoint_sets.find(mapping.splitMappingId);
      }

      // Sort by merged ID
      std::sort(readMappings.begin(), readMappings.end(),
          [](const MappingResult &a, const MappingResult &b) {
              return std::tie(a.splitMappingId, a.queryStartPos, a.refStartPos)
                   < std::tie(b.splitMappingId, b.queryStartPos, b.refStartPos);
          });

      // Create maximally merged mappings
      MappingResultsVector_t maximallyMergedMappings;
      for (auto it = readMappings.begin(); it != readMappings.end();) {
          auto it_end = std::find_if(it, readMappings.end(), 
              [splitId = it->splitMappingId](const MappingResult &e) {
                  return e.splitMappingId != splitId;
              });

          // Process chain into fragments
          auto fragment_start = it;
          while (fragment_start != it_end) {
              MappingResult mergedMapping = *fragment_start;
              auto fragment_end = fragment_start;
              auto next = std::next(fragment_end);
             
              offset_t current_length = fragment_start->queryEndPos - fragment_start->queryStartPos;
              
              while (next != it_end) {
                  offset_t potential_length = next->queryEndPos - fragment_start->queryStartPos;
                  if (potential_length > param.max_mapping_length) break;
                  fragment_end = next;
                  current_length = potential_length;
                  next = std::next(next);
              }

              // Set merged mapping properties
              mergedMapping.queryStartPos = fragment_start->queryStartPos;
              mergedMapping.queryEndPos = fragment_end->queryEndPos;

              if (mergedMapping.strand == strnd::FWD) {
                  mergedMapping.refStartPos = fragment_start->refStartPos;
                  mergedMapping.refEndPos = fragment_end->refEndPos;
              } else {
                  mergedMapping.refStartPos = fragment_end->refStartPos;
                  mergedMapping.refEndPos = fragment_start->refEndPos;
              }
              
              mergedMapping.blockLength = std::max(
                  mergedMapping.refEndPos - mergedMapping.refStartPos,
                  mergedMapping.queryEndPos - mergedMapping.queryStartPos
              );

              // Calculate averages
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
              fragment_start = std::next(fragment_end);
          }

          // Process chain with splits
          processChainWithSplits(it, it_end, param);
          it = it_end;
      }

      // Remove discarded mappings
      readMappings.erase(
          std::remove_if(readMappings.begin(), readMappings.end(), 
                         [](const MappingResult& e) { return e.discard == 1; }),
          readMappings.end()
      );

      return maximallyMergedMappings;
    }

    /**
     * @brief Process a chain of mappings with splits
     */
    static void processChainWithSplits(std::vector<MappingResult>::iterator begin, 
                                      std::vector<MappingResult>::iterator end,
                                      const Parameters& param)
    {
      if (begin == end) return;

      std::vector<bool> is_cuttable(std::distance(begin, end), true);
      
      // Mark positions that are not cuttable
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
              processMappingFragment(fragment_start, std::next(it));
              fragment_start = std::next(it);
              accumulate_length = 0;
          }
      }

      if (fragment_start != end) {
          processMappingFragment(fragment_start, end);
      }

      computeChainStatistics(begin, end);
    }

    /**
     * @brief Compute chain statistics
     */
    static void computeChainStatistics(std::vector<MappingResult>::iterator begin, 
                                      std::vector<MappingResult>::iterator end)
    {
      if (begin == end) return;

      offset_t chain_start_query = std::numeric_limits<offset_t>::max();
      offset_t chain_end_query = std::numeric_limits<offset_t>::min();
      offset_t chain_start_ref = std::numeric_limits<offset_t>::max();
      offset_t chain_end_ref = std::numeric_limits<offset_t>::min();
      double accumulate_nuc_identity = 0.0;
      double accumulate_kmer_complexity = 0.0;
      int total_conserved_sketches = 0;
      int total_sketch_size = 0;
      int n_in_full_chain = std::distance(begin, end);

      for (auto it = begin; it != end; ++it) {
          chain_start_query = std::min(chain_start_query, it->queryStartPos);
          chain_end_query = std::max(chain_end_query, it->queryEndPos);
          chain_start_ref = std::min(chain_start_ref, it->refStartPos);
          chain_end_ref = std::max(chain_end_ref, it->refEndPos);
          accumulate_nuc_identity += it->nucIdentity;
          accumulate_kmer_complexity += it->kmerComplexity;
          total_conserved_sketches += it->conservedSketches;
          total_sketch_size += it->sketchSize;
      }

      auto chain_nuc_identity = accumulate_nuc_identity / n_in_full_chain;
      auto chain_kmer_complexity = accumulate_kmer_complexity / n_in_full_chain;
      auto block_length = std::max(chain_end_query - chain_start_query, chain_end_ref - chain_start_ref);

      for (auto it = begin; it != end; ++it) {
          it->n_merged = n_in_full_chain;
          it->blockLength = block_length;
          it->blockNucIdentity = chain_nuc_identity;
          it->kmerComplexity = chain_kmer_complexity;
          it->conservedSketches = total_conserved_sketches;
          it->sketchSize = total_sketch_size;
          it->approxMatches = std::round(chain_nuc_identity * block_length / 100.0);
      }
    }

    // Scaffold filtering structures and methods
    struct RotatedEnvelope {
        double u_start;
        double u_end;
        double v_min;
        double v_max;
        bool antidiagonal;
    };

    // Event types for sweep line algorithm
    enum EventType { START, END };
    enum MappingType { SCAFFOLD, RAW };

    struct Event {
        double u;
        EventType type;
        MappingType mappingType;
        double v_min, v_max;
        size_t id;
        
        bool operator<(const Event& other) const {
            if (u != other.u) return u < other.u;
            if (type != other.type) return type < other.type;
            return mappingType < other.mappingType;
        }
    };

    // Interval tree for v-coordinate ranges
    struct Interval {
        double low, high;
        size_t id;
        
        bool operator<(const Interval& other) const {
            return low < other.low;
        }
    };

    class IntervalTree {
        std::set<Interval> intervals;

    public:
        void insert(double low, double high, size_t id) {
            intervals.insert({low, high, id});
        }

        void remove(double low, double high, size_t id) {
            intervals.erase({low, high, id});
        }

        std::vector<size_t> findOverlapping(double low, double high) const {
            std::vector<size_t> result;
            auto it = intervals.upper_bound({low, 0, 0});
            if (it != intervals.begin()) --it;
            
            while (it != intervals.end() && it->low <= high) {
                if (!(it->high < low || it->low > high)) {
                    result.push_back(it->id);
                }
                ++it;
            }
            return result;
        }

        bool hasOverlap(double low, double high) const {
            auto it = intervals.upper_bound({low, 0, 0});
            if (it != intervals.begin()) --it;
            
            while (it != intervals.end() && it->low <= high) {
                if (!(it->high < low || it->low > high)) {
                    return true;
                }
                ++it;
            }
            return false;
        }
    };

    /**
     * @brief Compute orientation score for a mapping
     */
    static double computeOrientationScore(const MappingResult& m) {
        int64_t q_span = m.queryEndPos - m.queryStartPos;
        int64_t r_span = m.refEndPos - m.refStartPos;
        double diag_proj = (r_span + q_span) / std::sqrt(2.0);
        double anti_proj = (r_span - q_span) / std::sqrt(2.0);
        return std::abs(anti_proj / diag_proj);
    }

    /**
     * @brief Determine if group should use antidiagonal projection
     */
    static bool shouldUseAntidiagonal(const std::vector<MappingResult>& mappings) {
        double total_weight = 0.0;
        double weighted_score = 0.0;
        for (const auto& m : mappings) {
            double weight = m.queryEndPos - m.queryStartPos;
            total_weight += weight;
            weighted_score += weight * computeOrientationScore(m);
        }
        return (weighted_score / total_weight) > 1.0;
    }

    /**
     * @brief Compute rotated envelope for a mapping
     */
    static RotatedEnvelope computeRotatedEnvelope(const MappingResult& m, 
                                                  bool use_antidiagonal,
                                                  const Parameters& param) {
        const double invSqrt2 = 1.0 / std::sqrt(2.0);
        double u_start, u_end, v1, v2;
             
        if (!use_antidiagonal) {
            u_start = (m.queryStartPos + m.refStartPos) * invSqrt2;
            u_end   = (m.queryEndPos   + m.refEndPos)   * invSqrt2;
            v1 = (m.refStartPos - m.queryStartPos) * invSqrt2;
            v2 = (m.refEndPos   - m.queryEndPos)   * invSqrt2;
        } else {
            u_start = (m.refStartPos - m.queryStartPos) * invSqrt2;
            u_end   = (m.refEndPos   - m.queryEndPos)   * invSqrt2;
            v1 = (m.queryStartPos + m.refStartPos) * invSqrt2;
            v2 = (m.queryEndPos   + m.refEndPos)   * invSqrt2;
        }
             
        double u_min = std::min(u_start, u_end);
        double u_max = std::max(u_start, u_end);
        double v_min = std::min(v1, v2) - param.scaffold_max_deviation;
        double v_max = std::max(v1, v2) + param.scaffold_max_deviation;
        return { u_min, u_max, v_min, v_max, use_antidiagonal };
    }

    /**
     * @brief Filter mappings by scaffolds
     */
    static void filterByScaffolds(MappingResultsVector_t& readMappings,
                                 const MappingResultsVector_t& mergedMappings,
                                 const Parameters& param,
                                 const SequenceIdManager& idManager,
                                 progress_meter::ProgressMeter& progress)
    {
        if (param.scaffold_gap == 0 && param.scaffold_min_length == 0 && param.scaffold_max_deviation == 0) {
            return;
        }

        // Build scaffold mappings
        MappingResultsVector_t scaffoldMappings = mergedMappings;
        Parameters scaffoldParam = param;
        scaffoldParam.chain_gap *= 2;
        auto superChains = mergeMappingsInRange(scaffoldMappings, scaffoldParam.chain_gap, scaffoldParam, progress);
        
        // Filter and expand scaffolds
        filterWeakMappings(superChains, std::floor(param.scaffold_min_length / param.segLength), param, idManager);
        
        for (auto& mapping : superChains) {
            int64_t expansion = param.scaffold_max_deviation / 2;
            mapping.refStartPos = std::max<int64_t>(0, mapping.refStartPos - expansion);
            mapping.refEndPos = std::min<int64_t>(idManager.getSequenceLength(mapping.refSeqId),
                                                 mapping.refEndPos + expansion);
            mapping.queryStartPos = std::max<int64_t>(0, mapping.queryStartPos - expansion);
            mapping.queryEndPos = std::min<int64_t>(mapping.queryLen, mapping.queryEndPos + expansion);
            mapping.blockLength = std::max(mapping.refEndPos - mapping.refStartPos,
                                         mapping.queryEndPos - mapping.queryStartPos);
        }

        // Group mappings
        struct GroupKey {
            seqno_t querySeqId;
            seqno_t refSeqId;
            bool operator==(const GroupKey &other) const {
                return querySeqId == other.querySeqId && refSeqId == other.refSeqId;
            }
        };
        
        struct GroupKeyHash {
            std::size_t operator()(const GroupKey &k) const {
                auto h1 = std::hash<seqno_t>()(k.querySeqId);
                auto h2 = std::hash<seqno_t>()(k.refSeqId);
                return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
            }
        };

        std::unordered_map<GroupKey, std::vector<MappingResult>, GroupKeyHash> rawGroups;
        std::unordered_map<GroupKey, std::vector<MappingResult>, GroupKeyHash> scafGroups;
        
        for (const auto& m : readMappings) {
            GroupKey key { m.querySeqId, m.refSeqId };
            rawGroups[key].push_back(m);
        }
        for (const auto& m : superChains) {
            GroupKey key { m.querySeqId, m.refSeqId };
            scafGroups[key].push_back(m);
        }

        // Process each group with 2D sweep
        MappingResultsVector_t filteredMappings;
        for (const auto& kv : rawGroups) {
            const GroupKey& key = kv.first;
            const auto& groupRaw = kv.second;
            
            if (scafGroups.find(key) == scafGroups.end()) {
                continue;
            }
            const auto& groupScaf = scafGroups[key];

            bool use_antidiagonal = shouldUseAntidiagonal(groupScaf);

            // Generate events
            std::vector<Event> events;
            std::vector<bool> keep(groupRaw.size(), false);

            // Helper to compute rotated coordinates
            auto computeRotatedCoords = [use_antidiagonal](const MappingResult& m) {
                const double invSqrt2 = 1.0 / std::sqrt(2.0);
                double u_start, u_end, v1, v2;
                
                if (!use_antidiagonal) {
                    u_start = (m.queryStartPos + m.refStartPos) * invSqrt2;
                    u_end = (m.queryEndPos + m.refEndPos) * invSqrt2;
                    v1 = (m.refStartPos - m.queryStartPos) * invSqrt2;
                    v2 = (m.refEndPos - m.queryEndPos) * invSqrt2;
                } else {
                    u_start = (m.refStartPos - m.queryStartPos) * invSqrt2;
                    u_end = (m.refEndPos - m.queryEndPos) * invSqrt2;
                    v1 = (m.queryStartPos + m.refStartPos) * invSqrt2;
                    v2 = (m.queryEndPos + m.refEndPos) * invSqrt2;
                }
                
                return std::make_tuple(
                    std::min(u_start, u_end),
                    std::max(u_start, u_end),
                    std::min(v1, v2),
                    std::max(v1, v2)
                );
            };

            // Generate events for scaffolds
            for (size_t i = 0; i < groupScaf.size(); i++) {
                auto [u_min, u_max, v_min, v_max] = computeRotatedCoords(groupScaf[i]);
                v_min -= param.scaffold_max_deviation;
                v_max += param.scaffold_max_deviation;
                events.push_back(Event{u_min, START, SCAFFOLD, v_min, v_max, i});
                events.push_back(Event{u_max, END, SCAFFOLD, v_min, v_max, i});
            }

            // Generate events for raw mappings
            for (size_t i = 0; i < groupRaw.size(); i++) {
                auto [u_min, u_max, v_min, v_max] = computeRotatedCoords(groupRaw[i]);
                events.push_back(Event{u_min, START, RAW, v_min, v_max, i});
                events.push_back(Event{u_max, END, RAW, v_min, v_max, i});
            }

            std::sort(events.begin(), events.end());

            // Process events
            IntervalTree activeScaffolds;
            IntervalTree activeRaws;

            for (const auto& event : events) {
                if (event.type == START && event.mappingType == SCAFFOLD) {
                    activeScaffolds.insert(event.v_min, event.v_max, event.id);
                    auto overlapping = activeRaws.findOverlapping(event.v_min, event.v_max);
                    for (auto raw_id : overlapping) {
                        keep[raw_id] = true;
                    }
                } else if (event.type == END && event.mappingType == SCAFFOLD) {
                    activeScaffolds.remove(event.v_min, event.v_max, event.id);
                } else if (event.type == START && event.mappingType == RAW) {
                    if (!keep[event.id]) {
                        if (activeScaffolds.hasOverlap(event.v_min, event.v_max)) {
                            keep[event.id] = true;
                        } else {
                            activeRaws.insert(event.v_min, event.v_max, event.id);
                        }
                    }
                } else if (event.type == END && event.mappingType == RAW) {
                    if (!keep[event.id]) {
                        activeRaws.remove(event.v_min, event.v_max, event.id);
                    }
                }
            }

            // Collect kept mappings
            for (size_t i = 0; i < groupRaw.size(); i++) {
                if (keep[i]) {
                    filteredMappings.push_back(groupRaw[i]);
                }
            }
        }

        readMappings = std::move(filteredMappings);
    }
  };
}

#endif // MAPPING_FILTER_HPP