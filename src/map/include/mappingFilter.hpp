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
#include <limits>
#include <fstream>
#include <iostream>
#include "map/include/base_types.hpp"
#include "map/include/filter.hpp"
#include "map/include/sequenceIds.hpp"
#include "common/dset64.hpp"
#include "common/progress.hpp"

namespace skch
{
  /**
   * @struct MappingAuxData
   * @brief A sidecar data structure holding all transient state for filtering/merging.
   */
  struct MappingAuxData {
      offset_t splitMappingId{0};
      double   chainPairScore{std::numeric_limits<double>::max()};
      int64_t  chainPairId{std::numeric_limits<int64_t>::min()};
      seqno_t  querySeqId{-1};
      offset_t queryLen{0};
      
      // We also need a place to store mutable end coordinates for gap-filling
      uint32_t tempQueryEndPos;
      uint32_t tempRefEndPos;
  };

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
                                  const SequenceIdManager& idManager,
                                  offset_t queryLen) // Pass queryLen
    {
      readMappings.erase(
          std::remove_if(readMappings.begin(),
                       readMappings.end(),
                       [&](MappingResult &e){
                           bool is_boundary_mapping = 
                               e.queryStartPos < param.segLength ||
                               e.queryEndPos() > queryLen - param.segLength || // Use passed-in queryLen
                               e.refStartPos < param.segLength ||
                               e.refEndPos() > idManager.getSequenceLength(e.refSeqId) - param.segLength;
                           
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
                             int64_t q_l = (int64_t)e.queryEndPos() - (int64_t)e.queryStartPos;
                             int64_t r_l = (int64_t)e.refEndPos() - (int64_t)e.refStartPos;
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
              auto a_strand = a.strand();
              auto b_strand = b.strand();
              return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos, a_strand) 
                  < std::tie(b.queryStartPos, b.refSeqId, b.refStartPos, b_strand);
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

          if (prev.refSeqId != curr.refSeqId || prev.strand() != curr.strand()) continue;

          int query_gap = curr.queryStartPos - prev.queryEndPos();
          int ref_gap = curr.refStartPos - prev.refEndPos();

          if (query_gap > 0 && ref_gap > 0 && query_gap <= threshold && ref_gap <= threshold) {
              int query_mid = (prev.queryEndPos() + curr.queryStartPos) / 2;
              int ref_mid = (prev.refEndPos() + curr.refStartPos) / 2;

              // TODO: This function needs to be rewritten to work with compact mappings
              // The algorithm tries to modify end positions which are now calculated
              // from blockLength. Need to adjust blockLength and start positions instead.
              /*
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
              */
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
          // TODO: Can't modify calculated fields
          //fragment.queryEndPos = std::max(fragment.queryEndPos(), it->queryEndPos());
          //fragment.refEndPos = std::max(fragment.refEndPos(), it->refEndPos());
      }

      // TODO: Recalculate blockLength based on the merged fragment
      //fragment.blockLength = std::max(fragment.refEndPos() - fragment.refStartPos, 
      //                               fragment.queryEndPos() - fragment.queryStartPos);
      //fragment.approxMatches = std::round(fragment.getNucIdentity() * fragment.blockLength / 100.0);
      fragment.n_merged = std::distance(start, end);

      // Average the scaled identity values
      uint16_t avgIdentity = std::accumulate(start, end, 0,
          [](int sum, const MappingResult& e) { return sum + e.nucIdentity; }
      ) / fragment.n_merged;
      fragment.nucIdentity = avgIdentity;

      // Average the scaled complexity values
      uint8_t avgComplexity = std::accumulate(start, end, 0,
          [](int sum, const MappingResult& e) { return sum + e.kmerComplexity; }
      ) / fragment.n_merged;
      fragment.kmerComplexity = avgComplexity;

      std::for_each(std::next(start), end, [](MappingResult& e) { e.setDiscard(true); });
    }

    /**
     * @brief Merge mappings within specified range
     */
    template <typename VecIn>
    static VecIn mergeMappingsInRange(VecIn &readMappings,
                                     int max_dist,
                                     const Parameters& param,
                                     progress_meter::ProgressMeter& progress,
                                     seqno_t querySeqId,
                                     offset_t queryLen)
    {
        if (!param.split || readMappings.size() < 2) return readMappings;

        // 1. Create and Populate the "Sidecar" Auxiliary Data
        std::vector<MappingAuxData> aux_data(readMappings.size());
        for (size_t i = 0; i < readMappings.size(); ++i) {
            aux_data[i].splitMappingId = i; // Assign a unique, stable ID for chaining
            aux_data[i].querySeqId = querySeqId;
            aux_data[i].queryLen = queryLen;
        }

        // 2. Correct Lexicographic Sort (using indices for efficiency)
        std::vector<uint32_t> p(readMappings.size());
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
            [&](uint32_t i, uint32_t j) {
                const auto& a = readMappings[i];
                const auto& b = readMappings[j];
                auto a_strand = a.strand();
                auto b_strand = b.strand();
                return std::tie(a.refSeqId, a_strand, a.queryStartPos, a.refStartPos)
                     < std::tie(b.refSeqId, b_strand, b.queryStartPos, b.refStartPos);
            });

        readMappings = reorder(readMappings, p);
        aux_data = reorder(aux_data, p);
        
        // 3. Restore the Union-Find and Geometric Chaining Logic
        std::vector<dsets::DisjointSets::Aint> ufv(readMappings.size());
        auto disjoint_sets = dsets::DisjointSets(ufv.data(), ufv.size());

        size_t group_start_idx = 0;
        while (group_start_idx < readMappings.size()) {
            const auto& start_map = readMappings[group_start_idx];
            size_t group_end_idx = group_start_idx + 1;
            while (group_end_idx < readMappings.size() &&
                   readMappings[group_end_idx].refSeqId == start_map.refSeqId &&
                   readMappings[group_end_idx].strand() == start_map.strand()) {
                group_end_idx++;
            }

            for (size_t i = group_start_idx; i < group_end_idx; ++i) {
                if (aux_data[i].chainPairScore != std::numeric_limits<double>::max()) {
                    disjoint_sets.unite(aux_data[i].splitMappingId, aux_data[i].chainPairId);
                }
                
                double best_score = std::numeric_limits<double>::max();
                size_t best_j = group_end_idx;

                for (size_t j = i + 1; j < group_end_idx; ++j) {
                    if (readMappings[j].queryStartPos > readMappings[i].queryEndPos() + max_dist) break;

                    int64_t q_dist = readMappings[j].queryStartPos - readMappings[i].queryEndPos();
                    if (q_dist < 0) q_dist = 0; // Prevent negative distance from small overlaps
                    
                    int64_t r_dist = (readMappings[i].strand() == strnd::FWD) ? 
                                     (readMappings[j].refStartPos - readMappings[i].refEndPos()) : 
                                     (readMappings[i].refStartPos - readMappings[j].refEndPos());
                    
                    if (q_dist <= max_dist && r_dist >= -param.segLength/5 && r_dist <= max_dist) {
                        double dist_sq = (double)q_dist * q_dist + (double)r_dist * r_dist;
                        if (dist_sq < best_score && dist_sq < aux_data[j].chainPairScore) {
                            best_score = dist_sq;
                            best_j = j;
                        }
                    }
                }

                if (best_j != group_end_idx) {
                    aux_data[best_j].chainPairScore = best_score;
                    aux_data[best_j].chainPairId = aux_data[i].splitMappingId;
                }
            }
            group_start_idx = group_end_idx;
        }

        // 4. Group by Chain and Construct Merged Fragments
        for (size_t i = 0; i < readMappings.size(); ++i) {
            if (aux_data[i].chainPairScore != std::numeric_limits<double>::max()) {
                 disjoint_sets.unite(aux_data[i].splitMappingId, aux_data[i].chainPairId);
            }
        }
        for (size_t i = 0; i < readMappings.size(); ++i) {
            aux_data[i].splitMappingId = disjoint_sets.find(aux_data[i].splitMappingId);
        }
        
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
            [&](uint32_t i, uint32_t j) {
                return std::tie(aux_data[i].splitMappingId, readMappings[i].queryStartPos, readMappings[i].refStartPos)
                     < std::tie(aux_data[j].splitMappingId, readMappings[j].queryStartPos, readMappings[j].refStartPos);
            });
        
        readMappings = reorder(readMappings, p);
        aux_data = reorder(aux_data, p);

        MappingResultsVector_t maximallyMergedMappings;
        size_t i = 0;
        while (i < readMappings.size()) {
            size_t j = i;
            while (j + 1 < readMappings.size() && aux_data[j+1].splitMappingId == aux_data[i].splitMappingId) {
                ++j;
            }

            size_t frag_start = i;
            while(frag_start <= j) {
                size_t frag_end = frag_start;
                while(frag_end + 1 <= j && 
                      (readMappings[frag_end+1].queryEndPos() - readMappings[frag_start].queryStartPos < param.max_mapping_length)) {
                    frag_end++;
                }
                
                MappingResult merged = readMappings[frag_start];
                uint32_t q_start = readMappings[frag_start].queryStartPos;
                uint32_t q_end = readMappings[frag_end].queryEndPos();
                uint32_t r_start_fwd = readMappings[frag_start].refStartPos;
                uint32_t r_end_fwd = readMappings[frag_end].refEndPos();
                
                double total_id = 0, total_comp = 0;
                uint32_t total_conserved = 0;
                for(size_t k = frag_start; k <= frag_end; ++k) {
                    total_id += readMappings[k].getNucIdentity();
                    total_comp += readMappings[k].getKmerComplexity();
                    total_conserved += readMappings[k].conservedSketches;
                    if (merged.strand() == strnd::REV) {
                        r_start_fwd = std::min(r_start_fwd, readMappings[k].refStartPos);
                        r_end_fwd = std::max(r_end_fwd, (uint32_t)readMappings[k].refEndPos());
                    }
                }

                merged.queryStartPos = q_start;
                merged.refStartPos = (merged.strand() == strnd::FWD) ? r_start_fwd : readMappings[frag_end].refStartPos;
                merged.blockLength = std::max(q_end - q_start, r_end_fwd - r_start_fwd);
                
                merged.n_merged = frag_end - frag_start + 1;
                merged.setNucIdentity(total_id / merged.n_merged);
                merged.setKmerComplexity(total_comp / merged.n_merged);
                merged.conservedSketches = total_conserved;
                
                maximallyMergedMappings.push_back(merged);
                frag_start = frag_end + 1;
            }
            i = j + 1;
        }
        
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
          if (it->queryStartPos - prev->queryEndPos() > param.segLength / 5 ||
              it->refStartPos - prev->refEndPos() > param.segLength / 5) {
              is_cuttable[std::distance(begin, prev)] = false;
              is_cuttable[std::distance(begin, it)] = false;
          }
      }

      adjustConsecutiveMappings(begin, end, param.segLength);

      auto fragment_start = begin;
      offset_t accumulate_length = 0;

      for (auto it = begin; it != end; ++it) {
          accumulate_length += it->queryEndPos() - it->queryStartPos;
          
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
          chain_start_query = std::min(chain_start_query, static_cast<offset_t>(it->queryStartPos));
          chain_end_query = std::max(chain_end_query, static_cast<offset_t>(it->queryEndPos()));
          chain_start_ref = std::min(chain_start_ref, static_cast<offset_t>(it->refStartPos));
          chain_end_ref = std::max(chain_end_ref, static_cast<offset_t>(it->refEndPos()));
          accumulate_nuc_identity += it->getNucIdentity();
          accumulate_kmer_complexity += it->getKmerComplexity();
          total_conserved_sketches += it->conservedSketches;
          //total_sketch_size += it->sketchSize; // Not in compact struct
      }

      auto chain_nuc_identity = accumulate_nuc_identity / n_in_full_chain;
      auto chain_kmer_complexity = accumulate_kmer_complexity / n_in_full_chain;
      auto block_length = std::max(chain_end_query - chain_start_query, chain_end_ref - chain_start_ref);

      for (auto it = begin; it != end; ++it) {
          it->n_merged = n_in_full_chain;
          it->blockLength = block_length;
          // blockNucIdentity is a method that returns getNucIdentity()
          it->setKmerComplexity(chain_kmer_complexity / 100.0f); // Convert back to 0-1 range
          it->conservedSketches = total_conserved_sketches;
          //it->sketchSize = total_sketch_size; // Not in compact struct
          //it->approxMatches = std::round(chain_nuc_identity * block_length / 100.0); // Not in compact struct
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
        int64_t q_span = m.queryEndPos() - m.queryStartPos;
        int64_t r_span = m.refEndPos() - m.refStartPos;
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
            double weight = m.queryEndPos() - m.queryStartPos;
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
            u_end   = (m.queryEndPos() + m.refEndPos()) * invSqrt2;
            v1 = (m.refStartPos - m.queryStartPos) * invSqrt2;
            v2 = (m.refEndPos() - m.queryEndPos()) * invSqrt2;
        } else {
            u_start = (m.refStartPos - m.queryStartPos) * invSqrt2;
            u_end   = (m.refEndPos() - m.queryEndPos()) * invSqrt2;
            v1 = (m.queryStartPos + m.refStartPos) * invSqrt2;
            v2 = (m.queryEndPos() + m.refEndPos()) * invSqrt2;
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
                                 progress_meter::ProgressMeter& progress,
                                 seqno_t querySeqId,   // Pass in context
                                 offset_t queryLen)    // Pass in context
    {
        if (param.scaffold_gap == 0) return; // Simplified check

        // 1. Build scaffolds using the now-correct merge function
        MappingResultsVector_t scaffoldMappings = readMappings; // Use raw mappings to build scaffolds
        Parameters scaffoldParam = param;
        scaffoldParam.chain_gap = param.scaffold_gap;
        auto superChains = mergeMappingsInRange(scaffoldMappings, scaffoldParam.chain_gap, scaffoldParam, progress, querySeqId, queryLen);
        
        // 2. Filter scaffolds by length
        superChains.erase(
            std::remove_if(superChains.begin(), superChains.end(),
                [&](const MappingResult& m) { return m.blockLength < param.scaffold_min_length; }),
            superChains.end());

        // 3. Expand scaffold boundaries for tolerance (RESTORED LOGIC)
        for (auto& mapping : superChains) {
            int64_t expansion = param.scaffold_max_deviation;
            
            uint32_t newQueryStart = (mapping.queryStartPos > expansion) ? mapping.queryStartPos - expansion : 0;
            uint32_t newRefStart = (mapping.refStartPos > expansion) ? mapping.refStartPos - expansion : 0;
            
            uint32_t newQueryEnd = std::min<offset_t>(queryLen, mapping.queryEndPos() + expansion);
            uint32_t newRefEnd = std::min<offset_t>(idManager.getSequenceLength(mapping.refSeqId), mapping.refEndPos() + expansion);

            mapping.blockLength = std::max(newQueryEnd - newQueryStart, newRefEnd - newRefStart);
            mapping.queryStartPos = newQueryStart;
            mapping.refStartPos = newRefStart;
        }

        // Output scaffolds to file if requested
        if (!param.scaffold_output_file.empty()) {
            static std::ofstream scaffoldOut;
            static bool firstWrite = true;
            
            if (firstWrite) {
                scaffoldOut.open(param.scaffold_output_file);
                if (!scaffoldOut.is_open()) {
                    std::cerr << "[wfmash] WARNING: Unable to open scaffold output file: " 
                              << param.scaffold_output_file << std::endl;
                } else {
                    firstWrite = false;
                }
            }
            
            if (scaffoldOut.is_open()) {
                for (const auto& scaffold : superChains) {
                    // Output in PAF format
                    scaffoldOut << idManager.getSequenceName(querySeqId)
                               << "\t" << queryLen
                               << "\t" << scaffold.queryStartPos
                               << "\t" << scaffold.queryEndPos()
                               << "\t" << (scaffold.strand() == strnd::FWD ? "+" : "-")
                               << "\t" << idManager.getSequenceName(scaffold.refSeqId)
                               << "\t" << idManager.getSequenceLength(scaffold.refSeqId)
                               << "\t" << scaffold.refStartPos
                               << "\t" << scaffold.refEndPos()
                               << "\t" << scaffold.conservedSketches
                               << "\t" << scaffold.blockLength
                               << "\t" << 60  // fake mapQ
                               << "\t" << "tp:A:S"  // tag to indicate this is a scaffold
                               << "\t" << "id:f:" << scaffold.getNucIdentity()
                               << "\t" << "kc:f:" << scaffold.getKmerComplexity()
                               << "\n";
                }
                scaffoldOut.flush();
            }
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
            GroupKey key { querySeqId, static_cast<seqno_t>(m.refSeqId) }; // FIX: Use the passed-in querySeqId
            rawGroups[key].push_back(m);
        }
        for (const auto& m : superChains) {
            GroupKey key { querySeqId, static_cast<seqno_t>(m.refSeqId) }; // FIX: Use the passed-in querySeqId
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
                    u_end = (m.queryEndPos() + m.refEndPos()) * invSqrt2;
                    v1 = (m.refStartPos - m.queryStartPos) * invSqrt2;
                    v2 = (m.refEndPos() - m.queryEndPos()) * invSqrt2;
                } else {
                    u_start = (m.refStartPos - m.queryStartPos) * invSqrt2;
                    u_end = (m.refEndPos() - m.queryEndPos()) * invSqrt2;
                    v1 = (m.queryStartPos + m.refStartPos) * invSqrt2;
                    v2 = (m.queryEndPos() + m.refEndPos()) * invSqrt2;
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

    /**
     * @brief Helper function to reorder a vector based on indices
     */
    template<typename T>
    static std::vector<T> reorder(const std::vector<T>& in, const std::vector<uint32_t>& p) {
        std::vector<T> out(in.size());
        for (size_t i = 0; i < p.size(); ++i) {
            out[i] = in[p[i]];
        }
        return out;
    }
  };
}

#endif // MAPPING_FILTER_HPP