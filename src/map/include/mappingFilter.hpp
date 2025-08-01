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
#include <atomic>
#include <thread>
#include <queue>
#include <memory>
#include "map/include/base_types.hpp"
#include "map/include/filter.hpp"
#include "map/include/sequenceIds.hpp"
#include "common/dset64.hpp"
#include "common/progress.hpp"

namespace skch
{
  /**
   * @struct Point2D
   * @brief Simple 2D point for spatial indexing
   */
  struct Point2D {
    float x, y;
    uint32_t idx;
    
    Point2D(float x_, float y_, uint32_t idx_) : x(x_), y(y_), idx(idx_) {}
  };

  /**
   * @class SimpleKDTree
   * @brief Simple 2D KD-tree for nearest neighbor search
   */
  class SimpleKDTree {
  private:
    struct Node {
      Point2D point;
      std::unique_ptr<Node> left;
      std::unique_ptr<Node> right;
      
      Node(const Point2D& p) : point(p) {}
    };
    
    std::unique_ptr<Node> root;
    std::vector<Point2D> points;
    
    std::unique_ptr<Node> buildTree(std::vector<Point2D>& pts, int start, int end, int depth) {
      if (start >= end) return nullptr;
      
      int axis = depth % 2;
      int mid = (start + end) / 2;
      
      // Sort by current axis
      if (axis == 0) {
        std::nth_element(pts.begin() + start, pts.begin() + mid, pts.begin() + end,
                        [](const Point2D& a, const Point2D& b) { return a.x < b.x; });
      } else {
        std::nth_element(pts.begin() + start, pts.begin() + mid, pts.begin() + end,
                        [](const Point2D& a, const Point2D& b) { return a.y < b.y; });
      }
      
      auto node = std::make_unique<Node>(pts[mid]);
      node->left = buildTree(pts, start, mid, depth + 1);
      node->right = buildTree(pts, mid + 1, end, depth + 1);
      
      return node;
    }
    
    void findKNearest(const Node* node, const Point2D& target, int k, int depth,
                      std::priority_queue<std::pair<float, uint32_t>>& heap) const {
      if (!node) return;
      
      float dist = std::sqrt((node->point.x - target.x) * (node->point.x - target.x) +
                            (node->point.y - target.y) * (node->point.y - target.y));
      
      if (node->point.idx != target.idx) {  // Don't include self
        if (heap.size() < k) {
          heap.push({dist, node->point.idx});
        } else if (dist < heap.top().first) {
          heap.pop();
          heap.push({dist, node->point.idx});
        }
      }
      
      int axis = depth % 2;
      float diff = (axis == 0) ? (target.x - node->point.x) : (target.y - node->point.y);
      
      Node* first = (diff < 0) ? node->left.get() : node->right.get();
      Node* second = (diff < 0) ? node->right.get() : node->left.get();
      
      findKNearest(first, target, k, depth + 1, heap);
      
      if (heap.size() < k || std::abs(diff) < heap.top().first) {
        findKNearest(second, target, k, depth + 1, heap);
      }
    }
    
  public:
    void build(std::vector<Point2D> pts) {
      points = std::move(pts);
      root = buildTree(points, 0, points.size(), 0);
    }
    
    std::vector<std::pair<uint32_t, float>> findKNearest(const Point2D& target, int k) const {
      std::priority_queue<std::pair<float, uint32_t>> heap;
      findKNearest(root.get(), target, k, 0, heap);
      
      std::vector<std::pair<uint32_t, float>> result;
      while (!heap.empty()) {
        result.push_back({heap.top().second, heap.top().first});
        heap.pop();
      }
      std::reverse(result.begin(), result.end());
      return result;
    }
  };
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
                               e.queryStartPos < param.windowLength ||
                               e.queryEndPos() > queryLen - param.windowLength || // Use passed-in queryLen
                               e.refStartPos < param.windowLength ||
                               e.refEndPos() > idManager.getSequenceLength(e.refSeqId) - param.windowLength;
                           
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
     * @brief Merge mappings within specified range (with chain info)
     */
    template <typename VecIn>
    static MappingsWithChains mergeMappingsInRangeWithChains(VecIn &readMappings,
                                     int max_dist,
                                     const Parameters& param,
                                     progress_meter::ProgressMeter& progress,
                                     seqno_t querySeqId,
                                     offset_t queryLen)
    {
        MappingsWithChains result;
        
        if (!param.split || readMappings.size() < 2) {
            result.mappings = readMappings;
            // For unmerged mappings, each is its own chain
            result.chainInfo.resize(readMappings.size());
            for (size_t i = 0; i < readMappings.size(); ++i) {
                result.chainInfo[i] = {static_cast<uint32_t>(i), 1, 1};
            }
            return result;
        }

        // Continue with existing merge logic...
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
                    
                    if (q_dist <= max_dist && r_dist >= -param.windowLength/5 && r_dist <= max_dist) {
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
            // Update progress for union-find processing
            progress.increment(1);
        }
        
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
            [&](uint32_t i, uint32_t j) {
                return std::tie(aux_data[i].splitMappingId, readMappings[i].queryStartPos, readMappings[i].refStartPos)
                     < std::tie(aux_data[j].splitMappingId, readMappings[j].queryStartPos, readMappings[j].refStartPos);
            });
        
        readMappings = reorder(readMappings, p);
        aux_data = reorder(aux_data, p);

        // Track chain assignments
        std::map<uint32_t, uint32_t> chainIdMap; // splitMappingId -> sequential chainId
        uint32_t nextChainId = 0;

        size_t i = 0;
        while (i < readMappings.size()) {
            size_t j = i;
            while (j + 1 < readMappings.size() && aux_data[j+1].splitMappingId == aux_data[i].splitMappingId) {
                ++j;
            }

            // Assign sequential chain ID if not seen before
            uint32_t chainId = nextChainId;
            auto it = chainIdMap.find(aux_data[i].splitMappingId);
            if (it == chainIdMap.end()) {
                chainIdMap[aux_data[i].splitMappingId] = nextChainId++;
                chainId = chainIdMap[aux_data[i].splitMappingId];
            } else {
                chainId = it->second;
            }

            // Count total mappings in this chain
            uint16_t chainLen = j - i + 1;
            uint16_t chainPos = 1;

            size_t frag_start = i;
            while(frag_start <= j) {
                size_t frag_end = frag_start;
                while(frag_end + 1 <= j) {
                    // Check both query and reference spans to prevent huge mappings
                    offset_t query_span = readMappings[frag_end+1].queryEndPos() - readMappings[frag_start].queryStartPos;
                    offset_t ref_span = readMappings[frag_end+1].refEndPos() - readMappings[frag_start].refStartPos;
                    if (std::max(query_span, ref_span) >= param.max_mapping_length) {
                        break;
                    }
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
                
                result.mappings.push_back(merged);
                result.chainInfo.push_back({chainId, chainPos++, chainLen});
                
                frag_start = frag_end + 1;
            }
            i = j + 1;
        }
        
        return result;
    }

    /**
     * @brief Merge mappings within specified range (legacy interface)
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
                    
                    if (q_dist <= max_dist && r_dist >= -param.windowLength/5 && r_dist <= max_dist) {
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
            // Update progress for union-find processing
            progress.increment(1);
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
                while(frag_end + 1 <= j) {
                    // Check both query and reference spans to prevent huge mappings
                    offset_t query_span = readMappings[frag_end+1].queryEndPos() - readMappings[frag_start].queryStartPos;
                    offset_t ref_span = readMappings[frag_end+1].refEndPos() - readMappings[frag_start].refStartPos;
                    if (std::max(query_span, ref_span) >= param.max_mapping_length) {
                        break;
                    }
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
          if (it->queryStartPos - prev->queryEndPos() > param.windowLength / 5 ||
              it->refStartPos - prev->refEndPos() > param.windowLength / 5) {
              is_cuttable[std::distance(begin, prev)] = false;
              is_cuttable[std::distance(begin, it)] = false;
          }
      }

      adjustConsecutiveMappings(begin, end, param.windowLength);

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



    /**
     * @brief Filter mappings by scaffolds using Direct Gravitational Model
     * Keep mappings that are within a direct distance threshold from any anchor
     */
    static void filterByScaffolds(MappingResultsVector_t& readMappings,
                                 const MappingResultsVector_t& mergedMappings,
                                 const Parameters& param,
                                 const SequenceIdManager& idManager,
                                 progress_meter::ProgressMeter& progress,
                                 seqno_t querySeqId,
                                 offset_t queryLen,
                                 std::shared_ptr<progress_meter::ProgressMeter> scaffold_progress = nullptr,
                                 std::shared_ptr<std::atomic<size_t>> scaffold_total_work = nullptr,
                                 std::shared_ptr<std::atomic<size_t>> scaffold_completed_work = nullptr)
    {
        if (param.scaffold_gap <= 0) return;
        
        // Add timing for expensive scaffold operations
        auto scaffold_start = std::chrono::high_resolution_clock::now();
        size_t initial_mappings = readMappings.size();

        // Phase 1: Anchor Identification - Modified approach
        MappingResultsVector_t scaffoldMappings = readMappings;
        Parameters scaffoldParam = param;
        scaffoldParam.chain_gap = param.scaffold_gap;
        
        // Step 1: Keep a copy of original mappings before merging
        MappingResultsVector_t originalMappings = scaffoldMappings;
        
        // Step 2: Merge to identify chains (this is expensive)
        auto merge_start = std::chrono::high_resolution_clock::now();
        auto mergedChains = mergeMappingsInRange(scaffoldMappings, scaffoldParam.chain_gap, scaffoldParam, progress, querySeqId, queryLen);
        auto merge_end = std::chrono::high_resolution_clock::now();
        auto merge_duration = std::chrono::duration_cast<std::chrono::milliseconds>(merge_end - merge_start);
        
        if (merge_duration.count() > 500) {
            std::cerr << "[wfmash::scaffold] Merging " << initial_mappings 
                      << " mappings took " << std::fixed << std::setprecision(1) 
                      << merge_duration.count() / 1000.0 << "s" << std::endl;
        }
        
        // Step 3: Filter merged chains by length
        mergedChains.erase(
            std::remove_if(mergedChains.begin(), mergedChains.end(),
                [&](const MappingResult& m) { return m.blockLength < param.scaffold_min_length; }),
            mergedChains.end());
        
        // Step 4: Apply plane sweep filter to the merged chains to remove off-diagonal/weaker scaffold chains
        if (!mergedChains.empty() && (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE)) {
            MappingResultsVector_t filteredChains;
            filterByGroup(mergedChains, filteredChains, param.numMappingsForSegment - 1, 
                         false, idManager, param, progress);
            mergedChains = std::move(filteredChains);
        }

        // Step 5: Collect original mappings that fall within the bounds of our filtered, high-quality scaffold chains.
        // These will become the anchors.
        MappingResultsVector_t allAnchorMappings;
        for (const auto& chain : mergedChains) {
            for (const auto& orig : originalMappings) {
                if (orig.refSeqId == chain.refSeqId &&
                    orig.strand() == chain.strand() &&
                    orig.queryStartPos >= chain.queryStartPos &&
                    orig.queryEndPos() <= chain.queryEndPos() &&
                    orig.refStartPos >= chain.refStartPos &&
                    orig.refEndPos() <= chain.refEndPos()) {
                    allAnchorMappings.push_back(orig);
                }
            }
        }
        
        // Use allAnchorMappings as anchors for the rest of the function
        auto& anchors = allAnchorMappings;
        
        if (readMappings.empty()) return;
        if (anchors.empty()) {
            // No anchors, filter out everything
            readMappings.clear();
            return;
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
                // Output the merged chains, not individual anchors
                for (const auto& scaffold : mergedChains) {
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
                               << "\t" << 60
                               << "\t" << "tp:A:S"
                               << "\t" << "id:f:" << scaffold.getNucIdentity()
                               << "\t" << "kc:f:" << scaffold.getKmerComplexity()
                               << "\n";
                }
                scaffoldOut.flush();
            }
        }

        // Phase 2: Build spatial index of ANCHORS ONLY
        std::vector<Point2D> anchor_points;
        anchor_points.reserve(anchors.size());
        for (size_t i = 0; i < anchors.size(); ++i) {
            const auto& m = anchors[i];
            float x = m.queryStartPos + m.blockLength * 0.5f;
            float y = m.refStartPos + m.blockLength * 0.5f;
            anchor_points.emplace_back(x, y, i);
        }
        
        SimpleKDTree anchor_kdtree;
        anchor_kdtree.build(anchor_points);
        
        // Phase 3: Compute direct distance to nearest anchor for each mapping
        std::vector<float> dist_to_nearest_anchor(readMappings.size());
        
        // Parallel computation of distances
        int num_threads = param.threads > 0 ? param.threads : std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        std::atomic<size_t> work_index(0);
        
        auto compute_anchor_distance = [&]() {
            size_t i;
            while ((i = work_index.fetch_add(1)) < readMappings.size()) {
                const auto& m = readMappings[i];
                Point2D p(m.queryStartPos + m.blockLength * 0.5f, 
                         m.refStartPos + m.blockLength * 0.5f, -1);
                
                auto nn = anchor_kdtree.findKNearest(p, 1);
                if (nn.empty()) {
                    dist_to_nearest_anchor[i] = std::numeric_limits<float>::infinity();
                } else {
                    dist_to_nearest_anchor[i] = nn[0].second;
                }
            }
        };
        
        for (int t = 0; t < num_threads; ++t) {
            threads.emplace_back(compute_anchor_distance);
        }
        for (auto& t : threads) {
            t.join();
        }
        
        // Phase 4: Apply distance filter
        MappingResultsVector_t keepers;
        float max_dist = static_cast<float>(param.scaffold_max_deviation);
        
        for (size_t i = 0; i < readMappings.size(); ++i) {
            if (dist_to_nearest_anchor[i] <= max_dist) {
                keepers.push_back(readMappings[i]);
            }
        }
        
        readMappings = std::move(keepers);
        
        // Report total scaffold filtering time if significant
        auto scaffold_end = std::chrono::high_resolution_clock::now();
        auto scaffold_duration = std::chrono::duration_cast<std::chrono::milliseconds>(scaffold_end - scaffold_start);
        
        if (scaffold_duration.count() > 1000) {
            std::cerr << "[wfmash::scaffold] Total scaffold filtering: " 
                      << initial_mappings << " → " << readMappings.size() 
                      << " mappings in " << std::fixed << std::setprecision(1) 
                      << scaffold_duration.count() / 1000.0 << "s" << std::endl;
        }
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