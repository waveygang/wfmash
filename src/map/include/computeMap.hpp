
/**
 * @file    computeMap.hpp
 * @brief   implements the sequence mapping logic
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
#include <cstring>
#include <iostream>
#include <unistd.h>
#include <queue>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/slidingMap.hpp"
#include "map/include/MIIteratorL2.hpp"
#include "map/include/ThreadPool.hpp"
#include "map/include/filter.hpp"

//External includes
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "map_stats.hpp"
#include "robin-hood-hashing/robin_hood.h"
// if we ever want to do the union-find chaining in parallel
//#include "common/dset64-gccAtomic.hpp"
// this is for single-threaded use, but is more portable
#include "common/dset64.hpp"
//#include "assert.hpp"
#include "gsl/gsl_randist.h"

namespace skch
{
  /**
   * @brief  LSD radix sort of interval points in [start, end) by (seqId, pos, side),
   *         reproducing IntervalPoint::operator< order. Replaces std::sort, which
   *         dominates mapping self-time at scale (all-vs-all makes the point count grow).
   *         (seqId, pos, side) are packed into an order-preserving uint64 key; ties among
   *         equal keys (same seqId/pos/side, differing hash) are order-insensitive downstream
   *         (verified byte-identical), so a stable radix is safe. Falls back to std::sort for
   *         small ranges or keys outside the packable range. thread_local scratch is reused to
   *         avoid per-call allocation (called once per query fragment).
   */
  template <typename Vec>
  inline void radixSortIntervalPoints(Vec& ip, std::size_t start)
  {
    const std::size_t n = ip.size() - start;
    if (n < 128) {                         // radix setup not worth it for tiny ranges
      std::sort(ip.begin() + start, ip.end());
      return;
    }

    thread_local std::vector<uint64_t> keys;
    keys.resize(n);
    // Pack: [ seqId : bits 33..63 ][ pos : bits 1..32 ][ sideOpen : bit 0 ]
    // side::CLOSE(-1)->0, side::OPEN(1)->1, so CLOSE sorts before OPEN (matches operator<).
    for (std::size_t i = 0; i < n; ++i) {
      const IntervalPoint& p = ip[start + i];
      if (p.seqId < 0 || p.pos < 0 ||
          (uint64_t)p.seqId >= (UINT64_C(1) << 31) ||
          (uint64_t)p.pos   >= (UINT64_C(1) << 32)) {   // out of packable range: bail
        std::sort(ip.begin() + start, ip.end());
        return;
      }
      const uint64_t sideOpen = (p.side == side::OPEN) ? 1u : 0u;
      keys[i] = ((uint64_t)(uint32_t)p.seqId << 33) | ((uint64_t)p.pos << 1) | sideOpen;
    }

    // 11-bit radix digits: 6 passes over a 64-bit key instead of 8 byte-passes, so
    // ~40% fewer index-scatter passes (this sort is movement-bound). 2048-bucket
    // histograms fit in L2. constant-digit passes are skipped, so in practice only
    // the ~3 populated digits (small seqId/pos) are scattered.
    constexpr int RB = 11;
    constexpr int RN = 1 << RB;          // 2048 buckets
    constexpr uint64_t RM = RN - 1;
    constexpr int NP = (64 + RB - 1) / RB;   // 6 passes
    thread_local std::vector<std::size_t> histbuf;
    histbuf.assign((std::size_t)NP * RN, 0);
    for (std::size_t i = 0; i < n; ++i) {
      const uint64_t k = keys[i];
      for (int d = 0; d < NP; ++d) histbuf[(std::size_t)d * RN + ((k >> (d * RB)) & RM)]++;
    }

    thread_local std::vector<uint32_t> ordA, ordB;
    ordA.resize(n); ordB.resize(n);
    for (uint32_t i = 0; i < (uint32_t)n; ++i) ordA[i] = i;
    uint32_t* src = ordA.data();
    uint32_t* dst = ordB.data();

    for (int d = 0; d < NP; ++d) {
      std::size_t* h = histbuf.data() + (std::size_t)d * RN;
      if (h[(keys[src[0]] >> (d * RB)) & RM] == n) continue;   // constant digit: skip pass
      std::size_t sum = 0;
      for (int c = 0; c < RN; ++c) { std::size_t t = h[c]; h[c] = sum; sum += t; }
      for (std::size_t i = 0; i < n; ++i) {
        const uint32_t idx = src[i];
        const uint32_t c = (uint32_t)((keys[idx] >> (d * RB)) & RM);
        dst[h[c]++] = idx;
      }
      std::swap(src, dst);
    }

    // Gather the permutation into scratch, then write back.
    thread_local std::vector<typename Vec::value_type> tmp;
    tmp.resize(n);
    for (std::size_t i = 0; i < n; ++i) tmp[i] = ip[start + src[i]];
    std::copy(tmp.begin(), tmp.end(), ip.begin() + start);
  }

  // encodePackedIP / decodePackedIP / PackedIPCursor live in base_types.hpp
  // (winSketch.hpp uses them to build the CSR position-lookup arena).

  // In-place LSD radix sort of packed uint64 keys in [start,end): 11-bit digits,
  // constant-digit skip, std::sort fallback for tiny ranges. The keys themselves
  // are the payload -- no index array, no struct gather.
  inline void radixSortPackedKeys(std::vector<uint64_t>& keys, std::size_t start) {
    const std::size_t n = keys.size() - start;
    if (n < 128) { std::sort(keys.begin() + start, keys.end()); return; }
    thread_local std::vector<uint64_t> buf;
    buf.resize(n);
    uint64_t* a = keys.data() + start;
    uint64_t* b = buf.data();

    constexpr int RB = 11, RN = 1 << RB, NP = (64 + RB - 1) / RB;
    constexpr uint64_t RM = RN - 1;
    thread_local std::vector<std::size_t> histbuf;
    histbuf.assign((std::size_t)NP * RN, 0);
    for (std::size_t i = 0; i < n; ++i) {
      const uint64_t k = a[i];
      for (int d = 0; d < NP; ++d) histbuf[(std::size_t)d * RN + ((k >> (d * RB)) & RM)]++;
    }
    uint64_t* src = a;
    uint64_t* dst = b;
    for (int d = 0; d < NP; ++d) {
      std::size_t* h = histbuf.data() + (std::size_t)d * RN;
      if (h[(src[0] >> (d * RB)) & RM] == n) continue;   // constant digit: skip
      std::size_t sum = 0;
      for (int c = 0; c < RN; ++c) { std::size_t t = h[c]; h[c] = sum; sum += t; }
      for (std::size_t i = 0; i < n; ++i) {
        const uint64_t k = src[i];
        dst[h[(k >> (d * RB)) & RM]++] = k;
      }
      std::swap(src, dst);
    }
    if (src != a) std::copy(src, src + n, a);   // odd #passes: result is in buf
  }

  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages
   */
  class Map
  {
    public:

      //Type for Stage L1's predicted candidate location
      struct L1_candidateLocus_t
      {
        seqno_t seqId;                    //sequence id where read is mapped

        /* read could be mapped with its begin location
         * from [rangeStartPos, rangeEndPos]
         */
        offset_t rangeStartPos;
        offset_t rangeEndPos;
        int intersectionSize;
      };

      static constexpr auto L1_locus_intersection_cmp = [](L1_candidateLocus_t& a, L1_candidateLocus_t& b)
      {
        return a.intersectionSize < b.intersectionSize;
      };

      //Type for Stage L2's predicted mapping coordinate within each L1 candidate
      struct L2_mapLocus_t
      {
        seqno_t seqId;                    //sequence id where read is mapped
        offset_t meanOptimalPos;          //Among multiple consecutive optimal positions, save the avg.
        offset_t optimalStart;            //optimal start mapping position (begin iterator)
        offset_t optimalEnd;              //optimal end mapping position (end iterator)
        int sharedSketchSize;             //count of shared sketch elements
        strand_t strand;
      };

    private:

      //algorithm parameters
      const skch::Parameters &param;

      //reference sketch
      const skch::Sketch &refSketch;

      //Container type for saving read sketches during L1 and L2 both
      typedef Sketch::MI_Type MinVec_Type;

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

      //Vector for obtaining group from refId
      //if refIdGroup[i] == refIdGroup[j], then sequence i and j have the same prefix;
      std::vector<int> refIdGroup;

      // True if every reference position fits in 32 bits and seqId in 31 bits, so
      // interval points can use the compact packed-uint64 representation (IP-1).
      bool packed_ip_ok = false;

      // Reference sequence name -> seqId, so the per-interval-point skip_self test can be
      // an integer compare (queryRefId != ip.seqId) instead of a std::string comparison.
      ankerl::unordered_dense::map<std::string, seqno_t> refNameToId;

      // Allowed (query, target) pairs from --pairs-file
      std::unordered_set<std::string> allowed_pairs;
      std::unordered_set<std::string> allowed_queries_from_pairs;

      // Memoized Stat::estimateMinimumHitsRelaxed for every possible Q.sketchSize;
      // all its other arguments are run constants, so this avoids the per-fragment
      // GSL binomial-CDF loops.
      std::vector<int> minHitsCache;

      // Per-seqId start offsets into refSketch.minmerIndex (globally sorted by
      // (seqId, wpos)), so the per-candidate lower_bound searches one contig
      // instead of the whole index.
      std::vector<size_t> minmerIndexSeqStart;

    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] refSketch   reference sketch
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      Map(const skch::Parameters &p, const skch::Sketch &refsketch,
          PostProcessResultsFn_t f = nullptr) :
        param(p),
        refSketch(refsketch),
        processMappingResults(f),
        sketchCutoffs(std::min<double>(p.sketchSize, skch::fixed::ss_table_max) + 1, 1),
        refIdGroup(refsketch.metadata.size())
    {
      if (p.stage1_topANI_filter) {
        this->setProbs();
      }
      if (p.skip_prefix)
      {
        this->setRefGroups();
      }
      if (!p.pairs_file.empty()) {
        this->loadPairsFile(p.pairs_file);
      }
      // The compact packed interval-point path is usable exactly when the Sketch
      // flattened its position lookup into the packed CSR arena.
      this->packed_ip_ok = refSketch.packed_ok;
      // Build reference name -> seqId once (for the integer skip_self test).
      this->refNameToId.reserve(refSketch.metadata.size());
      for (seqno_t i = 0; i < (seqno_t)refSketch.metadata.size(); ++i) {
        this->refNameToId[refSketch.metadata[i].name] = i;
      }
      // Memoize estimateMinimumHitsRelaxed over all possible sketch sizes.
      this->minHitsCache.resize(param.sketchSize + 1);
      for (int s = 0; s <= param.sketchSize; ++s) {
        this->minHitsCache[s] = Stat::estimateMinimumHitsRelaxed(s, param.kmerSize, param.percentageIdentity, skch::fixed::confidence_interval);
      }
      // Partition offsets of the (seqId, wpos)-sorted minmerIndex by seqId;
      // empty seqIds point at the next sequence's start.
      {
        const auto& mi = refSketch.minmerIndex;
        minmerIndexSeqStart.assign(refSketch.metadata.size() + 1, mi.size());
        for (size_t i = mi.size(); i-- > 0; ) {
          minmerIndexSeqStart[mi[i].seqId] = i;
        }
        for (size_t s = refSketch.metadata.size(); s-- > 0; ) {
          if (minmerIndexSeqStart[s] > minmerIndexSeqStart[s + 1])
            minmerIndexSeqStart[s] = minmerIndexSeqStart[s + 1];
        }
      }
      this->mapQuery();
    }

    private:

      // Sets the groups of reference contigs based on prefix
      void setRefGroups()
      {
        int group = 0;
        int start_idx = 0;
        int idx = 0;
        while (start_idx < this->refSketch.metadata.size())
        {
          const auto currPrefix = prefix(this->refSketch.metadata[start_idx].name, param.prefix_delim);
          idx = start_idx;
          while (idx < this->refSketch.metadata.size()
              && currPrefix == prefix(this->refSketch.metadata[idx].name, param.prefix_delim))
          {
            this->refIdGroup[idx++] = group;
          }
          group++;
          start_idx = idx;
        }
      }

      // Load allowed pairs from a TSV file (query<TAB>target per line)
      void loadPairsFile(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
          std::cerr << "[mashmap::skch::Map::loadPairsFile] ERROR: cannot open pairs file: " << filename << std::endl;
          exit(1);
        }
        std::string line;
        while (std::getline(in, line)) {
          // Skip empty lines and comments
          if (line.empty() || line[0] == '#') continue;
          // Each line should be query<TAB>target
          auto tab_pos = line.find('\t');
          if (tab_pos != std::string::npos) {
            std::string query_name = line.substr(0, tab_pos);
            std::string target_name = line.substr(tab_pos + 1);
            // Trim trailing whitespace/carriage return
            while (!target_name.empty() && (target_name.back() == '\r' || target_name.back() == ' ' || target_name.back() == '\t')) {
              target_name.pop_back();
            }
            allowed_pairs.insert(query_name + "\t" + target_name);
            allowed_queries_from_pairs.insert(query_name);
          }
        }
        std::cerr << "[mashmap::skch::Map::loadPairsFile] Loaded " << allowed_pairs.size() << " allowed pairs from " << filename << std::endl;
      }

      // Reference seqId whose name equals seqName, or -1 if none (for skip_self).
      seqno_t queryRefSeqId(const std::string& seqName) const
      {
        const auto it = refNameToId.find(seqName);
        return it != refNameToId.end() ? it->second : (seqno_t)-1;
      }

      // Gets the ref group of a query based on the prefix
      int getRefGroup(const std::string& seqName)
      {
        const auto queryPrefix = prefix(seqName, param.prefix_delim);
        for (int i = 0; i < this->refSketch.metadata.size(); i++)
        {
          const auto currPrefix = prefix(this->refSketch.metadata[i].name, param.prefix_delim);
          if (queryPrefix == currPrefix)
          {
            return this->refIdGroup[i];
          }
        }
        // Doesn't belong to any ref group
        return -1;
      }
      void setProbs() 
      {

        float deltaANI = param.ANIDiff;
        float min_p = 1 - param.ANIDiffConf;
        int ss = std::min<double>(param.sketchSize, skch::fixed::ss_table_max);

        // Cache hg pmf results
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
        
        // Return true iff Pr(ANI_i >= ANI_max - deltaANI) >= min_p
        const auto distDiff = [this, &sketchProbs, deltaANI, min_p, ss] (int cmax, int ci) {
          double prAboveCutoff = 0;
          for (double ymax = 0; ymax <= cmax; ymax++) {
            // Pr (Ymax = ymax)
            double pymax = sketchProbs[cmax][ymax];

            // yi_cutoff is minimum jaccard numerator required to be within deltaANI of ymax
            double yi_cutoff = deltaANI == 0 ? ymax : (std::floor(skch::Stat::md2j(
                skch::Stat::j2md(ymax / ss, param.kmerSize) + deltaANI, 
                param.kmerSize
            ) * ss));

            // Pr Y_i < yi_cutoff
            //std::cerr << "CMF " << yi_cutoff - 1 << " " << ss << " " << ss-ci << " " << ci << std::endl;
            double pi_acc = (yi_cutoff - 1) >= 0 ? gsl_cdf_hypergeometric_P(yi_cutoff-1, ss, ss-ci, ci) : 0;

            // Pr Y_i >= yi_cutoff
            pi_acc = 1-pi_acc;

            // Pr that mash score from cj leads to an ANI at least deltaJ less than the ANI from cmax
            prAboveCutoff += pymax * pi_acc;
            if (prAboveCutoff > min_p)
            {
              return true;
            }
          }
          return prAboveCutoff > min_p; 
        };

        // Helper vector for binary search
        std::vector<int> ss_range(ss+1);
        std::iota (ss_range.begin(), ss_range.end(), 0);

        for (auto cmax = 1; cmax <= ss; cmax++) 
        {
          // Binary search to find the lowest acceptable ci
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

          // For really high min_p values and some values of cmax, there are no values of
          // ci that satisfy the cutoff, so we just set to the max
          if (sketchCutoffs[cmax] == 0) {
            sketchCutoffs[cmax] = 1;
          }
        }
        //for (auto overlap = 1; overlap <= ss; overlap++) 
        //{
          //DEBUG_ASSERT(sketchCutoffs[overlap] <= overlap);
        //}
      }

      /**
       * @brief   parse over sequences in query file and map each on the reference
       */
      void mapQuery()
      {
        //Count of reads mapped by us
        //Some reads are dropped because of short length
        seqno_t totalReadsPickedForMapping = 0;
        seqno_t totalReadsMapped = 0;
        seqno_t seqCounter = 0;

        std::ofstream outstrm(param.outFileName);
        MappingResultsVector_t allReadMappings;  //Aggregate mapping results for the complete run

        //Create the thread pool
        ThreadPool<InputSeqProgContainer, MapModuleOutput> threadPool( [this](InputSeqProgContainer* e){return mapModule(e);}, param.threads);

		// allowed set of queries
		std::unordered_set<std::string> allowed_query_names;
		if (!param.query_list.empty()) {
			std::ifstream filter_list(param.query_list);
			std::string name;
			while (getline(filter_list, name)) {
				allowed_query_names.insert(name); 
			}
		}

		// Count the total number of sequences and sequence length
		uint64_t total_seqs = 0;
		uint64_t total_seq_length = 0;
		for (const auto& fileName : param.querySequences) {
			// Check if there is a .fai file
			std::string fai_name = fileName + ".fai";
			if ((access((fai_name).c_str(), F_OK) == 0)) {
				std::string line;
				std::ifstream in(fai_name.c_str());
                while (std::getline(in, line)) {
                    auto line_split = CommonFunc::split(line, '\t');
					auto seq_name = line_split[0];
					// Apply same filtering logic as for_each_seq_in_file_filtered:
					// skip only when a filter is active and doesn't match
					if (!param.query_prefix.empty()) {
						bool prefix_match = false;
						for (const auto& prefix : param.query_prefix) {
							if (seq_name.substr(0, prefix.size()) == prefix) {
								prefix_match = true;
								break;
							}
						}
						if (!prefix_match) continue;
					}
					if (!allowed_query_names.empty()
						&& allowed_query_names.find(seq_name) == allowed_query_names.end()) {
						continue;
					}
					total_seqs++;
					total_seq_length += std::stoul(line_split[1]);
				}
			}
#ifdef WFMASH_HAVE_AGC
			else if (agcidx::is_agc_file(fileName)) {
				// AGC: count sequences / sum lengths from metadata (GetCtgLen), no decompression
				agcidx::AgcIndex agc;
				if (!agc.open(fileName)) {
					std::cerr << "[mashmap::skch::Map::mapQuery] ERROR: could not open AGC archive " << fileName << std::endl;
					exit(1);
				}
				for (const auto& nl : agc.names_and_lengths()) {
					const std::string& seq_name = nl.first;
					if (!param.query_prefix.empty()) {
						bool prefix_match = false;
						for (const auto& prefix : param.query_prefix) {
							if (seq_name.substr(0, prefix.size()) == prefix) { prefix_match = true; break; }
						}
						if (!prefix_match) continue;
					}
					if (!allowed_query_names.empty()
						&& allowed_query_names.find(seq_name) == allowed_query_names.end()) continue;
					++total_seqs;
					total_seq_length += nl.second;
				}
			}
#endif
			else {
				// If .fai file doesn't exist, warn and use the for_each_seq_in_file_filtered function
				std::cerr << "[mashmap::skch::Map::mapQuery] WARNING, no .fai index found for " << fileName << ", reading the file to filter query sequences (slow)" << std::endl;
				seqiter::for_each_seq_in_file_filtered(
					fileName,
					param.query_prefix,
					allowed_query_names,
					[&](const std::string& seq_name, const std::string& seq) {
						++total_seqs;
						total_seq_length += seq.size();
					});
			}
		}
		
        progress_meter::ProgressMeter progress(total_seq_length, "[mashmap::skch::Map::mapQuery] mapped");

        for(const auto &fileName : param.querySequences)
        {

#ifdef DEBUG
            std::cerr << "[mashmap::skch::Map::mapQuery] mapping reads in " << fileName << std::endl;
#endif

			seqiter::for_each_seq_in_file_filtered(
				fileName,
				param.query_prefix,
				allowed_query_names,
                [&](const std::string& seq_name,
                    const std::string& seq) {
                    // todo: offset_t is an 32-bit integer, which could cause problems
                    offset_t len = seq.length();
					if (param.skip_self
						&& param.target_prefix != ""
						&& seq_name.substr(0, param.target_prefix.size()) == param.target_prefix) {
						// skip
					} else if (!allowed_queries_from_pairs.empty()
						&& allowed_queries_from_pairs.find(seq_name) == allowed_queries_from_pairs.end()) {
						// skip: query not in any allowed pair
					} else {
						if (param.filterMode == filter::ONETOONE)
							qmetadata.push_back( ContigInfo{seq_name, len} );
						//Is the read too short?
						if(len < param.kmerSize)
						{
//#ifdef DEBUG
							// TODO Should we somehow revert to < windowSize?
							std::cerr << std::endl
									  << "WARNING, skch::Map::mapQuery, read "
									  << seq_name << " of " << len << "bp "
									  << " is not long enough for mapping at segment length "
									  << param.segLength << std::endl;
//#endif
						}
						else
						{
							totalReadsPickedForMapping++;
							//Dispatch input to thread
							threadPool.runWhenThreadAvailable(new InputSeqProgContainer(seq, seq_name, seqCounter, progress));

							//Collect output if available
							while ( threadPool.outputAvailable() ) {
								mapModuleHandleOutput(threadPool.popOutputWhenAvailable(), allReadMappings, totalReadsMapped, outstrm, progress);
							}
						}
						//progress.increment(seq.size()/2);
						seqCounter++;
					}
                }); //Finish reading query input file

        }

        //Collect remaining output objects
        while ( threadPool.running() )
            mapModuleHandleOutput(threadPool.popOutputWhenAvailable(), allReadMappings, totalReadsMapped, outstrm, progress);

        //Filter over reference axis and report the mappings
        if (param.filterMode == filter::ONETOONE)
        {
          // how many secondary mappings to keep
          int n_mappings = param.numMappingsForSegment - 1;

          // Group sequences by query prefix, then pass to ref filter
          auto subrange_begin = allReadMappings.begin();
          auto subrange_end = allReadMappings.begin();
          MappingResultsVector_t tmpMappings;
          MappingResultsVector_t filteredMappings;

          // Precompute each query's group once (getRefGroup is O(#ref contigs))
          std::vector<int> queryGroup;
          if (param.skip_prefix)
          {
            queryGroup.resize(qmetadata.size());
            for (size_t i = 0; i < qmetadata.size(); i++)
              queryGroup[i] = this->getRefGroup(qmetadata[i].name);
          }

          while (subrange_end != allReadMappings.end())
          {
            if (param.skip_prefix)
            {
              int currGroup = queryGroup[subrange_begin->querySeqId];
              subrange_end = std::find_if_not(subrange_begin, allReadMappings.end(), [&queryGroup, currGroup] (const auto& allReadMappings_candidate) {
                  return currGroup == queryGroup[allReadMappings_candidate.querySeqId];
              });
            }
            else
            {
              subrange_end = allReadMappings.end();
            }
            tmpMappings.insert(
                tmpMappings.end(), 
                std::make_move_iterator(subrange_begin), 
                std::make_move_iterator(subrange_end));

            // tmpMappings now contains mappings from one group of query sequences to all reference groups
            // we now run filterByGroup, which filters based on the reference group.
            filterByGroup(tmpMappings, filteredMappings, n_mappings, true);
            tmpMappings.clear();
            subrange_begin = subrange_end;
          }
          allReadMappings = std::move(filteredMappings);

          //Re-sort mappings by input order of query sequences
          //This order may be needed for any post analysis of output
          std::sort(
              allReadMappings.begin(), allReadMappings.end(),
              [](const MappingResult &a, const MappingResult &b) {
                  return std::tie(a.querySeqId, a.queryStartPos, a.refSeqId, a.refStartPos) 
                    < std::tie(b.querySeqId, b.queryStartPos, b.refSeqId, b.refStartPos);
              });

          reportReadMappings(allReadMappings, "", outstrm);
        }

        progress.finish();

        std::cerr << "[mashmap::skch::Map::mapQuery] "
                  << "count of mapped reads = " << totalReadsMapped
                  << ", reads qualified for mapping = " << totalReadsPickedForMapping
                  << ", total input reads = " << seqCounter
                  << ", total input bp = " << total_seq_length << std::endl;

      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings with fewer than the target number of merged base mappings
       * @param[in]   input   mappings
       * @return              void
       */
      void filterWeakMappings(MappingResultsVector_t &readMappings, int64_t min_count)
      {
          readMappings.erase(
              std::remove_if(readMappings.begin(),
                             readMappings.end(),
                             [&](MappingResult &e){
                                 return e.queryLen > e.blockLength
                                     && e.n_merged < min_count;
                             }),
              readMappings.end());
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings whose identity and query/ref length don't agree
       * @param[in]   input   mappings
       * @return              void
       */
      void filterFalseHighIdentity(MappingResultsVector_t &readMappings)
      {
          readMappings.erase(
              std::remove_if(readMappings.begin(),
                             readMappings.end(),
                             [&](MappingResult &e){
                                 int64_t q_l = (int64_t)e.queryEndPos - (int64_t)e.queryStartPos;
                                 int64_t r_l = (int64_t)e.refEndPos + 1 - (int64_t)e.refStartPos;
                                 uint64_t delta = std::abs(r_l - q_l);
                                 float len_id_bound = (1.0 - (float)delta/(float)q_l);
                                 return len_id_bound < std::min(0.7, std::pow(param.percentageIdentity,3));
                             }),
              readMappings.end());
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings whose split ids aren't to be kept
       * @param[in]   input   mappings
       * @param[in]   input
       * @return              void
       */
      void filterFailedSubMappings(MappingResultsVector_t &readMappings,
                                   const robin_hood::unordered_set<offset_t>& kept_chains)
      {
          readMappings.erase(
              std::remove_if(readMappings.begin(),
                             readMappings.end(),
                             [&](MappingResult &e){
                                 return kept_chains.count(e.splitMappingId) == 0;
                             }),
              readMappings.end());
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings by hash value
       * @param[in]   input   mappings
       * @param[in]   input
       * @return              void
       */
      void sparsifyMappings(MappingResultsVector_t &readMappings)
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
       * @brief                    helper to main filtering function
       * @details                  filters mappings by group
       * @param[in]   input        unfiltered mappings
       * @param[in]   output       filtered mappings
       * @param[in]   n_mappings   num mappings per segment
       * @param[in]   filter_ref   use Filter::ref instead of Filter::query
       * @return                   void
       */
      void filterByGroup(
          MappingResultsVector_t &unfilteredMappings,
          MappingResultsVector_t &filteredMappings,
          int n_mappings,
          bool filter_ref)
      {
        filteredMappings.reserve(unfilteredMappings.size());

        std::sort(unfilteredMappings.begin(), unfilteredMappings.end(), [](const auto& a, const auto& b) 
            { return std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos); });
        auto subrange_begin = unfilteredMappings.begin();
        auto subrange_end = unfilteredMappings.begin();
        if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) 
        {
          std::vector<skch::MappingResult> tmpMappings;
          while (subrange_end != unfilteredMappings.end())
          {
            if (param.skip_prefix)
            {
              int currGroup = this->refIdGroup[subrange_begin->refSeqId];
              subrange_end = std::find_if_not(subrange_begin, unfilteredMappings.end(), [this, currGroup] (const auto& unfilteredMappings_candidate) {
                  return currGroup == this->refIdGroup[unfilteredMappings_candidate.refSeqId];
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
            std::sort(tmpMappings.begin(), tmpMappings.end(), [](const auto& a, const auto& b) 
                { return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos) < std::tie(b.queryStartPos, b.refSeqId, b.refStartPos); });
            if (filter_ref)
            {
              skch::Filter::ref::filterMappings(tmpMappings, this->refSketch, n_mappings);
            }
            else
            {
              skch::Filter::query::filterMappings(tmpMappings, n_mappings, param.dropRand);
            }
            filteredMappings.insert(
                filteredMappings.end(), 
                std::make_move_iterator(tmpMappings.begin()), 
                std::make_move_iterator(tmpMappings.end()));
            tmpMappings.clear();
            subrange_begin = subrange_end;
          }
        }
        //Sort the mappings by query (then reference) position
        std::sort(
            filteredMappings.begin(), filteredMappings.end(),
            [](const MappingResult &a, const MappingResult &b) {
                return std::tie(a.queryStartPos, a.refSeqId, a.refStartPos) 
                  < std::tie(b.queryStartPos, b.refSeqId, b.refStartPos);
                //return std::tie(a.refSeqId, a.refStartPos, a.queryStartPos)
                    //< std::tie(b.refSeqId, b.refStartPos, b.queryStartPos);
            });
      }


      /**
       * @brief               main mapping function given an input read
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @return              output object containing the mappings
       */
      MapModuleOutput* mapModule (InputSeqProgContainer* input)
      {
        MapModuleOutput* output = new MapModuleOutput();

        //save query sequence name and length
        output->qseqName = input->seqName;
        output->qseqLen = input->len;
        bool split_mapping = true;
        std::vector<IntervalPoint> intervalPoints;
        std::vector<L1_candidateLocus_t> l1Mappings;
        MappingResultsVector_t l2Mappings;
        MappingResultsVector_t unfilteredMappings;
        int refGroup = param.skip_prefix ? this->getRefGroup(input->seqName) : -1;

        if(! param.split || input->len <= param.segLength)
        {
          QueryMetaData <MinVec_Type> Q;
          Q.seq = &(input->seq)[0u];
          Q.len = input->len;
          Q.fullLen = input->len;
          Q.seqCounter = input->seqCounter;
          Q.seqName = input->seqName;
          Q.refGroup = refGroup;

          //Map this sequence
          mapSingleQueryFrag(Q, intervalPoints, l1Mappings, l2Mappings);

          // save the output
          unfilteredMappings.insert(unfilteredMappings.end(), l2Mappings.begin(), l2Mappings.end());

          // indicate that we mapped full length
          split_mapping = false;

          input->progress.increment(input->len);
        }
        else  //Split read mapping
        {
          int noOverlapFragmentCount = input->len / param.segLength;

          //Map individual non-overlapping fragments in the read
          for (int i = 0; i < noOverlapFragmentCount; i++)
          {
            //Prepare fragment sequence object
            QueryMetaData <MinVec_Type> Q;
            Q.seq = &(input->seq)[0u] + i * param.segLength;
            Q.len = param.segLength;
            Q.fullLen = input->len;
            Q.seqCounter = input->seqCounter;
            Q.seqName = input->seqName;
            Q.refGroup = refGroup;

            intervalPoints.clear();
            l1Mappings.clear();
            l2Mappings.clear();

            //Map this fragment
            mapSingleQueryFrag(Q, intervalPoints, l1Mappings, l2Mappings);

            //Adjust query coordinates and length in the reported mapping
            std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){
                e.queryLen = input->len;
                e.queryStartPos = i * param.segLength;
                e.queryEndPos = i * param.segLength + Q.len;
                });

            // save the output
            unfilteredMappings.insert(unfilteredMappings.end(), l2Mappings.begin(), l2Mappings.end());
            input->progress.increment(param.segLength);
          }

          //Map last overlapping fragment to cover the whole read
          if (noOverlapFragmentCount >= 1 && input->len % param.segLength != 0)
          {
            //Prepare fragment sequence object
            QueryMetaData <MinVec_Type> Q;
            Q.seq = &(input->seq)[0u] + input->len - param.segLength;
            Q.len = param.segLength;
            Q.seqCounter = input->seqCounter;
            Q.seqName = input->seqName;
            Q.refGroup = refGroup;

            intervalPoints.clear();
            l1Mappings.clear();
            l2Mappings.clear();

            //Map this fragment
            mapSingleQueryFrag(Q, intervalPoints, l1Mappings, l2Mappings);

            //Adjust query coordinates and length in the reported mapping
            std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){
                e.queryLen = input->len;
                e.queryStartPos = input->len - param.segLength;
                e.queryEndPos = input->len;
                });

            unfilteredMappings.insert(unfilteredMappings.end(), l2Mappings.begin(), l2Mappings.end());

            input->progress.increment(input->len % param.segLength);
          }
        }

        // how many mappings to keep
        int n_mappings = (input->len < param.segLength ?
                          param.numMappingsForShortSequence
                          : param.numMappingsForSegment) - 1;

        if (split_mapping) 
        {
          if (param.mergeMappings) 
          {
            // hardcore merge using the chain gap
            mergeMappingsInRange(unfilteredMappings, param.chain_gap);
            //mergeMappings(unfilteredMappings);

            // remove short chains that didn't exceed block length
            filterWeakMappings(unfilteredMappings, std::floor(param.block_length / param.segLength));
          }
        }

        if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {                      
          MappingResultsVector_t tmpMappings;
          tmpMappings.reserve(output->readMappings.size());
          filterByGroup(unfilteredMappings, tmpMappings, n_mappings, false);
          unfilteredMappings = std::move(tmpMappings);
        }

        output->readMappings = std::move(unfilteredMappings);

        //Make sure mapping boundary don't exceed sequence lengths
        this->mappingBoundarySanityCheck(input, output->readMappings);

        // remove alignments where the ratio between query and target length is < our identity threshold
        if (param.filterLengthMismatches)
        {
          this->filterFalseHighIdentity(output->readMappings);
        }

        // sparsify the mappings, if requested
        this->sparsifyMappings(output->readMappings);

        return output;
      }

      /**
       * @brief                       routine to handle mapModule's output of mappings
       * @param[in] output            mapping output object
       * @param[in] allReadMappings   vector to store mappings of all reads (optional use depending on filter)
       * @param[in] totalReadsMapped  counter to track count of reads mapped
       * @param[in] outstrm           outstream stream object
       */
      template <typename Vec>
      void mapModuleHandleOutput(MapModuleOutput* output,
                                 Vec &allReadMappings,
                                 seqno_t &totalReadsMapped,
                                 std::ofstream &outstrm,
                                 progress_meter::ProgressMeter& progress)
        {
          if(output->readMappings.size() > 0)
            totalReadsMapped++;

          if (param.filterMode == filter::ONETOONE)
          {
            //Save for another filtering round
            allReadMappings.insert(allReadMappings.end(), output->readMappings.begin(), output->readMappings.end());
          }
          else
          {
            //Report mapping
            reportReadMappings(output->readMappings, output->qseqName, outstrm);
          }

          //progress.increment(output->qseqLen/2 + (output->qseqLen % 2 != 0));

          delete output;
        }

      /**
       * @brief                   map the parsed query sequence (L1 and L2 mapping)
       * @param[in]   Q           metadata about query sequence
       * @param[in]   outstrm     outstream stream where mappings will be reported
       * @param[out]  l2Mappings  Mapping results in the L2 stage
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
              int currGroup = this->refIdGroup[l1_begin->seqId];
              l1_end = std::find_if_not(l1_begin, l1Mappings.end(), [this, currGroup] (const auto& candidate) {
                  return currGroup == this->refIdGroup[candidate.seqId];
              });
            }
            else
            {
              l1_end = l1Mappings.end();
            }

            //Sort L1 windows based on intersection size if using hg filter
            if (param.stage1_topANI_filter)
            {
              std::make_heap(l1_begin, l1_end, L1_locus_intersection_cmp);
            }
            doL2Mapping(Q, l1_begin, l1_end, l2Mappings);

            // Set beginning of next range
            l1_begin = l1_end;
          }

          // Sort output mappings
          std::sort(l2Mappings.begin(), l2Mappings.end(), [](const auto& a, const auto& b) 
              { return std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos); });

#ifdef ENABLE_TIME_PROFILE_L1_L2
          {
            std::chrono::duration<double> timeSpentL2 = skch::Time::now() - t1;
            std::chrono::duration<double> timeSpentMappingFragment = skch::Time::now() - t0;

            std::cerr << Q.seqCounter << " " << Q.len
              << " " << timeSpentL1.count()
              << " " << timeSpentL2.count()
              << " " << timeSpentMappingFragment.count()
              << "\n";
          }
#endif
        }

      template <typename Q_Info>
        void getSeedHits(Q_Info &Q)
        {
          Q.minmerTableQuery.reserve(param.sketchSize + 1);
          CommonFunc::sketchSequence(Q.minmerTableQuery, Q.seq, Q.len, param.kmerSize, param.alphabetSize, param.sketchSize, Q.seqCounter);
          if(Q.minmerTableQuery.size() == 0) {
            Q.sketchSize = 0;
            return;
          }

#ifdef DEBUG
          int orig_len = Q.minmerTableQuery.size();
#endif
          const double max_hash_01 = (long double)(Q.minmerTableQuery.back().hash) / std::numeric_limits<hash_t>::max();
          Q.kmerComplexity = (double(Q.minmerTableQuery.size()) / max_hash_01) / ((Q.len - param.kmerSize + 1)*2);

          // TODO remove them from the original sketch instead of removing for each read
          auto new_end = std::remove_if(Q.minmerTableQuery.begin(), Q.minmerTableQuery.end(), [&](auto& mi) {
            return refSketch.isFreqSeed(mi.hash);
          });
          Q.minmerTableQuery.erase(new_end, Q.minmerTableQuery.end());

          Q.sketchSize = Q.minmerTableQuery.size();
#ifdef DEBUG
          std::cerr << "INFO, skch::Map::getSeedHits, read id " << Q.seqCounter << ", minmer count = " << Q.minmerTableQuery.size() << ", bad minmers = " << orig_len - Q.sketchSize << "\n";
#endif
        } 


      /**
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is 
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an 
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details 
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename Vec>
        void getSeedIntervalPoints(Q_Info &Q, Vec& intervalPoints)
        {

#ifdef DEBUG
          std::cerr<< "INFO, skch::Map::getSeedHits, read id " << Q.seqCounter << ", minmer count = " << Q.minmerTableQuery.size() << " " << Q.len << "\n";
#endif

          //For invalid query (example : just NNNs), we may be left with 0 sketch size
          //Ignore the query in this case
          if(Q.minmerTableQuery.size() == 0)
            return;

          // Reserve the "expected" number of interval points
          // (lazy: the packed path never reaches this function)
          if (intervalPoints.capacity() == 0)
            intervalPoints.reserve(
                2 * param.sketchSize * refSketch.minmerIndex.size() / refSketch.nUniqueMinmers);

          // Gather matched interval points directly during the reference lookup
          // (no separate priority-queue pass; radixSortIntervalPoints sorts afterwards).
          const size_t ip_start = intervalPoints.size();
          const seqno_t queryRefId = param.skip_self ? this->queryRefSeqId(Q.seqName) : (seqno_t)-1;
          const bool doSelf = param.skip_self;
          const bool doPref = param.skip_prefix;
          const bool doLT = param.lower_triangular;
          const bool anyPairs = !allowed_pairs.empty();
          const auto* refGroupData = this->refIdGroup.data();
          const std::string pairPrefix = anyPairs ? Q.seqName + "\t" : std::string();
          for(auto it = Q.minmerTableQuery.begin(); it != Q.minmerTableQuery.end(); it++)
          {
            //Check if hash value exists in the reference lookup index
            if (refSketch.packed_ok) {
              // CSR arena: decode each packed key back to the identical
              // IntervalPoint (hash is the lookup key itself).
              const auto seedFind = refSketch.posLookupCSR.find(it->hash);
              if(seedFind == refSketch.posLookupCSR.end())
                continue;
              const uint64_t* runB = refSketch.ipArena.data() + seedFind->second.off;
              const uint64_t* runE = runB + seedFind->second.cnt;
              for (const uint64_t* k = runB; k != runE; ++k)
              {
                IntervalPoint ip = decodePackedIP(*k);
                ip.hash = it->hash;
                if ((!doSelf || queryRefId != ip.seqId)
                    && (!doPref || refGroupData[ip.seqId] != Q.refGroup)
                    && (!doLT || Q.seqCounter > ip.seqId)
                    && (!anyPairs
                        || allowed_pairs.count(pairPrefix + this->refSketch.metadata[ip.seqId].name))
                ) {
                  intervalPoints.push_back(ip);
                }
              }
              continue;
            }
            const auto seedFind = refSketch.minmerPosLookupIndex.find(it->hash);
            if(seedFind == refSketch.minmerPosLookupIndex.end())
              continue;

            for (const auto& ip : seedFind->second)
            {
              if ((!doSelf || queryRefId != ip.seqId)
                  && (!doPref || refGroupData[ip.seqId] != Q.refGroup)
                  && (!doLT || Q.seqCounter > ip.seqId)
                  && (!anyPairs
                      || allowed_pairs.count(pairPrefix + this->refSketch.metadata[ip.seqId].name))
              ) {
                intervalPoints.push_back(ip);
              }
            }
          }
          radixSortIntervalPoints(intervalPoints, ip_start);

#ifdef DEBUG
          std::cerr << "INFO, skch::Map:getSeedHits, read id " << Q.seqCounter << ", Count of seed hits in the reference = " << intervalPoints.size() / 2 << "\n";
#endif
        }

      /**
       * @brief  Packed-key variant of getSeedIntervalPoints (IP-1): produces the same
       *         filtered, order-preserving-key-sorted set of interval points as
       *         encodePackedIP(uint64) instead of 24-byte IntervalPoint structs.
       *         Used only when windowLen == 0, where IntervalPoint::hash is unused.
       */
      template <typename Q_Info>
        void getSeedIntervalPointsPacked(Q_Info &Q, std::vector<uint64_t>& packed)
        {
          if(Q.minmerTableQuery.size() == 0)
            return;

          const std::size_t start = packed.size();
          const seqno_t queryRefId = param.skip_self ? this->queryRefSeqId(Q.seqName) : (seqno_t)-1;
          const bool doSelf = param.skip_self;
          const bool doPref = param.skip_prefix;
          const bool doLT = param.lower_triangular;
          const bool anyPairs = !allowed_pairs.empty();
          const auto* refGroupData = this->refIdGroup.data();
          const std::string pairPrefix = anyPairs ? Q.seqName + "\t" : std::string();
          // When no per-point predicate can reject anything, whole runs are
          // appended with one insert (the arena is already encodePackedIP keys).
          const bool noFilter = (!doSelf || queryRefId == (seqno_t)-1) && !doPref && !doLT && !anyPairs;
          const uint64_t* arena = refSketch.ipArena.data();
          for(auto it = Q.minmerTableQuery.begin(); it != Q.minmerTableQuery.end(); it++)
          {
            const auto seedFind = refSketch.posLookupCSR.find(it->hash);
            if(seedFind == refSketch.posLookupCSR.end())
              continue;

            const uint64_t* runB = arena + seedFind->second.off;
            const uint64_t* runE = runB + seedFind->second.cnt;
            if (noFilter) {
              packed.insert(packed.end(), runB, runE);
              continue;
            }
            for (const uint64_t* k = runB; k != runE; ++k)
            {
              const seqno_t sid = (seqno_t)(*k >> 33);
              if ((!doSelf || queryRefId != sid)
                  && (!doPref || refGroupData[sid] != Q.refGroup)
                  && (!doLT || Q.seqCounter > sid)
                  && (!anyPairs
                      || allowed_pairs.count(pairPrefix + this->refSketch.metadata[sid].name))
              ) {
                packed.push_back(*k);
              }
            }
          }
          radixSortPackedKeys(packed, start);
        }


      // One recorded pos-group of the fused L1 counting sweep: the group's
      // coordinate (seqId of its first point, pos) and the overlap count right
      // after consuming the group.
      struct SweepStep {
        seqno_t seqId;
        offset_t pos;
        int overlapAfter;
      };

      // End of the run of interval points with seqId == sid starting at runStart.
      template <typename IP_iter>
      static IP_iter findIPRunEnd(IP_iter runStart, IP_iter ip_end, seqno_t sid)
      {
        if constexpr (std::is_same_v<IP_iter, PackedIPCursor>) {
          const uint64_t bound = (uint64_t)(uint32_t)(sid + 1) << 33;
          return PackedIPCursor{ std::lower_bound(runStart.p, ip_end.p, bound) };
        } else {
          return std::partition_point(runStart, ip_end,
              [sid](const auto& p) { return p.seqId <= sid; });
        }
      }

      template <typename IP_iter>
      static bool ipBefore(const IP_iter& a, const IP_iter& b)
      {
        if constexpr (std::is_same_v<IP_iter, PackedIPCursor>) {
          return a.p < b.p;
        } else {
          return a < b;
        }
      }

      template <typename Q_Info, typename IP_iter, typename Vec2>
        void computeL1CandidateRegions(
            Q_Info &Q, 
            IP_iter ip_begin, 
            IP_iter ip_end, 
            int minimumHits, 
            Vec2 &l1Mappings)
        {
#ifdef DEBUG
          std::cerr << "INFO, skch::Map:computeL1CandidateRegions, read id " << Q.seqCounter << std::endl;
#endif

          if (ip_begin == ip_end)
            return;

          int overlapCount = 0;
          int bestIntersectionSize = 0;
          thread_local std::vector<L1_candidateLocus_t> localOpts;
          localOpts.clear();

          // Keep track of all minmer windows that intersect with [i, i+windowLen]
          int windowLen = std::max<offset_t>(0, Q.len - param.segLength);
          auto trailingIt = ip_begin;
          auto leadingIt = ip_begin;

          // Group together local sketch intersection maximums that are within clusterLen of eachother
          //
          // Since setting up the L2 window [i, j] requires aggregating minmer windows over
          // [i-segLength, i), we might as well group L2 windows together which are closer than
          // segLength  
          int clusterLen = param.segLength;

          // Used to keep track of how many minmer windows for a particular hash are currently "open"
          // Only necessary when windowLen != 0.
          std::unordered_map<hash_t, int> hash_to_freq;

          bool in_candidate = false;
          L1_candidateLocus_t l1_out = {};

          // Candidate-emission state machine, applied once per pos-group with the
          // overlap count and coordinate of the PREVIOUS group. Shared verbatim by
          // the fused replay (stage1_topANI_filter) and the plain sweep below.
          auto emit = [&](int prevOverlap, const SeqCoord& prevPos) {
          if ( prevOverlap >= minimumHits
              //&& prevOverlap > overlapCount && prevOverlap >= prevPrevOverlap)
          ) {
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
            // Fused counting sweep: one traversal records each pos-group's
            // coordinate and post-group overlap; the emission machine is replayed
            // from that compact record after minimumHits is raised. seqId runs
            // that cannot reach minimumHits collapse to a single zero-overlap
            // step: their opens/closes balance to zero, they can never emit, and
            // they cannot own bestIntersectionSize in a way that changes the
            // early return or the raise (their best < minimumHits <= any
            // qualifying run's best).
            thread_local std::vector<SweepStep> steps;
            steps.clear();
            auto runStart = ip_begin;
            bool cleanStart = true;   // leading did not overshoot into this run
            while (runStart != ip_end)
            {
              const seqno_t runSeqId = runStart->seqId;
              const auto runEnd = findIPRunEnd(runStart, ip_end, runSeqId);
              if (windowLen == 0 && cleanStart)
              {
                std::size_t runLen;
                offset_t lastPos;
                if constexpr (std::is_same_v<IP_iter, PackedIPCursor>) {
                  runLen = (std::size_t)(runEnd.p - runStart.p);
                  lastPos = (offset_t)((*(runEnd.p - 1) >> 1) & 0xFFFFFFFFULL);
                } else {
                  runLen = (std::size_t)std::distance(runStart, runEnd);
                  lastPos = std::prev(runEnd)->pos;
                }
                // Prune only when the run also does not share a position with the
                // next run (pos-only grouping would merge such points into one
                // group whose overlap must include both runs' opens).
                if ((int)(runLen / 2) < minimumHits
                    && !(runEnd != ip_end && runEnd->pos == lastPos))
                {
                  steps.push_back(SweepStep{runSeqId, runStart->pos, 0});
                  trailingIt = runEnd;
                  leadingIt = runEnd;
                  runStart = runEnd;
                  continue;
                }
              }
              while (ipBefore(leadingIt, runEnd))
              {
                // Catch the trailing iterator up to the leading iterator - windowLen
                while (
                    trailingIt != ip_end
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
                const seqno_t groupSeqId = leadingIt->seqId;
                const offset_t groupPos = leadingIt->pos;
                while (leadingIt != ip_end && leadingIt->pos == groupPos) {
                  if (leadingIt->side == side::OPEN) {
                    if (windowLen == 0 || hash_to_freq[leadingIt->hash] == 0) {
                      overlapCount++;
                    }
                    if (windowLen != 0)
                      hash_to_freq[leadingIt->hash]++;
                  }
                  leadingIt++;
                }
                //Is this sliding window the best we have so far?
                bestIntersectionSize = std::max(bestIntersectionSize, overlapCount);
                steps.push_back(SweepStep{groupSeqId, groupPos, overlapCount});
              }
              cleanStart = (leadingIt == runEnd);
              runStart = leadingIt;
            }

            // Only go back through to find local opts if we know that there are some that are 
            // large enough
            if (bestIntersectionSize < minimumHits) 
            {
              return;
            } else 
            {
              minimumHits = std::max(
                  sketchCutoffs[
                    int(std::min(bestIntersectionSize, Q.sketchSize)
                      / std::max<double>(1, param.sketchSize / skch::fixed::ss_table_max))
                  ],
                  minimumHits);
            }

            // Replay the emission machine over the recorded pos-groups: step t's
            // emission sees the overlap/coordinate of step t-1 (the last step's
            // overlap is never read, exactly like the plain sweep's last group).
            int prevOverlap = 0;
            SeqCoord prevPos{};
            for (const auto& st : steps) {
              emit(prevOverlap, prevPos);
              prevOverlap = st.overlapAfter;
              prevPos = SeqCoord{st.seqId, st.pos};
            }
          }
          
#ifdef WFMASH_SWEEP_VERIFY
          // Verification mode: finalize and stash the fused result, then let the
          // reference sweep below recompute into localOpts for comparison.
          std::vector<L1_candidateLocus_t> fusedOpts;
          if (param.stage1_topANI_filter) {
            if (in_candidate) localOpts.push_back(l1_out);
            in_candidate = false;
            l1_out = {};
            fusedOpts = localOpts;
            localOpts.clear();
          }
          if (true)
#else
          if (!param.stage1_topANI_filter)
#endif
          {
          // Clear freq dict, as there will be left open CLOSE points at the end of the last seq
          // that we never got to
          hash_to_freq.clear();

          // Since there can be more than sketchSize windows that overlap w/ [i, i+windowLen]
          // cap the best intersection size
          bestIntersectionSize = std::min(bestIntersectionSize, Q.sketchSize);

          trailingIt = ip_begin;
          leadingIt = ip_begin;

          // Keep track of 3 consecutive points so that we can track local optimums
          overlapCount = 0;
          int prevOverlap = 0;
          int prevPrevOverlap = 0;

          // Need to keep track of two positions, as the previous one will be the local optimum
          SeqCoord prevPos;
          SeqCoord currentPos{leadingIt->seqId, leadingIt->pos};


          while (leadingIt != ip_end)
          {
            prevPrevOverlap = prevOverlap;
            prevOverlap = overlapCount;

            //TODO LEADING it should only hit opens
            // We should only iterate through a new window when we come across an OPEN,
            // right now, this basically happens since every CLOSE should be an OPEN.
            // This doesn't invalidate the logic, just potentially wastes time
            while (
                trailingIt != ip_end 
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
#ifdef WFMASH_SWEEP_VERIFY
        if (param.stage1_topANI_filter) {
          bool same = fusedOpts.size() == localOpts.size();
          for (std::size_t i = 0; same && i < fusedOpts.size(); ++i) {
            same = fusedOpts[i].seqId == localOpts[i].seqId
                && fusedOpts[i].rangeStartPos == localOpts[i].rangeStartPos
                && fusedOpts[i].rangeEndPos == localOpts[i].rangeEndPos
                && fusedOpts[i].intersectionSize == localOpts[i].intersectionSize;
          }
          if (!same) {
            std::cerr << "[wfmash] WFMASH_SWEEP_VERIFY mismatch, query " << Q.seqCounter
                      << " fused=" << fusedOpts.size() << " ref=" << localOpts.size() << std::endl;
            std::abort();
          }
        }
#endif
        

        // Join together proximal local opts
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
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename IPVec, typename L1Vec>
        void doL1Mapping(Q_Info &Q, IPVec& intervalPoints, L1Vec& l1Mappings)
        {
          //1. Compute the minmers
          getSeedHits(Q);

          //Catch all NNNNNN case
          if (Q.sketchSize == 0 || Q.kmerComplexity < param.kmerComplexityThreshold) {
            return;
          }

          //3. Compute L1 windows
          int minimumHits = (size_t)Q.sketchSize < minHitsCache.size()
              ? minHitsCache[Q.sketchSize]
              : Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity, skch::fixed::confidence_interval);

          // Fast packed-key path: windowLen == max(0, Q.len - segLength) is 0 here
          // (Q.len <= segLength, the default split fragmentation), so IntervalPoint::hash
          // is unused and interval points can be compact uint64 keys -- no 24-byte struct,
          // no scattered gather in the sort. Byte-identical to the struct path below.
          if (this->packed_ip_ok && Q.len <= param.segLength)
          {
            thread_local std::vector<uint64_t> packedPoints;
            packedPoints.clear();
            getSeedIntervalPointsPacked(Q, packedPoints);

            const std::size_t np = packedPoints.size();
            std::size_t b = 0;
            while (b < np)
            {
              std::size_t e;
              if (param.skip_prefix)
              {
                const int currGroup = this->refIdGroup[(seqno_t)(packedPoints[b] >> 33)];
                e = b;
                while (e < np && this->refIdGroup[(seqno_t)(packedPoints[e] >> 33)] == currGroup) ++e;
              }
              else
              {
                e = np;
              }
              computeL1CandidateRegions(Q, PackedIPCursor{packedPoints.data() + b},
                                            PackedIPCursor{packedPoints.data() + e},
                                            minimumHits, l1Mappings);
              b = e;
            }
            return;
          }

          //2. Compute windows and sort (struct path; also handles windowLen != 0)
          getSeedIntervalPoints(Q, intervalPoints);

          // For each "group"
          auto ip_begin = intervalPoints.begin();
          auto ip_end = intervalPoints.begin();
          while (ip_end != intervalPoints.end())
          {
            if (param.skip_prefix)
            {
              int currGroup = this->refIdGroup[ip_begin->seqId];
              ip_end = std::find_if_not(ip_begin, intervalPoints.end(), [this, currGroup] (const auto& ip) {
                  return currGroup == this->refIdGroup[ip.seqId];
              });
            }
            else
            {
              ip_end = intervalPoints.end();
            }
            computeL1CandidateRegions(Q, ip_begin, ip_end, minimumHits, l1Mappings);

            ip_begin = ip_end;
          }
        }


      // helper to get the prefix of a string
      const std::string prefix(const std::string& s, const char c) {
          //std::cerr << "prefix of " << s << " by " << c << " is " << s.substr(0, s.find_last_of(c)) << std::endl;
          return s.substr(0, s.find_last_of(c));
      }

      /**
       * @brief                                 Revise L1 candidate regions to more precise locations
       * @param[in]   Q                         query sequence information
       * @param[in]   l1Mappings                candidate regions for query sequence found at L1
       * @param[out]  l2Mappings                Mapping results in the L2 stage
       */
      template <typename Q_Info, typename L1_Iter, typename VecOut>
        void doL2Mapping(Q_Info &Q, L1_Iter l1_begin, L1_Iter l1_end, VecOut &l2Mappings)
        {
          ///2. Walk the read over the candidate regions and compute the jaccard similarity with minimum s sketches
          thread_local std::vector<L2_mapLocus_t> l2_vec;
          double bestJaccardNumerator = 0;
          auto loc_iterator = l1_begin;
          while (loc_iterator != l1_end)
          {
            L1_candidateLocus_t& candidateLocus = *loc_iterator;

            if (param.stage1_topANI_filter)
            {
              // If using HG filter, don't consider any mappings which have no chance of being 
              // within param.ANIDiff of the best mapping seen so far
              double cutoff_ani = std::max(0.0, double((1 - Stat::j2md(bestJaccardNumerator / Q.sketchSize, param.kmerSize)) - param.ANIDiff));
              double cutoff_j = Stat::md2j(1 - cutoff_ani, param.kmerSize);
              if (double(candidateLocus.intersectionSize) / Q.sketchSize < cutoff_j) 
              {
                break;
              }
            }


            l2_vec.clear();
            computeL2MappedRegions(Q, candidateLocus, l2_vec);

            for (auto& l2 : l2_vec) 
            {
              //Compute mash distance using calculated jaccard
              float mash_dist = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);

              float nucIdentity = (1 - mash_dist);
              //float nucIdentityUpperBound = getANIUBfromJaccardNum(Q.sketchSize, l2.sharedSketchSize);
              float nucIdentityUpperBound = 1 - Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, skch::fixed::confidence_interval);

              //Report the alignment if it passes our identity threshold and,
              // if we are in all-vs-all mode, it isn't a self-mapping,
              // and if we are self-mapping, the query is shorter than the target
              const auto& ref = this->refSketch.metadata[l2.seqId];
              if((param.keep_low_pct_id && nucIdentityUpperBound >= param.percentageIdentity)
                  || nucIdentity >= param.percentageIdentity)
              {
                //Track the best jaccard numerator
                bestJaccardNumerator = std::max<double>(bestJaccardNumerator, l2.sharedSketchSize);

                MappingResult res;

                //Save the output
                {
                  res.queryLen = Q.len;
                  res.refStartPos = l2.meanOptimalPos;
                  res.refEndPos = l2.meanOptimalPos + Q.len;
                  res.queryStartPos = 0;
                  res.queryEndPos = Q.len;
                  res.refSeqId = l2.seqId;
                  res.querySeqId = Q.seqCounter;
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
              l1_end--; //"Pop back" 
            }
            else 
            {
              loc_iterator++;
            }
          }
          //std::cerr << "For an segment with " << l1Mappings.size()
            //<< " L1 mappings "
            //<< " there were " << l2Mappings.size() << " L2 mappings\n";
        }

      /**
       * @brief                                 Find optimal mapping within an L1 candidate
       * @param[in]   Q                         query sequence information
       * @param[in]   candidateLocus            L1 candidate location
       * @param[out]  l2_out                    L2 mapping inside L1 candidate
       */
      template <typename Q_Info, typename Vec>
        void computeL2MappedRegions(Q_Info &Q,
            L1_candidateLocus_t &candidateLocus,
            Vec &l2_vec_out)
        {
#ifdef DEBUG
          //std::cerr << "INFO, skch::Map:computeL2MappedRegions, read id " << Q.seqName << "_" << Q.startPos << std::endl; 
#endif
           
          auto& minmerIndex = refSketch.minmerIndex;

          //candidateLocus.rangeStartPos -= param.segLength;
          //candidateLocus.rangeEndPos += param.segLength;
          
          // Get first potential mashimizer: search only within this seqId's slice
          // of the (seqId, wpos)-sorted index (equivalent to the former whole-index
          // lower_bound on {seqId, wpos}).
          const offset_t windowStartPos = candidateLocus.rangeStartPos - param.segLength - 1;
          auto firstOpenIt = std::lower_bound(
              minmerIndex.begin() + minmerIndexSeqStart[candidateLocus.seqId],
              minmerIndex.begin() + minmerIndexSeqStart[candidateLocus.seqId + 1],
              windowStartPos,
              [](const MinmerInfo& mi, offset_t w) { return mi.wpos < w; });

          // Keeps track of the lowest end position
          thread_local std::vector<skch::MinmerInfo> slidingWindow;
          slidingWindow.clear();
          slidingWindow.reserve(Q.sketchSize);

          // Used to make a min-heap
          constexpr auto heap_cmp = [](const skch::MinmerInfo& l, const skch::MinmerInfo& r) {return l.wpos_end > r.wpos_end;};

          // windowIt keeps track of the end of window
          auto windowIt = firstOpenIt;

          // Keep track of all minmer windows that intersect with [i, i+windowLen]
          int windowLen = std::max<offset_t>(0, Q.len - param.segLength);

          // Used to keep track of how many minmer windows for a particular hash are currently "open"
          // Only necessary when windowLen != 0.
          std::unordered_map<hash_t, int> hash_to_freq;
          
          // slideMap tracks the S(A or B) and S(A) and S(B)
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

          while (windowIt != minmerIndex.end() && windowIt->seqId == candidateLocus.seqId && windowIt->wpos <= candidateLocus.rangeEndPos + windowLen) 
          {
            int prev_strand_votes = slideMap.strand_votes;
            bool inserted = false;
            while (!slidingWindow.empty() && slidingWindow.front().wpos_end <= windowIt->wpos - windowLen) {

              // Remove minmer from end-ordered heap
              if (windowLen > 0) 
              {
                hash_to_freq[slidingWindow.front().hash]--;
              }
              if (windowLen == 0 || hash_to_freq[slidingWindow.front().hash] == 0) {
                // Remove minmer from  sorted window
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

            //Is this sliding window the best we have so far?
            if (slideMap.sharedSketchElements > bestSketchSize)
            {
              // Get rid of all candidates seen so far
              l2_vec_out.clear();

              in_candidate = true;
              bestSketchSize = slideMap.sharedSketchElements;
              l2_out.sharedSketchSize = slideMap.sharedSketchElements;

              //Save the position
              l2_out.optimalStart = windowIt->wpos - windowLen;
              l2_out.optimalEnd = windowIt->wpos - windowLen;
            }
            else if(slideMap.sharedSketchElements == bestSketchSize)
            {
              if (!in_candidate) {
                l2_out.sharedSketchSize = slideMap.sharedSketchElements;

                //Save the position
                l2_out.optimalStart = windowIt->wpos - windowLen;
              }

              in_candidate = true;
              //Still save the position
              l2_out.optimalEnd = windowIt->wpos - windowLen;
            } else {
              if (in_candidate) {
                // Save and reset
                l2_out.meanOptimalPos =  (l2_out.optimalStart + l2_out.optimalEnd) / 2;
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
            // Save and reset
            l2_out.meanOptimalPos =  (l2_out.optimalStart + l2_out.optimalEnd) / 2;
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


      /**
       * @brief                       Merge the consecutive fragment mappings reported in each query
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void expandMappings(VecIn &readMappings, int expansion)
        {
            for (auto& m : readMappings) {
                m.refStartPos -= expansion;
                m.refEndPos += expansion;
                m.queryStartPos -= expansion;
                m.queryEndPos += expansion;
            }
        }

      /**
       * @brief                       Merge the consecutive fragment mappings reported in each query
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void mergeMappings(VecIn &readMappings)
        {
          assert(param.split == true);

          if(readMappings.size() < 2)
            return;

          //Sort the mappings by reference position
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)
              {
              return std::tie(a.refSeqId, a.refStartPos, a.queryStartPos) < std::tie(b.refSeqId, b.refStartPos, b.queryStartPos);
              });

          //First assign a unique id to each split mapping in the sorted order
          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            it->splitMappingId = std::distance(readMappings.begin(), it);
          }

          //Start the procedure to identify the chains
          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            //Which fragment is this wrt. the complete read
            auto currMappingFragno = std::ceil(it->queryStartPos * 1.0/param.segLength);

            for(auto it2 = std::next(it); it2 != readMappings.end(); it2++)
            {
              auto thisMappingFragno = std::ceil(it2->queryStartPos * 1.0/ param.segLength);

              //If this mapping is too far from current mapping being evaluated, stop finding a merge
              if(
                  it2->refSeqId != it->refSeqId 
                  || std::abs(it2->refStartPos - it->refEndPos) > param.chain_gap
                  )
                break;

              //If the next mapping is within range, check if it is consecutive query fragment and strand matches
              if( it2->strand == it->strand
                  //&& std::abs(it2->queryStartPos - it->queryEndPos) <= param.chain_gap
                  && thisMappingFragno == currMappingFragno + (it->strand == strnd::FWD ? 1 : -1)
              )
              {
                it2->splitMappingId = it->splitMappingId;   //merge
                continue;
              }
            }
          }
          //Keep single mapping for each chain and discard others

          //Sort the mappings by post-merge split mapping id
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)
              {
              return a.splitMappingId < b.splitMappingId;
              });

          for(auto it = readMappings.begin(); it != readMappings.end();)
          {
            //Bucket by each chain
            auto it_end = std::find_if(it, readMappings.end(), [&](const MappingResult &e){return e.splitMappingId != it->splitMappingId;} );

            //[it -- it_end) represents same chain

            //Incorporate chain information into first mapping

            //compute chain length
            std::for_each(it, it_end, [&](MappingResult &e)
            {
              it->queryStartPos = std::min( it->queryStartPos, e.queryStartPos);
              it->refStartPos = std::min( it->refStartPos, e.refStartPos);

              it->queryEndPos = std::max( it->queryEndPos, e.queryEndPos);
              it->refEndPos = std::max( it->refEndPos, e.refEndPos);

              it->blockLength = std::max(it->refEndPos - it->refStartPos, it->queryEndPos - it->queryStartPos);
              it->approxMatches = std::round(it->nucIdentity * it->blockLength / 100.0);
            });

            it->n_merged = std::distance(it, it_end);

            //Mean identity of all mappings in the chain
            it->nucIdentity = (   std::accumulate(it, it_end, 0.0,
                                  [](double x, MappingResult &e){ return x + e.nucIdentity; })     )/ std::distance(it, it_end);

            //Mean sequence complexity of all mappings in the chain
            it->kmerComplexity = (   std::accumulate(it, it_end, 0.0,
                                  [](double x, MappingResult &e){ return x + e.kmerComplexity; })     )/ std::distance(it, it_end);

            //Discard other mappings of this chain
            std::for_each( std::next(it), it_end, [&](MappingResult &e){ e.discard = 1; });

            //advance the iterator
            it = it_end;
          }

          readMappings.erase(
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());
       }


      /**
       * @brief                       Merge fragment mappings by convolution of a 2D range over the alignment matrix
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       * @param[in]     max_dist      Distance to look in target and query
       */
      template <typename VecIn>
      void mergeMappingsInRange(VecIn &readMappings,
                                int max_dist) {
          assert(param.split == true);

          if(readMappings.size() < 2) return;

          //Sort the mappings by reference (then query) position
          std::sort(
              readMappings.begin(), readMappings.end(),
              [](const MappingResult &a, const MappingResult &b) {
                  return std::tie(a.refSeqId, a.refStartPos, a.queryStartPos)
                      < std::tie(b.refSeqId, b.refStartPos, b.queryStartPos);
              });

          //First assign a unique id to each split mapping in the sorted order
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              it->splitMappingId = std::distance(readMappings.begin(), it);
              it->discard = 0;
          }

          // set up our union find data structure to track merges
          std::vector<dsets::DisjointSets::Aint> ufv(readMappings.size());
          // this initializes everything
          auto disjoint_sets = dsets::DisjointSets(ufv.data(), ufv.size());

          //Start the procedure to identify the chains
          std::vector<std::pair<double, uint64_t>> distances;
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              distances.clear();
              for (auto it2 = std::next(it); it2 != readMappings.end(); it2++) {
                  //If this mapping is for the same segment, ignore
                  if (it2->refSeqId == it->refSeqId && it2->queryStartPos == it->queryStartPos) {
                    continue;
                  }
                  //If this mapping is too far from current mapping being evaluated, stop finding a merge
                  if (it2->refSeqId != it->refSeqId || it2->refStartPos > it->refEndPos + max_dist) {
                      break;
                  }
                  //If the next mapping is within range, check if it's in range and
                  if (it2->strand == it->strand) {
                      int ref_dist = it2->refStartPos - it->refEndPos;
                      int query_dist = 0;
                      auto dist = std::numeric_limits<double>::max();
                      auto score = std::numeric_limits<double>::max();
                      if (it->strand == strnd::FWD && it->queryStartPos <= it2->queryStartPos) {
                          query_dist = it2->queryStartPos - it->queryEndPos;
                          dist = std::sqrt(std::pow(query_dist,2) + std::pow(ref_dist,2));
                          score = std::pow(query_dist - ref_dist, 2);
                      } else if (it->strand != strnd::FWD && it->queryEndPos >= it2->queryEndPos) {
                          query_dist = it->queryStartPos - it2->queryEndPos;
                          dist = std::sqrt(std::pow(query_dist,2) + std::pow(ref_dist,2));
                          score = std::pow(query_dist - ref_dist, 2);
                      }
                      int query_mapping_len = std::min((it->queryEndPos - it->queryStartPos),
                                                       (it2->queryEndPos - it2->queryStartPos));
                      if (dist < max_dist) {
                          distances.push_back(std::make_pair(dist + score, it2->splitMappingId));
                      }
                  }
              }
              if (distances.size()) {
                  disjoint_sets.unite(it->splitMappingId,
                      std::min_element(distances.begin(), distances.end())->second);
              }
          }

          //Assign the merged mapping ids
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              it->splitMappingId = disjoint_sets.find(it->splitMappingId);
          }

          //Sort the mappings by post-merge split mapping id
          std::sort(
              readMappings.begin(),
              readMappings.end(),
              [](const MappingResult &a, const MappingResult &b) {
                  return a.splitMappingId < b.splitMappingId;
              });

          for(auto it = readMappings.begin(); it != readMappings.end();) {

              //Bucket by each chain
              auto it_end = std::find_if(it, readMappings.end(), [&](const MappingResult &e){return e.splitMappingId != it->splitMappingId;} );

              //std::cerr << "Got chain with " <<

              //[it -- it_end) represents same chain

              //Incorporate chain information into first mapping

              //compute chain length
              std::for_each(it, it_end, [&](MappingResult &e)
                  {
                      it->queryStartPos = std::min( it->queryStartPos, e.queryStartPos);
                      it->refStartPos = std::min( it->refStartPos, e.refStartPos);

                      it->queryEndPos = std::max( it->queryEndPos, e.queryEndPos);
                      it->refEndPos = std::max( it->refEndPos, e.refEndPos);

                      it->blockLength = std::max(it->refEndPos - it->refStartPos, it->queryEndPos - it->queryStartPos);
                      it->approxMatches = std::round(it->nucIdentity * it->blockLength / 100.0);
                  });

              it->n_merged = std::distance(it, it_end);

              //Mean identity of all mappings in the chain
              it->nucIdentity = ( std::accumulate(
                                      it, it_end, 0.0,
                                      [](double x, MappingResult &e){ return x + e.nucIdentity; })
                  ) / it->n_merged; // this would scale directly by the number of mappings in the chain

              //Mean identity of all kmer complexities in the chain
              it->kmerComplexity = ( std::accumulate(
                                      it, it_end, 0.0,
                                      [](double x, MappingResult &e){ return x + e.kmerComplexity; })
                  ) / it->n_merged; // this would scale directly by the number of mappings in the chain

              //Discard other mappings of this chain
              std::for_each( std::next(it), it_end, [&](MappingResult &e){ e.discard = 1; });

              //advance the iterator
              it = it_end;
          }

          readMappings.erase(
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());

      }

     /**
       * @brief                       This routine is to make sure that all mapping boundaries
       *                              on query and reference are not outside total
       *                              length of sequeunces involved
       * @param[in]     input         input read details
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void mappingBoundarySanityCheck(InputSeqProgContainer* input, VecIn &readMappings)
        {
          for(auto &e : readMappings)
          {
            //reference start pos
            {
              if(e.refStartPos < 0)
                e.refStartPos = 0;
              if(e.refStartPos >= this->refSketch.metadata[e.refSeqId].len)
                e.refStartPos = this->refSketch.metadata[e.refSeqId].len - 1;
            }

            //reference end pos
            {
              if(e.refEndPos < e.refStartPos)
                e.refEndPos = e.refStartPos;
              if(e.refEndPos >= this->refSketch.metadata[e.refSeqId].len)
                e.refEndPos = this->refSketch.metadata[e.refSeqId].len - 1;
            }

            //query start pos
            {
              if(e.queryStartPos < 0)
                e.queryStartPos = 0;
              if(e.queryStartPos >= input->len)
                e.queryStartPos = input->len;
            }

            //query end pos
            {
              if(e.queryEndPos < e.queryStartPos)
                e.queryEndPos = e.queryStartPos;
              if(e.queryEndPos >= input->len)
                e.queryEndPos = input->len;
            }
          }
        }

      /**
       * @brief                         Report the final read mappings to output stream
       * @param[in]   readMappings      mapping results for single or multiple reads
       * @param[in]   queryName         input required if reporting one read at a time
       * @param[in]   outstrm           file output stream object
       */
      void reportReadMappings(MappingResultsVector_t &readMappings, const std::string &queryName,
          std::ofstream &outstrm)
      {
        //Print the results
        for(auto &e : readMappings)
        {
          assert(e.refSeqId < this->refSketch.metadata.size());

          float fakeMapQ = e.nucIdentity == 1 ? 255 : std::round(-10.0 * std::log10(1-(e.nucIdentity)));
          std::string sep = param.legacy_output ? " " : "\t";

          outstrm  << (param.filterMode == filter::ONETOONE ? qmetadata[e.querySeqId].name : queryName)
                   << sep << e.queryLen
                   << sep << e.queryStartPos
                   << sep << e.queryEndPos - (param.legacy_output ? 1 : 0)
                   << sep << (e.strand == strnd::FWD ? "+" : "-")
                   << sep << this->refSketch.metadata[e.refSeqId].name
                   << sep << this->refSketch.metadata[e.refSeqId].len
                   << sep << e.refStartPos
                   << sep << e.refEndPos - (param.legacy_output ? 1 : 0);

          if (!param.legacy_output) 
          {
            outstrm  << sep << e.conservedSketches
                     << sep << e.blockLength
                     << sep << fakeMapQ
                     << sep << "id:f:" << e.nucIdentity
                     << sep << "kc:f:" << e.kmerComplexity;
            if (!param.mergeMappings) 
            {
              outstrm << sep << "jc:f:" << float(e.conservedSketches) / e.sketchSize;
            }
          } else
          {
            outstrm << sep << e.nucIdentity * 100.0;
          }

#ifdef DEBUG
          outstrm << std::endl;
#else
          outstrm << "\n";
#endif

          //User defined processing of the results
          if(processMappingResults != nullptr)
            processMappingResults(e);
        }
      }

    public:

      /**
       * @brief     An optional utility function to save the
       *            reported results by the L2 stage into a vector
       */
      static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
      {
        v.push_back(reportedL2Result);
      }

  };

}

#endif
