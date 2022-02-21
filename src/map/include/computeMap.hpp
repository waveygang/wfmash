/**
 * @file    computeMap.hpp
 * @brief   implements the sequence mapping logic
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_MAP_HPP 
#define SKETCH_MAP_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <zlib.h>  
#include <cassert>
#include <numeric>

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
#include "common/filesystem.hpp"
#include "robin-hood-hashing/robin_hood.h"
// if we ever want to do the union-find chaining in parallel
//#include "common/dset64-gccAtomic.hpp"
// this is for single-threaded use, but is more portable
#include "common/dset64.hpp"

namespace skch
{
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
      };

      //Type for Stage L2's predicted mapping coordinate within each L1 candidate
      struct L2_mapLocus_t 
      {
        seqno_t seqId;                    //sequence id where read is mapped
        offset_t meanOptimalPos;          //Among multiple consecutive optimal positions, save the avg.
        Sketch::MIIter_t optimalStart;    //optimal start mapping position (begin iterator)
        Sketch::MIIter_t optimalEnd;      //optimal end mapping position (end iterator) 
        int sharedSketchSize;             //count of shared sketch elements
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
        processMappingResults(f)
    {
      this->mapQuery();
    }

    private:

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
        ThreadPool<InputSeqContainer, MapModuleOutput> threadPool( [this](InputSeqContainer* e){return mapModule(e);}, param.threads);

        // kind've expensive if the fasta index is not available for the query sequences,
        // but it can help people know how long we're going to take
        uint64_t total_seqs = 0;
        uint64_t total_seq_length = 0;
        for(const auto &fileName : param.querySequences) {
            // check if there is a .fai
            std::string fai_name = fileName + ".fai";
            if (fs::file_exists(fai_name)) {
                // if so, process the .fai to determine our sequence length
                std::string line;
                std::ifstream in(fai_name.c_str());
                while (std::getline(in, line)) {
                    ++total_seqs;
                    auto line_split = CommonFunc::split(line, '\t');
                    total_seq_length += std::stoul(line_split[1]);
                }
            } else {
                // if not, warn that this is expensive
                std::cerr << "[wfmash::skch::Map::mapQuery] WARNING, no .fai index found for " << fileName << ", reading the file to sum sequence length (slow)" << std::endl;
                seqiter::for_each_seq_in_file(
                    fileName,
                    [&](const std::string& seq_name,
                        const std::string& seq) {
                        ++total_seqs;
                        total_seq_length += seq.size();
                    });
            }
        }

        progress_meter::ProgressMeter progress(total_seq_length, "[wfmash::skch::Map::mapQuery] mapped");

        for(const auto &fileName : param.querySequences)
        {

#ifdef DEBUG
            std::cerr << "[wfmash::skch::Map::mapQuery] mapping reads in " << fileName << std::endl;
#endif

            seqiter::for_each_seq_in_file(
                fileName,
                [&](const std::string& seq_name,
                    const std::string& seq) {
                    // todo: offset_t is an 32-bit integer, which could cause problems
                    offset_t len = seq.length();

                    if (param.filterMode == filter::ONETOONE)
                        qmetadata.push_back( ContigInfo{seq_name, len} );
                    //Is the read too short?
                    if(len < param.windowSize || len < param.kmerSize)
                    {
//#ifdef DEBUG
                        // TODO Should we somehow revert to < windowSize?
                        std::cerr << std::endl
                                  << "WARNING, skch::Map::mapQuery, read "
                                  << seq_name << " of " << len << "bp "
                                  << " is not long enough for mapping at window size "
                                  << param.windowSize << std::endl;
//#endif
                    }
                    else 
                    {
                        totalReadsPickedForMapping++;
                        //Dispatch input to thread
                        threadPool.runWhenThreadAvailable(new InputSeqContainer(seq, seq_name, seqCounter));

                        //Collect output if available
                        while ( threadPool.outputAvailable() ) {
                            mapModuleHandleOutput(threadPool.popOutputWhenAvailable(), allReadMappings, totalReadsMapped, outstrm, progress);
                        }
                    }
                    progress.increment(seq.size()/2);
                    seqCounter++;
                }); //Finish reading query input file

        }

        //Collect remaining output objects
        while ( threadPool.running() )
            mapModuleHandleOutput(threadPool.popOutputWhenAvailable(), allReadMappings, totalReadsMapped, outstrm, progress);

        //Filter over reference axis and report the mappings
        if (param.filterMode == filter::ONETOONE)
        {
          skch::Filter::ref::filterMappings(allReadMappings, this->refSketch,
                                            param.numMappingsForSegment - 1
                                           // (input->len < param.segLength ? param.shortSecondaryToKeep : param.secondaryToKeep)
                                            );

          //Re-sort mappings by input order of query sequences
          //This order may be needed for any post analysis of output
          std::sort(allReadMappings.begin(), allReadMappings.end(), [](const MappingResult &a, const MappingResult &b)  
          {
            return (a.querySeqId < b.querySeqId);
          });

          reportReadMappings(allReadMappings, "", outstrm);
        }

        progress.finish();

        std::cerr << "[wfmash::skch::Map::mapQuery] "
                  << "count of mapped reads = " << totalReadsMapped
                  << ", reads qualified for mapping = " << totalReadsPickedForMapping
                  << ", total input reads = " << seqCounter
                  << ", total input bp = " << total_seq_length << std::endl;

      }

      /**
       * @brief               helper to main mapping function
       * @details             filters mappings shorter than our minimum block length
       * @param[in]   input   mappings
       * @return              void
       */
      void filterShortMappings(MappingResultsVector_t &readMappings, int64_t length)
      {
          readMappings.erase(
              std::remove_if(readMappings.begin(),
                             readMappings.end(),
                             [&](MappingResult &e){
                                 return e.blockLength < length;
                             }),
              readMappings.end());
      }

      /**
       * @brief               helper to main mapping function
       * @details             filters long-to-short mappings if we're in an all-vs-all mode
       * @param[in]   input   mappings
       * @return              void
       */
      void filterSelfingLongToShorts(MappingResultsVector_t &readMappings)
      {
          if (param.skip_self || param.skip_prefix) {
              readMappings.erase(
                  std::remove_if(readMappings.begin(),
                                 readMappings.end(),
                                 [&](MappingResult &e){ return e.selfMapFilter == true; }),
                  readMappings.end());
          }
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
                                 return len_id_bound < param.percentageIdentity;
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
       * @brief               main mapping function given an input read
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @return              output object containing the mappings
       */
      MapModuleOutput* mapModule (InputSeqContainer* input)
      {
        MapModuleOutput* output = new MapModuleOutput();

        //save query sequence name and length
        output->qseqName = input->seqName;
        output->qseqLen = input->len;
        bool split_mapping = true;

        if(! param.split || input->len < param.segLength || input->len <= param.block_length_min)
        {
          QueryMetaData <MinVec_Type> Q;
          Q.seq = &(input->seq)[0u];
          Q.len = input->len;
          Q.fullLen = input->len;
          Q.seqCounter = input->seqCounter;
          Q.seqName = input->seqName;

          MappingResultsVector_t l2Mappings;   

          //Map this sequence
          mapSingleQueryFrag(Q, l2Mappings);

          // save the output
          output->readMappings.insert(output->readMappings.end(), l2Mappings.begin(), l2Mappings.end());

          // indicate that we mapped full length
          split_mapping = false;
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

            MappingResultsVector_t l2Mappings;   

            //Map this fragment
            mapSingleQueryFrag(Q, l2Mappings);

            //Adjust query coordinates and length in the reported mapping
            std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){ 
                e.queryLen = input->len;
                e.queryStartPos = i * param.segLength;
                e.queryEndPos = i * param.segLength + Q.len;
                });

            // save the output
            output->readMappings.insert(output->readMappings.end(), l2Mappings.begin(), l2Mappings.end());
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

            MappingResultsVector_t l2Mappings;   

            //Map this fragment
            mapSingleQueryFrag(Q, l2Mappings);

            //Adjust query coordinates and length in the reported mapping
            std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){ 
                e.queryLen = input->len;
                e.queryStartPos = input->len - param.segLength;
                e.queryEndPos = input->len;
                });

            output->readMappings.insert(output->readMappings.end(), l2Mappings.begin(), l2Mappings.end());
          }

          // merge mappings
          if (param.mergeMappings) {
              mergeMappingsInRange(output->readMappings, param.block_length_min);
          }
        }

        // remove self-mode don't-maps
        this->filterSelfingLongToShorts(output->readMappings);

        //filter mappings best over query sequence axis
        if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
            skch::Filter::query::filterMappings(output->readMappings,
                                                (input->len < param.segLength ?
                                                 param.numMappingsForShortSequence
                                                 : param.numMappingsForSegment) - 1);
        }

        // remove short merged mappings when we are merging
        if (split_mapping) {
            // filter mappings that didn't reach the min block length through merging
            filterShortMappings(output->readMappings, param.block_length_min);
            // merge mappings in range and get a copy of the unmerged, annotated with segment id
            auto unmerged = mergeMappingsInRange(output->readMappings, param.chain_gap);
            if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                skch::Filter::query::filterMappings(output->readMappings,
                                                    (input->len < param.segLength ?
                                                     param.numMappingsForShortSequence
                                                     : param.numMappingsForSegment) - 1);
            }
            // trivial: make a set of the kept segment ids
            // use this to filter the unmerged
            robin_hood::unordered_set<offset_t> x;
            // do this w/o the hash table
            for (auto& m : output->readMappings) {
                x.insert(m.splitMappingId);
            }
            filterFailedSubMappings(unmerged, x);
            output->readMappings = unmerged;
        }

        // remove alignments where the ratio between query and target length is < our identity threshold
        this->filterFalseHighIdentity(output->readMappings);

        //Make sure mapping boundary don't exceed sequence lengths
        this->mappingBoundarySanityCheck(input, output->readMappings);

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

          progress.increment(output->qseqLen/2 + (output->qseqLen % 2 != 0));

          delete output;
        }

      /**
       * @brief                   map the parsed query sequence (L1 and L2 mapping)
       * @param[in]   Q           metadata about query sequence
       * @param[in]   outstrm     outstream stream where mappings will be reported
       * @param[out]  l2Mappings  Mapping results in the L2 stage
       */
      template<typename Q_Info, typename VecOut>
        void mapSingleQueryFrag(Q_Info &Q, VecOut &l2Mappings)
        {
#ifdef ENABLE_TIME_PROFILE_L1_L2
          auto t0 = skch::Time::now();
#endif
          //L1 Mapping
          std::vector<L1_candidateLocus_t> l1Mappings; 
          doL1Mapping(Q, l1Mappings);

#ifdef ENABLE_TIME_PROFILE_L1_L2
          std::chrono::duration<double> timeSpentL1 = skch::Time::now() - t0;
          auto t1 = skch::Time::now();
#endif

          //L2 Mapping
          doL2Mapping(Q, l1Mappings, l2Mappings);


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
        void doL1Mapping(Q_Info &Q, Vec &l1Mappings)
        {
          //Vector of positions of all the hits 
          std::vector<MinimizerMetaData> seedHitsL1;

          ///1. Compute the minimizers

          if (param.spaced_seeds.empty()) {
              if (param.world_minimizers) {
                  CommonFunc::addWorldMinimizers(Q.minimizerTableQuery, Q.seq, Q.len, param.kmerSize, param.windowSize, param.alphabetSize, Q.seqCounter);
              } else {
                  CommonFunc::addMinimizers(Q.minimizerTableQuery, Q.seq, Q.len, param.kmerSize, param.windowSize, param.alphabetSize, Q.seqCounter);
              }
          } else {
            CommonFunc::addSpacedSeedMinimizers(Q.minimizerTableQuery, Q.seq, Q.len, param.kmerSize, param.windowSize, param.alphabetSize, Q.seqCounter, param.spaced_seeds);
          }

#ifdef DEBUG
          std::cerr << "[wfmash::skch::Map:doL1Mapping] read id " << Q.seqCounter << ", minimizer count = " << Q.minimizerTableQuery.size() << "\n";
#endif

          ///2. Find the hits in the reference, pick 's' unique minimizers as seeds, 

          std::sort(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::lessByHash);

          //note : unique preserves the original relative order of elements 
          auto uniqEndIter = std::unique(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::equalityByHash);

          //This is the sketch size for estimating jaccard
          Q.sketchSize = std::distance(Q.minimizerTableQuery.begin(), uniqEndIter);

          //For invalid query (example : just NNNs), we may be left with 0 sketch size
          //Ignore the query in this case
          if(Q.sketchSize == 0)
            return;

          int totalMinimizersPicked = 0;

          for(auto it = Q.minimizerTableQuery.begin(); it != uniqEndIter; it++)
          {
            //Check if hash value exists in the reference lookup index
            auto seedFind = refSketch.minimizerPosLookupIndex.find(it->hash);

            if(seedFind != refSketch.minimizerPosLookupIndex.end())
            {
              auto hitPositionList = seedFind->second;

              //Save the positions (Ignore high frequency hits)
              if(hitPositionList.size() < refSketch.getFreqThreshold())
              {
                seedHitsL1.insert(seedHitsL1.end(), hitPositionList.begin(), hitPositionList.end());
              }

            }
          }

          int minimumHits = Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity, skch::fixed::confidence_interval);

          this->computeL1CandidateRegions(Q, seedHitsL1, minimumHits, l1Mappings);

#ifdef DEBUG
          std::cerr << "[wfmash::skch::Map:doL1Mapping] read id " << Q.seqCounter << ", Count of L1 hits in the reference = " << seedHitsL1.size() << ", minimum hits required for a candidate = " << minimumHits << ", Count of L1 candidate regions = " << l1Mappings.size() << "\n";
#endif

        }

      /**
       * @brief                     Helper function to doL1Mapping()
       * @param[in]   Q             query
       * @param[in]   seedHitsL1    minimizer hits in the reference
       * @param[in]   minimumHits   estimated minimum hits required for significant match
       * @param[out]  l1Mappings    all the read mapping locations
       */
      template <typename Q_Info, typename Vec1, typename Vec2>
        void computeL1CandidateRegions(Q_Info &Q, Vec1 &seedHitsL1, int minimumHits, Vec2 &l1Mappings)
        {
          if(minimumHits < 1)
            minimumHits = 1;

          //Sort all the hit positions
          std::sort(seedHitsL1.begin(), seedHitsL1.end());

          for(auto it = seedHitsL1.begin(); it != seedHitsL1.end(); it++)
          {
            if(std::distance(it, seedHitsL1.end()) >= minimumHits)
            {
              auto it2 = it + minimumHits -1;
              //[it .. it2] are 'minimumHits' consecutive hits 

              //Check if consecutive hits are close enough
              //NOTE: hits may span more than a read length for a valid match, as we keep window positions 
              //      for each minimizer
              if(it2->seqId == it->seqId && it2->wpos - it->wpos < Q.len)
              {
                //Save <1st pos --- 2nd pos>
                L1_candidateLocus_t candidate{it->seqId, 
                    std::max(0, it2->wpos - Q.len + 1), it->wpos};

                //Check if this candidate overlaps with last inserted one
                auto lst = l1Mappings.end(); lst--;

                //match seq_no and see if this candidate begins before last element ends
                if( l1Mappings.size() > 0 
                    && candidate.seqId == lst->seqId 
                    && lst->rangeEndPos >= candidate.rangeStartPos)
                {
                  //Push the end pos of last candidate locus further out
                  lst->rangeEndPos = std::max(candidate.rangeEndPos, lst->rangeEndPos);
                }
                else
                  l1Mappings.push_back(candidate);
              }
            }
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
      template <typename Q_Info, typename VecIn, typename VecOut>
        void doL2Mapping(Q_Info &Q, VecIn &l1Mappings, VecOut &l2Mappings)
        {
          ///2. Walk the read over the candidate regions and compute the jaccard similarity with minimum s sketches
          for(auto &candidateLocus: l1Mappings)
          {
            L2_mapLocus_t l2 = {};
            computeL2MappedRegions(Q, candidateLocus, l2);

            //Compute mash distance using calculated jaccard
            float mash_dist = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);

            //Compute lower bound to mash distance within 95% confidence interval
            float mash_dist_lower_bound = Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, skch::fixed::confidence_interval);

            float nucIdentity = (1 - mash_dist);
            float nucIdentityUpperBound = (1 - mash_dist_lower_bound);



            //Report the alignment if it passes our identity threshold and,
            // if we are in all-vs-all mode, it isn't a self-mapping,
            // and if we are self-mapping, the query is shorter than the target
            const auto& ref = this->refSketch.metadata[l2.seqId];
            if((param.keep_low_pct_id && nucIdentityUpperBound >= param.percentageIdentity
                || nucIdentity >= param.percentageIdentity)
               && !(param.skip_self && Q.seqName == ref.name)
               && !(param.skip_prefix
                    && prefix(Q.seqName, param.prefix_delim)
                    == prefix(ref.name, param.prefix_delim)))
            {
              MappingResult res;

              //Save the output
              {
                res.queryLen = Q.len;
                res.refStartPos = l2.meanOptimalPos ;
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

                res.selfMapFilter = ((param.skip_self || param.skip_prefix) && Q.fullLen > ref.len);

                //Compute additional statistics -> strand, reference complexity
                {
                  SlideMapper<Q_Info> slidemap(Q);
                  slidemap.insert_ref(l2.optimalStart, l2.optimalEnd);
                  int strandVotes, uniqueRefHashes;
                  slidemap.computeStatistics(strandVotes, uniqueRefHashes);

                  res.strand = strandVotes > 0 ? strnd::FWD : strnd::REV;
                }

                l2Mappings.push_back(res);
              }
            }
          }

#ifdef DEBUG
          std::cerr << "[wfmash::skch::Map:doL2Mapping] read id " << Q.seqCounter << ", count of L1 candidates= " << l1Mappings.size() << ", count of L2 candidates= " << l2Mappings.size() << std::endl;
#endif
        }

      /**
       * @brief                                 Find optimal mapping within an L1 candidate
       * @param[in]   Q                         query sequence information
       * @param[in]   candidateLocus            L1 candidate location
       * @param[out]  l2_out                    L2 mapping inside L1 candidate 
       */
      template <typename Q_Info>
        void computeL2MappedRegions(Q_Info &Q, 
            L1_candidateLocus_t &candidateLocus, 
            L2_mapLocus_t &l2_out)
        {
          //Look up L1 candidate's begin in the index
          MIIter_t firstSuperWindowRangeStart = this->refSketch.searchIndex(candidateLocus.seqId, 
              candidateLocus.rangeStartPos);

          //If iterator points to index end or a different reference sequence, return as there are no candidates 
          if ( this->refSketch.isMinimizerIndexEnd(firstSuperWindowRangeStart) || firstSuperWindowRangeStart->seqId != candidateLocus.seqId)
            return;

          //Count of minimizer windows in a super-window
          offset_t countMinimizerWindows = Q.len - (param.windowSize-1) - (param.kmerSize-1); 

          //Look up the end of the first L2 super-window in the index
          MIIter_t firstSuperWindowRangeEnd = this->refSketch.searchIndex(candidateLocus.seqId, 
              firstSuperWindowRangeStart->wpos + countMinimizerWindows);

          //Look up L1 candidate's end in the index
          MIIter_t lastSuperWindowRangeEnd = this->refSketch.searchIndex(candidateLocus.seqId, 
              candidateLocus.rangeEndPos + Q.len);

          //Define std::map such that it contains only the query minimizers
          //Used to efficiently compute the jaccard similarity between qry and ref
          SlideMapper<Q_Info> slidemap(Q);

          //Initialize iterator over minimizerIndex
          MIIteratorL2 mi_L2iter( firstSuperWindowRangeStart, firstSuperWindowRangeEnd,
              countMinimizerWindows);

          //Insert all the minimizers in the first 'super-window'
          //  [ mi_L2iter.sw_beg, mi_L2iter.sw_end )
          slidemap.insert_ref(mi_L2iter.sw_beg, mi_L2iter.sw_end);

          auto prev_beg_iter = mi_L2iter.sw_beg;
          auto prev_end_iter = mi_L2iter.sw_end;

          int beginOptimalPos, lastOptimalPos;

          while ( std::distance(mi_L2iter.sw_end, lastSuperWindowRangeEnd) >= 0)
          {
            assert( std::distance(mi_L2iter.sw_beg, firstSuperWindowRangeStart) <= 0);
            assert( std::distance(mi_L2iter.sw_end, lastSuperWindowRangeEnd  ) >= 0);

            //Check if the previous first minimizer is out of current range
            if (prev_beg_iter != mi_L2iter.sw_beg)
              slidemap.delete_ref(prev_beg_iter);

            //Check if we have new minimizer in the current range
            if (prev_end_iter != mi_L2iter.sw_end)
              slidemap.insert_ref(prev_end_iter);
          
            //Is this sliding window the best we have so far?
            if (slidemap.sharedSketchElements > l2_out.sharedSketchSize)
            {
              l2_out.sharedSketchSize = slidemap.sharedSketchElements;
              l2_out.optimalStart = mi_L2iter.sw_beg;
              l2_out.optimalEnd = mi_L2iter.sw_end;

              //Save the position
              beginOptimalPos = mi_L2iter.sw_beg->wpos;
              lastOptimalPos = mi_L2iter.sw_beg->wpos;
            }
            else if(slidemap.sharedSketchElements == l2_out.sharedSketchSize)
            {
              //Still save the position
              lastOptimalPos = mi_L2iter.sw_beg->wpos; 
            }

            //Back up the current iterator values
            prev_beg_iter = mi_L2iter.sw_beg;
            prev_end_iter = mi_L2iter.sw_end;

            //If current sliding window touches end, we should stop further evaluation
            if(std::distance(mi_L2iter.sw_end, lastSuperWindowRangeEnd) == 0)
              break;

            //Advance the current super-window
            mi_L2iter.next();

          }//End of while loop

          //Save reference sequence id in the mapping output 
          l2_out.seqId = candidateLocus.seqId;
          l2_out.meanOptimalPos = (beginOptimalPos + lastOptimalPos)/2;
        }

      /**
       * @brief                                 Compute jaccard similarity of given query seq 
       *                                        at a single position on reference sequence
       * @param[in]   Q                         query sequence information
       * @param[in]   seqId                     reference sequence id
       * @param[in]   refStartPos               offset on reference sequence
       * return                                 jaccard similarity
       */
      template <typename Q_Info>
        double computeJaccardSinglePos(Q_Info &Q, seqno_t seqId, offset_t refStartPos)
        {
          //Look up L1 candidate's begin in the index
          MIIter_t superWindowRangeStart = this->refSketch.searchIndex(seqId, refStartPos);

          if (Q.sketchSize == 0)
            return 0;
          
          //If iterator points to index end or a different reference sequence, return zero 
          if ( this->refSketch.isMinimizerIndexEnd(superWindowRangeStart) || superWindowRangeStart->seqId != seqId)
            return 0;

          //Count of minimizer windows in a super-window
          offset_t countMinimizerWindows = (Q.len - (param.windowSize-1) - (param.kmerSize-1))
              * (param.spaced_seeds.empty() ? 1 : param.spaced_seeds.size());

          //Look up the end of the first L2 super-window in the index
          MIIter_t superWindowRangeEnd = this->refSketch.searchIndex(seqId, 
              superWindowRangeStart->wpos + countMinimizerWindows);

          //Define std::map and let it contain only the query minimizers
          SlideMapper<Q_Info> slidemap(Q);

          //Insert all the minimizers in the first super-window
          slidemap.insert_ref(superWindowRangeStart, superWindowRangeEnd);

          return slidemap.sharedSketchElements * 1.0 / Q.sketchSize;
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
              if(it2->refSeqId != it->refSeqId || it2->refStartPos - it->refEndPos > 2 * param.segLength)
                break;

              //If the next mapping is within range, check if it is consecutive query fragment and strand matches
              if( it2->strand == it->strand
                  && thisMappingFragno == currMappingFragno + (it->strand == strnd::FWD ? 1 : -1) )
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

            //Mean identity of all mappings in the chain
            it->nucIdentity = (   std::accumulate(it, it_end, 0.0, 
                                  [](double x, MappingResult &e){ return x + e.nucIdentity; })     )/ std::distance(it, it_end);

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
      VecIn mergeMappingsInRange(VecIn &readMappings, int max_dist) {
          assert(param.split == true);

          if(readMappings.size() < 2)
              return readMappings;

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
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              std::vector<std::pair<uint64_t, uint64_t>> distances;
              for (auto it2 = std::next(it); it2 != readMappings.end(); it2++) {
                  //If this mapping is too far from current mapping being evaluated, stop finding a merge
                  if (it2->refSeqId != it->refSeqId || it2->refStartPos > it->refEndPos + max_dist) {
                      break;
                  }

                  //If the next mapping is within range, check if it's in range and 
                  if (it2->strand == it->strand) {
                      auto ref_dist = it2->refStartPos - it->refEndPos;
                      auto score = std::numeric_limits<uint64_t>::max();
                      bool ok = false;
                      if (it->strand == strnd::FWD && it->queryStartPos < it2->queryStartPos) {
                          auto query_dist = it2->queryStartPos - it->queryEndPos;
                          ok = query_dist + ref_dist < max_dist;
                          score = ref_dist + query_dist;
                      } else if (it->strand != strnd::FWD && it->queryEndPos > it2->queryEndPos) {
                          auto query_dist = it->queryStartPos - it2->queryEndPos;
                          ok = query_dist + ref_dist < max_dist;
                          score = ref_dist + query_dist;
                      }
                      if (ok) distances.push_back(std::make_pair(score, it2->splitMappingId));
                      //if (ok) disjoint_sets.unite(it->splitMappingId, it2->splitMappingId);
                  }
              }
              if (distances.size()) {
                  std::sort(distances.begin(), distances.end());
                  disjoint_sets.unite(it->splitMappingId, distances.front().second);
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

          // copy for independent return of merged and unmerged
          auto unmergedReadMappings = readMappings;

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

              //Mean identity of all mappings in the chain
              it->nucIdentity = (   std::accumulate(it, it_end, 0.0, 
                                                    [](double x, MappingResult &e){ return x + e.nucIdentity; })     )/ std::distance(it, it_end);

              //Discard other mappings of this chain
              std::for_each( std::next(it), it_end, [&](MappingResult &e){ e.discard = 1; });

              //advance the iterator
              it = it_end;
          }

          readMappings.erase( 
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());

          return unmergedReadMappings;
      }

     /**
       * @brief                       This routine is to make sure that all mapping boundaries
       *                              on query and reference are not outside total 
       *                              length of sequeunces involved
       * @param[in]     input         input read details
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void mappingBoundarySanityCheck(InputSeqContainer* input, VecIn &readMappings)
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

          float fakeMapQ = std::round(-10.0 * std::log10(1-(e.nucIdentity)));
          if (std::isinf(fakeMapQ)) fakeMapQ = 255;

          outstrm  << (param.filterMode == filter::ONETOONE ? qmetadata[e.querySeqId].name : queryName)
                   << "\t" << e.queryLen 
                   << "\t" << e.queryStartPos
                   << "\t" << e.queryEndPos
                   << "\t" << (e.strand == strnd::FWD ? "+" : "-") 
                   << "\t" << this->refSketch.metadata[e.refSeqId].name
                   << "\t" << this->refSketch.metadata[e.refSeqId].len
                   << "\t" << e.refStartPos 
                   << "\t" << e.refEndPos
                   << "\t" << e.approxMatches
                   << "\t" << e.blockLength
                   << "\t" << fakeMapQ
                   << "\t" << "id:f:" << e.nucIdentity * 100.0;
              //<< "\t" << "nu:f:" << e.nucIdentityUpperBound;

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
