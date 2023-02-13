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
#include <iostream>

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
#include "map_stats.hpp"
#include "robin-hood-hashing/robin_hood.h"
// if we ever want to do the union-find chaining in parallel
//#include "common/dset64-gccAtomic.hpp"
// this is for single-threaded use, but is more portable
#include "common/dset64.hpp"
#include "assert.hpp"
#include "gsl/gsl_randist.h"

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
        int intersectionSize;
        strand_t strand;
      };

      //Type for Stage L2's predicted mapping coordinate within each L1 candidate
      struct L2_mapLocus_t
      {
        seqno_t seqId;                    //sequence id where read is mapped
        offset_t meanOptimalPos;          //Among multiple consecutive optimal positions, save the avg.
        //Sketch::MIIter_t optimalStart;    //optimal start mapping position (begin iterator)
        //Sketch::MIIter_t optimalEnd;      //optimal end mapping position (end iterator)
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
      std::vector<std::vector<int>> sketchCutoffs; 

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
        sketchCutoffs(p.sketchSize + 1, std::vector<int>(p.sketchSize + 1, 0)),
        refSketch(refsketch),
        processMappingResults(f)
    {
      this->setProbs();
      this->mapQuery();
    }

    private:

      void setProbs() 
      {
        std::vector<std::vector<std::vector<double>>> sketchProbs(
            param.sketchSize + 1,
            std::vector<std::vector<double>>(
              param.sketchSize + 1, 
              std::vector<double>( param.sketchSize + 1, 0))
        );
        
        for (int ss = 1; ss <= param.sketchSize; ss++) 
        {
          for (auto ci = 0; ci <= param.sketchSize; ci++) 
          {
            for (double y = 0; y <= ci; y++) {
              sketchProbs[ss][ci][y] = gsl_ran_hypergeometric_pdf(y, ss, ss-ci, ci);
            }
          }
          

          float deltaANI = param.ANIDiff;
          float min_p = 1 - param.ANIDiffConf;

          for (auto cmax = 1; cmax <= param.sketchSize; cmax++) 
          {
            for (auto ci = 1; ci <= cmax; ci++) 
            {
              double prAboveCutoff = 0;
              for (double ymax = 0; ymax <= cmax; ymax++) {
                // Pr (Ymax = ymax)
                double pymax = sketchProbs[ss][cmax][ymax];

                // yi_cutoff is minimum jaccard numerator required to be within deltaANI of ymax
                double yi_cutoff = deltaANI == 0 ? ymax : (std::floor(skch::Stat::md2j(
                    skch::Stat::j2md(ymax / param.sketchSize, param.kmerSize) + deltaANI, 
                    param.kmerSize
                ) * param.sketchSize));

                // Pr Y_i < yi_cutoff
                double pi_acc = (yi_cutoff - 1) >= 0 ? gsl_cdf_hypergeometric_P(yi_cutoff-1, ss, ss-ci, ci) : 0;
                // Pr Y_i >= yi_cutoff
                pi_acc = 1-pi_acc;

                // Pr that mash score from cj leads to an ANI at least deltaJ less than the ANI from cmax
                prAboveCutoff += pymax * pi_acc;
              }
              if (prAboveCutoff > min_p) 
              {
                sketchCutoffs[ss][cmax] = std::max(1, ci);
                break;
              }
            }
            // For really high min_p values and some values of cmax, there are no values of
            // ci that satisfy the cutoff, so we just set to the max
            if (sketchCutoffs[ss][cmax] == 0) {
              sketchCutoffs[ss][cmax] = cmax;
            }
          }
          for (auto overlap = 0; overlap <= ss; overlap++) {
            DEBUG_ASSERT(sketchCutoffs[ss][overlap] <= overlap);
          }
        }
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
                std::cerr << "[mashmap::skch::Map::mapQuery] WARNING, no .fai index found for " << fileName << ", reading the file to sum sequence length (slow)" << std::endl;
                seqiter::for_each_seq_in_file(
                    fileName,
                    [&](const std::string& seq_name,
                        const std::string& seq) {
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

            seqiter::for_each_seq_in_file(
                fileName,
                [&](const std::string& seq_name,
                    const std::string& seq) {
                    // todo: offset_t is an 32-bit integer, which could cause problems
                    offset_t len = seq.length();

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

        if(! param.split || input->len <= param.block_length)
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
        }

        // how many mappings to keep
        int n_mappings = (input->len < param.segLength ?
                          param.numMappingsForShortSequence
                          : param.numMappingsForSegment) - 1;

        if (split_mapping) {
            if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                skch::Filter::query::filterMappings(output->readMappings, n_mappings);
            }
            if (param.mergeMappings) {
              // hardcore merge using the chain gap
              mergeMappingsInRange(output->readMappings, param.chain_gap);
              if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                  skch::Filter::query::filterMappings(output->readMappings, n_mappings);
              }
              // remove short chains that didn't exceed block length
              filterWeakMappings(output->readMappings, std::floor(param.block_length / param.segLength));
            }
        } else {
            // filter the non-split mappings using plane sweep
            if (param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE) {
                skch::Filter::query::filterMappings(output->readMappings, n_mappings);
            }
        }

        // remove self-mode don't-maps
        this->filterSelfingLongToShorts(output->readMappings);

        // remove alignments where the ratio between query and target length is < our identity threshold
        this->filterFalseHighIdentity(output->readMappings);

        //Make sure mapping boundary don't exceed sequence lengths
        this->mappingBoundarySanityCheck(input, output->readMappings);

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
          if (l1Mappings.size() == 0) {
            return;
          }

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

      template <typename Q_Info>
        void getSeedHits(Q_Info &Q)
        {
          CommonFunc::sketchSequence(Q.minmerTableQuery, Q.seq, Q.len, param.kmerSize, param.alphabetSize, param.sketchSize, Q.seqCounter);
          if(Q.minmerTableQuery.size() == 0) {
            Q.sketchSize = 0;
            return;
          }

#ifdef DEBUG
          int orig_len = Q.minmerTableQuery.size();
#endif
          // TODO remove them from the original sketch instead of removing for each read
          auto new_end = std::remove_if(Q.minmerTableQuery.begin(), Q.minmerTableQuery.end(), [&](auto mi) {
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

          for(auto it = Q.minmerTableQuery.begin(); it != Q.minmerTableQuery.end(); it++)
          {
            //Check if hash value exists in the reference lookup index
            auto seedFind = refSketch.minmerPosLookupIndex.find(it->hash);

            if(seedFind != refSketch.minmerPosLookupIndex.end())
            {
              auto& hitPositionList = seedFind->second;

              // Let the strand of the hits denote wrt the reference. 
              // i.e. if - query mashimizer hits a - ref mashimizer, we mark the strand as fwd. 
              std::for_each(hitPositionList.begin(), hitPositionList.end(), [&](auto& mi) {
                  //mi.strand = mi.strand * it->strand;
                  // Only add hits w/ same name if !param.skip_self
                  if (!param.skip_self || Q.seqName != this->refSketch.metadata[mi.seqId].name) {
                    // Only add both points if not overlapping
                    // Since we always add the close interval second, we know we can check the last intervalpoint
                    if (intervalPoints.size() == 0 
                        || intervalPoints.back().hash != mi.hash 
                        || intervalPoints.back().pos != mi.wpos)
                    {
                      intervalPoints.push_back(IntervalPoint {side::OPEN, mi.seqId, mi.wpos, strand_t(mi.strand * it->strand), mi.hash});
                      intervalPoints.push_back(IntervalPoint {side::CLOSE, mi.seqId, mi.wpos_end, strand_t(mi.strand * it->strand), mi.hash});
                    } else {
                      intervalPoints.back().pos = mi.wpos_end;
                    }
                  }
              });
            }
          }

          //Sort all the hit positions
          std::sort(intervalPoints.begin(), intervalPoints.end());

#ifdef DEBUG
          std::cerr << "INFO, skch::Map:getSeedHits, read id " << Q.seqCounter << ", Count of seed hits in the reference = " << intervalPoints.size() / 2 << "\n";
#endif
        }


      template <typename Q_Info, typename Vec1, typename Vec2>
        void computeL1CandidateRegions(Q_Info &Q, Vec1 intervalPoints, int minimumHits, Vec2 &l1Mappings)
        {
#ifdef DEBUG
          std::cerr << "INFO, skch::Map:computeL1CandidateRegions, read id " << Q.seqCounter << std::endl;
#endif

          int overlapCount = 0;
          int strandCount = 0;
          int bestIntersectionSize = 0;
          std::vector<L1_candidateLocus_t> localOpts;

          // Keep track of all minmer windows that intersect with [i, i+windowLen]
          int windowLen = Q.len - param.segLength;
          auto trailingIt = intervalPoints.begin();
          auto leadingIt = intervalPoints.begin();

          // Group together local sketch intersection maximums that are within clusterLen of eachother
          //
          // Since setting up the L2 window [i, j] requires aggregating minmer windows over
          // [i-segLength, i), we might as well group L2 windows together which are closer than
          // segLength  
          int clusterLen = param.segLength;

          // Used to keep track of how many minmer windows for a particular hash are currently "open"
          // Only necessary when windowLen != 0.
          std::unordered_map<hash_t, int> hash_to_freq;
          
          if (param.stage1_topANI_filter) {
            while (leadingIt != intervalPoints.end())
            {
              // Catch the trailing iterator up to the leading iterator - windowLen
              while (
                  trailingIt != intervalPoints.end() 
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
              while (leadingIt != intervalPoints.end() && leadingIt->pos == currentPos) {
                if (leadingIt->side == side::OPEN) {
                  if (windowLen == 0 || hash_to_freq[leadingIt->hash] == 0) {
                    overlapCount++;
                  }
                  if (windowLen != 0)
                    hash_to_freq[leadingIt->hash]++;
                }
                leadingIt++;
              }

              DEBUG_ASSERT(overlapCount >= 0, windowLen, trailingIt->seqId, trailingIt->pos, leadingIt->seqId, leadingIt->pos);
              DEBUG_ASSERT(windowLen != 0 || overlapCount <= Q.sketchSize, windowLen, trailingIt->seqId, trailingIt->pos, leadingIt->seqId, leadingIt->pos);

              //Is this sliding window the best we have so far?
              bestIntersectionSize = std::max(bestIntersectionSize, overlapCount);
            }

            // Only go back through to find local opts if we know that there are some that are 
            // large enough
            if (bestIntersectionSize < minimumHits) 
            {
              return;
            }
          } 

          // Clear freq dict, as there will be left open CLOSE points at the end of the last seq
          // that we never got to
          hash_to_freq.clear();

          // Since there can be more than sketchSize windows that overlap w/ [i, i+windowLen]
          // cap the best intersection size 
          bestIntersectionSize = std::min(bestIntersectionSize, Q.sketchSize);

          bool in_candidate = false;
          L1_candidateLocus_t l1_out = {};
          trailingIt = intervalPoints.begin();
          leadingIt = intervalPoints.begin();

          // Keep track of 3 consecutive points so that we can track local optimums
          overlapCount = 0;
          int prevOverlap = 0;
          int prevPrevOverlap = 0;

          // Need to keep track of two positions, as the previous one will be the local optimum
          SeqCoord prevPos, currentPos;


          while (leadingIt != intervalPoints.end())
          {
            prevPrevOverlap = prevOverlap;
            prevOverlap = overlapCount;

            //TODO LEADING it should only hit opens
            // We should only iterate through a new window when we come across an OPEN,
            // right now, this basically happens since every CLOSE should be an OPEN.
            // This doesn't invalidate the logic, just potentially wastes time
            while (
                trailingIt != intervalPoints.end() 
                && ((trailingIt->seqId == leadingIt->seqId && trailingIt->pos <= leadingIt->pos - windowLen)
                  || trailingIt->seqId < leadingIt->seqId))
            {
              if (trailingIt->side == side::CLOSE) {
                if (windowLen != 0)
                  hash_to_freq[trailingIt->hash]--;
                if (windowLen == 0 || hash_to_freq[trailingIt->hash] == 0) {
                  overlapCount--;
                  strandCount -= trailingIt->strand;
                }
              }
              trailingIt++;
            }
            if (leadingIt->pos != currentPos.pos) {
              prevPos = currentPos;
              currentPos = SeqCoord{leadingIt->seqId, leadingIt->pos};
            }
            while (leadingIt != intervalPoints.end() && leadingIt->pos == currentPos.pos) 
            {
              if (leadingIt->side == side::OPEN) {
                if (windowLen == 0 || hash_to_freq[leadingIt->hash] == 0) {
                  overlapCount++;
                  strandCount += trailingIt->strand;
                }
                if (windowLen != 0)
                  hash_to_freq[leadingIt->hash]++;
              }
              leadingIt++;
            }
          if ((!param.stage1_topANI_filter && prevOverlap >= minimumHits)
              || prevOverlap >= bestIntersectionSize 
              || (prevOverlap >= sketchCutoffs[Q.sketchSize][bestIntersectionSize] 
              //&& prevOverlap > overlapCount && prevOverlap >= prevPrevOverlap)
             )
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
              l1_out.strand = (strandCount >= 0) ? strnd::FWD : strnd::REV;
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
                l1_out.strand = (strandCount >= 0) ? strnd::FWD : strnd::REV;
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
      template <typename Q_Info, typename Vec>
        void doL1Mapping(Q_Info &Q, Vec &l1Mappings)
        {
          //1. Compute the minmers
          getSeedHits(Q);

          //Catch all NNNNNN case
          if (Q.sketchSize == 0) {
            return;
          }

          //2. Compute windows and sort
          std::vector<skch::IntervalPoint> intervalPoints; 
          getSeedIntervalPoints(Q, intervalPoints);

          //3. Compute L1 windows
          int minimumHits = Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity, skch::fixed::confidence_interval);
          computeL1CandidateRegions(Q, intervalPoints, minimumHits, l1Mappings);
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
            std::vector<L2_mapLocus_t> l2_vec;
            if (true) {
              computeL2MappedRegions(Q, candidateLocus, l2_vec);
            } else {
              l2_vec.push_back(
                L2_mapLocus_t {
                  candidateLocus.seqId,
                  candidateLocus.rangeStartPos + ((int)(candidateLocus.rangeStartPos - candidateLocus.rangeEndPos) / 2),
                  (int) (((float)(2*param.sketchSize - candidateLocus.intersectionSize)) / param.sketchSize) * candidateLocus.intersectionSize,
                  candidateLocus.strand
                }
              );
            }
            for (auto& l2 : l2_vec) 
            {
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

                  res.selfMapFilter = ((param.skip_self || param.skip_prefix) && Q.fullLen > ref.len);

                  l2Mappings.push_back(res);
                } 
              }
            }
          }

#ifdef DEBUG
          std::cerr << "[mashmap::skch::Map:doL2Mapping] read id " << Q.seqCounter << ", count of L1 candidates= " << l1Mappings.size() << ", count of L2 candidates= " << l2Mappings.size() << std::endl;
#endif
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
          
          // Get first potential mashimizer
          const MinmerInfo first_minmer = MinmerInfo {0, candidateLocus.seqId, candidateLocus.rangeStartPos - param.segLength - 1, 0, 0};

          //const MinmerInfo first_minmer = MinmerInfo {0, candidateLocus.seqId, -1, 0, 0};
          auto firstOpenIt = std::lower_bound(minmerIndex.begin(), minmerIndex.end(), first_minmer); 

          // Keeps track of the lowest end position
          std::vector<skch::MinmerInfo> slidingWindow;

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


            //Is this sliding window the best we have so far?
            if (slideMap.sharedSketchElements > bestSketchSize)
            {
              in_candidate = true;
              bestSketchSize = slideMap.sharedSketchElements;
              l2_out.sharedSketchSize = slideMap.sharedSketchElements;

              //Save the position
              beginOptimalPos = windowIt->wpos - windowLen;
              lastOptimalPos = std::next(windowIt, windowIt != minmerIndex.end())->wpos - windowLen;
            }
            else if(slideMap.sharedSketchElements == bestSketchSize)
            {
              if (!in_candidate) {
                l2_out.sharedSketchSize = slideMap.sharedSketchElements;

                //Save the position
                beginOptimalPos = windowIt->wpos - windowLen;
              }

              in_candidate = true;
              //Still save the position
              lastOptimalPos = windowIt->wpos - windowLen;
            } else {
              if (in_candidate) {
                // Save and reset
                lastOptimalPos = windowIt->wpos - windowLen;
                l2_out.meanOptimalPos =  (beginOptimalPos + lastOptimalPos) / 2;
                l2_out.seqId = windowIt->seqId;
                l2_out.strand = prev_strand_votes >= 0 ? strnd::FWD : strnd::REV;
                l2_vec_out.push_back(l2_out);
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
            l2_out.meanOptimalPos =  (beginOptimalPos + lastOptimalPos) / 2;
            l2_out.seqId = std::prev(windowIt)->seqId;
            l2_out.strand = slideMap.strand_votes >= 0 ? strnd::FWD : strnd::REV;
            l2_vec_out.push_back(l2_out);
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
          for (auto it = readMappings.begin(); it != readMappings.end(); it++) {
              std::vector<std::pair<double, uint64_t>> distances;
              for (auto it2 = std::next(it); it2 != readMappings.end(); it2++) {
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
                  ) /// it->n_merged; // this would scale directly by the number of mappings in the chain
                  // this scales slightly by the amount of missing segments
                  / ( (double)it->n_merged
                      + std::pow(
                          std::log((double)it->blockLength / param.segLength),
                          0.01));


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
                   << "\t" << e.conservedSketches
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
