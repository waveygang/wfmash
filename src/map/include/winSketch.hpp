/**
 * @file    winSketch.hpp
 * @brief   routines to index the reference 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WIN_SKETCH_HPP 
#define WIN_SKETCH_HPP

#include <algorithm>
#include <cassert>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;

//#include <zlib.h>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/ThreadPool.hpp"

//External includes
#include "common/murmur3.h"
#include "common/prettyprint.hpp"
#include "csv.h"

//#include "common/sparsehash/dense_hash_map"
//#include "common/parallel-hashmap/parallel_hashmap/phmap.h"
//#include <abseil-cpp/absl/container/flat_hash_map.h>
//#include <common/sparse-map/include/tsl/sparse_map.h>
//#include <common/robin-hood-hashing/robin_hood.h>
#include "common/ankerl/unordered_dense.hpp"

#include "common/seqiter.hpp"

//#include "assert.hpp"

namespace skch
{
  /**
   * @class     skch::Sketch
   * @brief     sketches and indexes the reference (subject sequence)
   * @details  
   *            1.  Minmers are computed in streaming fashion
   *                Computing minmers is using double ended queue which gives
   *                O(reference size) complexity
   *                Algorithm described here:
   *                https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
   *
   *            2.  Index hashes into appropriate format to enable fast search at L1 mapping stage
   */
  class Sketch
    {
      //private members
    
      //algorithm parameters
      const skch::Parameters &param;

      //Minmers that occur this or more times will be ignored (computed based on percentageThreshold)
      uint64_t freqThreshold = std::numeric_limits<uint64_t>::max();

      //Set of frequent seeds to be ignored
      ankerl::unordered_dense::set<hash_t> frequentSeeds;

      //Make the default constructor private, non-accessible
      Sketch();

      public:

      using MI_Type = std::vector< MinmerInfo >;
      using MIIter_t = MI_Type::const_iterator;
      using HF_Map_t = ankerl::unordered_dense::map<hash_t, uint64_t>;

      // Frequency of each hash
      HF_Map_t hashFreq;

      //Keep sequence length, name that appear in the sequence (for printing the mappings later)
      std::vector< ContigInfo > metadata;

      /*
       * Keep the information of what sequences come from what file#
       * Example [a, b, c] implies 
       *  file 0 contains 0 .. a-1 sequences
       *  file 1 contains a .. b-1 
       *  file 2 contains b .. c-1
       */
      std::vector< int > sequencesByFileInfo;

      //Index for fast seed lookup (unordered_map)
      /*
       * [minmer #1] -> [pos1, pos2, pos3 ...]
       * [minmer #2] -> [pos1, pos2...]
       * ...
       */
      //using MI_Map_t = google::dense_hash_map< MinmerMapKeyType, MinmerMapValueType >;
      //using MI_Map_t = phmap::flat_hash_map< MinmerMapKeyType, MinmerMapValueType >;
      //using MI_Map_t = absl::flat_hash_map< MinmerMapKeyType, MinmerMapValueType >;
      //using MI_Map_t = tsl::sparse_map< MinmerMapKeyType, MinmerMapValueType >;
      using MI_Map_t = ankerl::unordered_dense::map< MinmerMapKeyType, MinmerMapValueType >;
      MI_Map_t minmerPosLookupIndex;
      MI_Type minmerIndex;

      private:

      /**
       * Keep list of minmers, sequence# , their position within seq , here while parsing sequence 
       * Note : position is local within each contig
       * Hashes saved here are non-unique, ordered as they appear in the reference
       */

      //Frequency histogram of minmers
      //[... ,x -> y, ...] implies y number of minmers occur x times
      std::map<uint64_t, uint64_t> minmerFreqHistogram;

      public:

      /**
       * @brief   constructor
       *          also builds, indexes the minmer table
       */
      Sketch(const skch::Parameters &p) 
        :
          param(p) {
            if (param.indexFilename.empty() 
                || !stdfs::exists(param.indexFilename)
                || param.overwrite_index)
            {
              this->build(true);
              this->computeFreqHist();
              this->computeFreqSeedSet();
              this->dropFreqSeedSet();
              this->hashFreq.clear();
              if (!param.indexFilename.empty())
              {
                this->writeIndex();
              }
            } else {
              this->build(false);
              this->readIndex();
            }
            std::cerr << "[mashmap::skch::Sketch] Unique minmer hashes after pruning = " << minmerPosLookupIndex.size() << std::endl;
            std::cerr << "[mashmap::skch::Sketch] Total minmer windows after pruning = " << minmerIndex.size() << std::endl;
          }

      private:

      /**
       * @brief     Get sequence metadata and optionally build the sketch table
       *
       * @details   Iterate through ref sequences to get metadata and
       *            optionally compute and save minmers from the reference sequence(s)
       *            assuming a fixed window size
       */
      void build(bool compute_seeds)
      {

        // allowed set of targets
        std::unordered_set<std::string> allowed_target_names;
        if (!param.target_list.empty()) {
                std::ifstream filter_list(param.target_list);
                std::string name;
                while (getline(filter_list, name)) {
                        allowed_target_names.insert(name); 
                }
        }
		

        //sequence counter while parsing file
        seqno_t seqCounter = 0;

        //Create the thread pool 
        ThreadPool<InputSeqContainer, MI_Type> threadPool( [this](InputSeqContainer* e) {return buildHelper(e);}, param.threads);

        for(const auto &fileName : param.refSequences)
        {

#ifdef DEBUG
        std::cerr << "[mashmap::skch::Sketch::build] building minmer index for " << fileName << std::endl;
#endif

        seqiter::for_each_seq_in_file(
            fileName,
            allowed_target_names,
            param.target_prefix,
            [&](const std::string& seq_name,
                const std::string& seq) {
                // todo: offset_t is an 32-bit integer, which could cause problems
                offset_t len = seq.length();

                //Save the sequence name
                metadata.push_back( ContigInfo{seq_name, len} );

                //Is the sequence too short?
                if(len < param.kmerSize)
                {
#ifdef DEBUG
                    std::cerr << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer" << std::endl;
#endif
                }
                else
                {
                  if (compute_seeds) {
                    threadPool.runWhenThreadAvailable(new InputSeqContainer(seq, seq_name, seqCounter));
                    
                    //Collect output if available
                    while ( threadPool.outputAvailable() )
                        this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
                  }
                }
                seqCounter++;
            });

          sequencesByFileInfo.push_back(seqCounter);
        }

        if (seqCounter == 0)
        {
          std::cerr << "[mashmap::skch::Sketch::build] ERROR: No sequences indexed!" << std::endl;
          exit(1);
        }

        if (compute_seeds) {
          //Collect remaining output objects
          while ( threadPool.running() )
            this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
          std::cerr << "[mashmap::skch::Sketch::build] Unique minmer hashes before pruning = " << minmerPosLookupIndex.size() << std::endl;
          std::cerr << "[mashmap::skch::Sketch::build] Total minmer windows before pruning = " << minmerIndex.size() << std::endl;
        }
      }

      /**
       * @brief               function to compute minmers given input sequence object
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @return              output object containing the mappings
       */
      MI_Type* buildHelper(InputSeqContainer *input)
      {
        MI_Type* thread_output = new MI_Type();

        //Compute minmers in reference sequence
        skch::CommonFunc::addMinmers(
                *thread_output, 
                &(input->seq[0u]), 
                input->len, 
                param.kmerSize, 
                param.segLength, 
                param.alphabetSize, 
                param.sketchSize,
                input->seqCounter);

        return thread_output;
      }

      /**
       * @brief                 routine to handle thread's local minmer index
       * @param[in] output      thread local minmer output
       */
      void buildHandleThreadOutput(MI_Type* contigMinmerIndex)
      {
        for (MinmerInfo& mi : *contigMinmerIndex)
        {
          this->hashFreq[mi.hash]++;
          if (minmerPosLookupIndex[mi.hash].size() == 0 
                  || minmerPosLookupIndex[mi.hash].back().hash != mi.hash 
                  || minmerPosLookupIndex[mi.hash].back().pos != mi.wpos)
            {
              minmerPosLookupIndex[mi.hash].push_back(IntervalPoint {mi.wpos, mi.hash, mi.seqId, side::OPEN});
              minmerPosLookupIndex[mi.hash].push_back(IntervalPoint {mi.wpos_end, mi.hash, mi.seqId, side::CLOSE});
            } else {
              minmerPosLookupIndex[mi.hash].back().pos = mi.wpos_end;
          }
        }

        this->minmerIndex.insert(
            this->minmerIndex.end(), 
            std::make_move_iterator(contigMinmerIndex->begin()), 
            std::make_move_iterator(contigMinmerIndex->end()));

        delete contigMinmerIndex;
      }


      /**
       * @brief  Write sketch as tsv. TSV indexing is slower but can be debugged easier
       */
      void writeSketchTSV() 
      {
        std::ofstream outStream;
        outStream.open(std::string(param.indexFilename) + ".tsv");
        outStream << "seqId" << "\t" << "strand" << "\t" << "start" << "\t" << "end" << "\t" << "hash\n";
        for (auto& mi : this->minmerIndex) {
          outStream << mi.seqId << "\t" << std::to_string(mi.strand) << "\t" << mi.wpos << "\t" << mi.wpos_end << "\t" << mi.hash << "\n";
        }
        outStream.close(); 
      }


      /**
       * @brief  Write sketch for quick loading
       */
      void writeSketchBinary(std::ofstream& outStream) 
      {
        typename MI_Type::size_type size = minmerIndex.size();
        outStream.write((char*)&size, sizeof(size));
        outStream.write((char*)&minmerIndex[0], minmerIndex.size() * sizeof(MinmerInfo));
      }

      /**
       * @brief  Write posList for quick loading
       */
      void writePosListBinary(std::ofstream& outStream) 
      {
        typename MI_Map_t::size_type size = minmerPosLookupIndex.size();
        outStream.write((char*)&size, sizeof(size));

        for (auto& [hash, ipVec] : minmerPosLookupIndex) 
        {
          MinmerMapKeyType key = hash;
          outStream.write((char*)&key, sizeof(key));
          typename MI_Type::size_type size = ipVec.size();
          outStream.write((char*)&size, sizeof(size));
          outStream.write((char*)&ipVec[0], ipVec.size() * sizeof(MinmerMapValueType::value_type));
        }
      }


      /**
       * @brief  Write posList for quick loading
       */
      void writeFreqKmersBinary(std::ofstream& outStream) 
      {
        typename MI_Map_t::size_type size = frequentSeeds.size();
        outStream.write((char*)&size, sizeof(size));

        for (hash_t kmerHash : frequentSeeds) 
        {
          MinmerMapKeyType key = kmerHash;
          outStream.write((char*)&key, sizeof(key));
        }
      }


      /**
       * @brief  Write all index data structures to disk
       */
      void writeIndex() 
      {
        fs::path freqListFilename = fs::path(param.indexFilename);
        std::ofstream outStream;
        outStream.open(freqListFilename, std::ios::binary);

        writeSketchBinary(outStream);
        writePosListBinary(outStream);
        writeFreqKmersBinary(outStream);
      }

      /**
       * @brief Read sketch from TSV file
       */
      void readSketchTSV() 
      {
        io::CSVReader<5, io::trim_chars<' '>, io::no_quote_escape<'\t'>> inReader(std::string(param.indexFilename) + ".tsv");
        inReader.read_header(io::ignore_missing_column, "seqId", "strand", "start", "end", "hash");
        hash_t hash;
        offset_t start, end;
        strand_t strand;
        seqno_t seqId;
        while (inReader.read_row(seqId, strand, start, end, hash))
        {
          this->minmerIndex.push_back(MinmerInfo {hash, start, end, seqId, strand});
        }
      }

      /**
       * @brief Read sketch from binary file
       */
      void readSketchBinary(std::ifstream& inStream) 
      {
        typename MI_Type::size_type size = 0;
        inStream.read((char*)&size, sizeof(size));
        minmerIndex.resize(size);
        inStream.read((char*)&minmerIndex[0], minmerIndex.size() * sizeof(MinmerInfo));
      }

      /**
       * @brief  Save posList for quick reading
       */
      void readPosListBinary(std::ifstream& inStream) 
      {
        typename MI_Map_t::size_type numKeys = 0;
        inStream.read((char*)&numKeys, sizeof(numKeys));
        minmerPosLookupIndex.reserve(numKeys);

        for (auto idx = 0; idx < numKeys; idx++) 
        {
          MinmerMapKeyType key = 0;
          inStream.read((char*)&key, sizeof(key));
          typename MinmerMapValueType::size_type size = 0;
          inStream.read((char*)&size, sizeof(size));

          minmerPosLookupIndex[key].resize(size);
          inStream.read((char*)&minmerPosLookupIndex[key][0], size * sizeof(MinmerMapValueType::value_type));
        }
      }


      /**
       * @brief  read frequent kmers from file
       */
      void readFreqKmersBinary(std::ifstream& inStream) 
      {
        typename MI_Map_t::size_type numKeys = 0;
        inStream.read((char*)&numKeys, sizeof(numKeys));
        frequentSeeds.reserve(numKeys);

        for (auto idx = 0; idx < numKeys; idx++) 
        {
          MinmerMapKeyType key = 0;
          inStream.read((char*)&key, sizeof(key));
          frequentSeeds.insert(key);
        }
      }


      /**
       * @brief  Read all index data structures from file
       */
      void readIndex() 
      { 
        fs::path indexFilename = fs::path(param.indexFilename);
        std::ifstream inStream;
        inStream.open(indexFilename, std::ios::binary);
        readSketchBinary(inStream);
        readPosListBinary(inStream);
        readFreqKmersBinary(inStream);
      }


      /**
       * @brief   report the frequency histogram of minmers using position lookup index
       *          and compute which high frequency minmers to ignore
       */
      void computeFreqHist()
      {
          if (!this->minmerPosLookupIndex.empty()) {
              //1. Compute histogram

              for (auto& [kmerHash, freq] : this->hashFreq)
                  this->minmerFreqHistogram[freq]++;

              std::cerr << "[mashmap::skch::Sketch::computeFreqHist] Frequency histogram of minmer intervals = "
                        << *this->minmerFreqHistogram.begin() << " ... " << *this->minmerFreqHistogram.rbegin()
                        << std::endl;

              //2. Compute frequency threshold to ignore most frequent minmers

              int64_t totalUniqueMinmers = this->hashFreq.size();
              int64_t minmerToIgnore = totalUniqueMinmers * param.kmer_pct_threshold / 100;

              int64_t sum = 0;

              //Iterate from highest frequent minmers
              for (auto it = this->minmerFreqHistogram.rbegin(); it != this->minmerFreqHistogram.rend(); it++) {
                  sum += it->second; //add frequency
                  if (sum < minmerToIgnore) {
                      this->freqThreshold = it->first;
                      //continue
                  } else if (sum == minmerToIgnore) {
                      this->freqThreshold = it->first;
                      break;
                  } else {
                      break;
                  }
              }

              if (this->freqThreshold != std::numeric_limits<uint64_t>::max())
                  std::cerr << "[mashmap::skch::Sketch::computeFreqHist] With threshold " << this->param.kmer_pct_threshold
                            << "\%, ignore minmers with more than >= " << this->freqThreshold << " intervals during mapping."
                            << std::endl;
              else
                  std::cerr << "[mashmap::skch::Sketch::computeFreqHist] With threshold " << this->param.kmer_pct_threshold
                            << "\%, consider all minmers during mapping." << std::endl;
          } else {
              std::cerr << "[mashmap::skch::Sketch::computeFreqHist] No minmers." << std::endl;
          }
      }

      public:

      /**
       * @brief               search hash associated with given position inside the index
       * @details             if MIIter_t iter is returned, than *iter's wpos >= winpos
       * @param[in]   seqId
       * @param[in]   winpos
       * @return              iterator to the minmer in the index
       */

      /**
       * @brief                 check if iterator points to index end
       * @param[in]   iterator
       * @return                boolean value
       */
      bool isMinmerIndexEnd(const MIIter_t &it) const
      {
        return it == this->minmerIndex.end();
      }

      /**
       * @brief     Return end iterator on minmerIndex
       */
      MIIter_t getMinmerIndexEnd() const
      {
        return this->minmerIndex.end();
      }

      int getFreqThreshold() const
      {
        return this->freqThreshold;
      }

      void computeFreqSeedSet()
      {
        for(auto& [kmerHash, frequency] : this->hashFreq) {
          if (frequency >= this->freqThreshold) {
            this->frequentSeeds.insert(kmerHash);
          }
        }
      }

      void dropFreqSeedSet()
      {
        this->minmerIndex.erase(
          std::remove_if(minmerIndex.begin(), minmerIndex.end(), [&] 
            (auto& mi) {return isFreqSeed(mi.hash);}
          ), minmerIndex.end()
        );
        for (hash_t kmer : this->frequentSeeds)
        {
          this->minmerPosLookupIndex.erase(kmer);
        }
      }

      bool isFreqSeed(hash_t h) const
      {
        return frequentSeeds.find(h) != frequentSeeds.end();
      }

    }; //End of class Sketch
} //End of namespace skch

#endif
