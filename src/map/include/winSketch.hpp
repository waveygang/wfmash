/**
 * @file    winSketch.hpp
 * @brief   routines to index the reference 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WIN_SKETCH_HPP 
#define WIN_SKETCH_HPP

#include <algorithm>
#include <cassert>
#include <cstring>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
static inline bool hasSuffix(const std::string &s, const std::string &suffix) {
  return s.size() >= suffix.size() &&
         s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

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
      int freqThreshold = std::numeric_limits<int>::max();

      //Set of frequent seeds to be ignored
      ankerl::unordered_dense::set<hash_t> frequentSeeds;

      //Make the default constructor private, non-accessible
      Sketch();

      public:

      typedef std::vector< MinmerInfo > MI_Type;
      using MIIter_t = MI_Type::const_iterator;

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

      // CSR-flattened position lookup, built by flattenPosLookup() at the end of
      // construction when every position fits the packed uint64 key (packed_ok):
      // per unique minmer an (offset,count) run into ipArena of encodePackedIP
      // keys, per-key point order preserved verbatim. Frequent seeds' runs are
      // dropped (their lookups can never happen: getSeedHits removes them from
      // the query table first). minmerPosLookupIndex is destroyed after
      // flattening; when !packed_ok it stays and consumers use it as before.
      bool packed_ok = false;
      std::vector<uint64_t> ipArena;
      struct PosLookupRun {
        uint64_t off;
        uint32_t cnt;
        uint32_t written;   // pass-2 fill cursor; equals cnt afterwards
      };
      ankerl::unordered_dense::map<MinmerMapKeyType, PosLookupRun> posLookupCSR;
      size_t nUniqueMinmers = 0;   // minmerPosLookupIndex.size() before flattening

      private:

      /**
       * Keep list of minmers, sequence# , their position within seq , here while parsing sequence 
       * Note : position is local within each contig
       * Hashes saved here are non-unique, ordered as they appear in the reference
       */

      //Frequency histogram of minmers
      //[... ,x -> y, ...] implies y number of minmers occur x times
      std::map<int, int> minmerFreqHistogram;

      public:

      /**
       * @brief   constructor
       *          also builds, indexes the minmer table
       */
      Sketch(const skch::Parameters &p) 
        :
          param(p) {
            this->build();
            this->index();
            if (!param.saveIndexFilename.empty()) {
              if (hasSuffix(param.saveIndexFilename, ".tsv")) {
                this->saveIndexTSV();
              } else {
                this->saveIndexBinary();
              }
              this->savePosListBinary();
            }
            this->computeFreqHist();
            this->computeFreqSeedSet();
            this->dropFreqSeedSet();
          }

      private:

      /**
       * @brief     build the sketch table
       * @details   compute and save minmers from the reference sequence(s)
       *            assuming a fixed window size
       */
      void build()
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
        if (!param.loadIndexFilename.empty()) {
          if (hasSuffix(param.loadIndexFilename, ".tsv")) {
            this->loadIndexTSV();
          } else {
            this->loadIndexBinary();
          }
        }

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
                  if (param.loadIndexFilename.empty()) {
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

        if (param.loadIndexFilename.empty()) {
          //Collect remaining output objects
          while ( threadPool.running() )
            this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
        }

        std::cerr << "[mashmap::skch::Sketch::build] minmer windows picked from reference = " << minmerIndex.size() << std::endl;

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
      void buildHandleThreadOutput(MI_Type* output)
      {
        this->minmerIndex.insert(this->minmerIndex.end(), output->begin(), output->end());
        delete output;
      }


      /**
       * @brief  Save index. TSV indexing is slower but can be debugged easier
       */
      void saveIndexTSV() 
      {
        std::ofstream outStream;
        outStream.open(param.saveIndexFilename);
        outStream << "seqId" << "\t" << "strand" << "\t" << "start" << "\t" << "end" << "\t" << "hash\n";
        for (auto& mi : this->minmerIndex) {
          outStream << mi.seqId << "\t" << std::to_string(mi.strand) << "\t" << mi.wpos << "\t" << mi.wpos_end << "\t" << mi.hash << "\n";
        }
        outStream.close(); 
      }

      /**
       * @brief  Save index for quick loading
       */
      void saveIndexBinary() 
      {
        std::string indexFilename = param.saveIndexFilename;
        indexFilename += ".index";
        std::ofstream outStream;
        outStream.open(indexFilename, std::ios::binary);
        typename MI_Type::size_type size = minmerIndex.size();
        outStream.write((char*)&size, sizeof(size));
        outStream.write((char*)&minmerIndex[0], minmerIndex.size() * sizeof(MinmerInfo));
      }

      /**
       * @brief  Save posList for quick loading
       */
      void savePosListBinary() 
      {
        std::string posListFilename = param.saveIndexFilename;
        posListFilename += ".map";
        std::ofstream outStream;
        outStream.open(posListFilename, std::ios::binary);
        typename MI_Map_t::size_type size = packed_ok ? posLookupCSR.size() : minmerPosLookupIndex.size();
        outStream.write((char*)&size, sizeof(size));

        if (packed_ok)
        {
          // Runs saved pre-pruning (this runs before computeFreqHist), decoded
          // back to the legacy on-disk record layout.
          std::vector<IntervalPoint> ipVec;
          for (auto& [hash, run] : posLookupCSR)
          {
            MinmerMapKeyType key = hash;
            outStream.write((char*)&key, sizeof(key));
            typename MI_Type::size_type sz = run.cnt;
            outStream.write((char*)&sz, sizeof(sz));
            ipVec.clear();
            for (uint32_t i = 0; i < run.cnt; ++i) {
              IntervalPoint ip = decodePackedIP(ipArena[run.off + i]);
              ip.hash = hash;
              ipVec.push_back(ip);
            }
            outStream.write((char*)&ipVec[0], ipVec.size() * sizeof(MinmerMapValueType::value_type));
          }
          return;
        }

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
       * @brief Load index from TSV file
       */
      void loadIndexTSV() 
      {
        io::CSVReader<5, io::trim_chars<' '>, io::no_quote_escape<'\t'>> inReader(param.loadIndexFilename);
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
       * @brief Load index from binary file
       */
      void loadIndexBinary() 
      {
        std::string indexFilename = param.loadIndexFilename;
        indexFilename += ".index";
        std::ifstream inStream;
        inStream.open(indexFilename, std::ios::binary);
        typename MI_Type::size_type size = 0;
        inStream.read((char*)&size, sizeof(size));
        minmerIndex.resize(size);
        inStream.read((char*)&minmerIndex[0], minmerIndex.size() * sizeof(MinmerInfo));
      }

      /**
       * @brief  Save posList for quick loading
       */
      void loadPosListBinary() 
      {
        std::string posListFilename = param.loadIndexFilename;
        posListFilename += ".map";
        std::ifstream inStream;
        inStream.open(posListFilename, std::ios::binary);
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
       * @brief   build the index for fast lookups using minmer table
       */
      void index()
      {
        // The packed CSR representation needs every position to fit 32 bits and
        // seqIds 31 bits; this depends only on metadata, so decide here.
        packed_ok = metadata.size() <= (size_t)0x7FFFFFFF;
        for (const auto& m : metadata) {
          if (m.len < 0 || (uint64_t)m.len >= (UINT64_C(1) << 32)) { packed_ok = false; break; }
        }
        //Parse all the minmers and push into the map
        //minmerPosLookupIndex.set_empty_key(0);
        if (param.loadIndexFilename.empty())
        {
          if (packed_ok)
          {
            this->buildPosLookupCSR();
          }
          else
          {
            for(auto &mi : minmerIndex)
            {
              // [hash value -> info about minmer]
              auto& ipVec = minmerPosLookupIndex[mi.hash];
              if (ipVec.size() == 0
                  || ipVec.back().hash != mi.hash
                  || ipVec.back().pos != mi.wpos)
              {
                ipVec.push_back(IntervalPoint {mi.wpos, mi.hash, mi.seqId, side::OPEN});
                ipVec.push_back(IntervalPoint {mi.wpos_end, mi.hash, mi.seqId, side::CLOSE});
              } else {
                ipVec.back().pos = mi.wpos_end;
              }
            }
          }
        }
        else
        {
          // Loaded legacy-format posList: keep the legacy map path.
          packed_ok = false;
          this->loadPosListBinary();
        }
        nUniqueMinmers = packed_ok ? posLookupCSR.size() : minmerPosLookupIndex.size();
        std::cerr << "[mashmap::skch::Sketch::index] unique minmers = " << nUniqueMinmers << std::endl;
      }

      // Build the position lookup directly as a CSR arena of packed uint64 keys
      // (no per-key vectors): pass 1 counts each key's points with index()'s
      // run-merging rule, pass 2 fills the runs in the identical per-key order
      // the legacy vectors would have had.
      void buildPosLookupCSR()
      {
        for (const auto& mi : minmerIndex) {
          auto [it, fresh] = posLookupCSR.try_emplace(mi.hash, PosLookupRun{0, 0, 0});
          auto& run = it->second;
          // run.off doubles as the last CLOSE position during this pass
          if (run.cnt == 0 || (offset_t)run.off != mi.wpos) {
            run.cnt += 2;
          }
          run.off = (uint64_t)mi.wpos_end;
        }
        uint64_t total = 0;
        for (auto& e : posLookupCSR) {
          const uint32_t c = e.second.cnt;
          e.second.off = total;
          e.second.written = 0;
          total += c;
        }
        ipArena.assign(total, 0);
        for (const auto& mi : minmerIndex) {
          auto& run = posLookupCSR[mi.hash];
          if (run.written == 0
              || decodePackedIP(ipArena[run.off + run.written - 1]).pos != mi.wpos) {
            ipArena[run.off + run.written++] = encodePackedIP(IntervalPoint{mi.wpos, mi.hash, mi.seqId, side::OPEN});
            ipArena[run.off + run.written++] = encodePackedIP(IntervalPoint{mi.wpos_end, mi.hash, mi.seqId, side::CLOSE});
          } else {
            // Merge extends the existing CLOSE point's position but keeps its
            // seqId (index() mutated only .pos; merges can chain across seqIds).
            uint64_t& last = ipArena[run.off + run.written - 1];
            last = (last & ~((uint64_t)0xFFFFFFFFULL << 1)) | ((uint64_t)mi.wpos_end << 1);
          }
        }
      }

      /**
       * @brief   report the frequency histogram of minmers using position lookup index
       *          and compute which high frequency minmers to ignore
       */
      void computeFreqHist()
      {
          if (packed_ok ? !posLookupCSR.empty() : !minmerPosLookupIndex.empty()) {
              //1. Compute histogram

              if (packed_ok) {
                for (auto &e : this->posLookupCSR)
                    this->minmerFreqHistogram[e.second.cnt] += 1;
              } else {
                for (auto &e : this->minmerPosLookupIndex)
                    this->minmerFreqHistogram[e.second.size()] += 1;
              }

              std::cerr << "[mashmap::skch::Sketch::computeFreqHist] Frequency histogram of minmer interval points = "
                        << *this->minmerFreqHistogram.begin() << " ... " << *this->minmerFreqHistogram.rbegin()
                        << std::endl;

              //2. Compute frequency threshold to ignore most frequent minmers

              int64_t totalUniqueMinmers = this->nUniqueMinmers;
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

              if (this->freqThreshold != std::numeric_limits<int>::max())
                  std::cerr << "[mashmap::skch::Sketch::computeFreqHist] With threshold " << this->param.kmer_pct_threshold
                            << "\%, ignore minmers occurring >= " << this->freqThreshold << " times during lookup."
                            << std::endl;
              else
                  std::cerr << "[mashmap::skch::Sketch::computeFreqHist] With threshold " << this->param.kmer_pct_threshold
                            << "\%, consider all minmers during lookup." << std::endl;
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
        if (packed_ok) {
          for(auto &e : this->posLookupCSR) {
            if (e.second.cnt >= (uint32_t)this->freqThreshold) {
              this->frequentSeeds.insert(e.first);
            }
          }
        } else {
          for(auto &e : this->minmerPosLookupIndex) {
            if (e.second.size() >= this->freqThreshold) {
              this->frequentSeeds.insert(e.first);
            }
          }
        }
      }

      void dropFreqSeedSet()
      {
        this->minmerIndex.erase(
          std::remove_if(minmerIndex.begin(), minmerIndex.end(), [&]
            (auto& mi) {return this->frequentSeeds.find(mi.hash) != this->frequentSeeds.end();}
          ), minmerIndex.end()
        );
        if (packed_ok && !frequentSeeds.empty()) {
          // Compact the arena over the surviving keys (iteration order equals
          // ascending offset order: the map is untouched since buildPosLookupCSR),
          // then drop the frequent keys. Their lookups can never happen anyway:
          // getSeedHits removes frequent hashes from the query table first.
          uint64_t w = 0;
          for (auto& e : posLookupCSR) {
            if (frequentSeeds.find(e.first) != frequentSeeds.end()) continue;
            if (e.second.off != w) {
              std::memmove(ipArena.data() + w, ipArena.data() + e.second.off,
                           (size_t)e.second.cnt * sizeof(uint64_t));
            }
            e.second.off = w;
            w += e.second.cnt;
          }
          ipArena.resize(w);
          for (const auto& h : frequentSeeds) posLookupCSR.erase(h);
        }
      }


      bool isFreqSeed(hash_t h) const
      {
        return frequentSeeds.find(h) != frequentSeeds.end();
      }

    }; //End of class Sketch
} //End of namespace skch

#endif
