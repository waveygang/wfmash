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
#include "common/atomic_queue/atomic_queue.h"
#include "sequenceIds.hpp"
#include "common/atomic_queue/atomic_queue.h"
#include "common/progress.hpp"
#include <thread>
#include <atomic>

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
      skch::Parameters param;

      //Minmers that occur this or more times will be ignored (computed based on percentageThreshold)
      uint64_t freqThreshold = std::numeric_limits<uint64_t>::max();

      //Set of frequent seeds to be ignored
      ankerl::unordered_dense::set<hash_t> frequentSeeds;

      //Make the default constructor protected, non-accessible
      protected:
      Sketch(SequenceIdManager& idMgr) : idManager(idMgr) {}

      public:

      //Flag to indicate if the Sketch is fully initialized
      bool isInitialized = false;

      using MI_Type = std::vector< MinmerInfo >;
      using MIIter_t = MI_Type::const_iterator;
      using HF_Map_t = ankerl::unordered_dense::map<hash_t, uint64_t>;

      // Frequency of each hash
      HF_Map_t hashFreq;

      public:
        uint64_t total_seq_length = 0;

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

      // Atomic queues for input and output
      using input_queue_t = atomic_queue::AtomicQueue<InputSeqContainer*, 1024>;
      using output_queue_t = atomic_queue::AtomicQueue<std::pair<uint64_t, MI_Type*>*, 1024>;
      input_queue_t input_queue;
      output_queue_t output_queue;

      private:

      /**
       * Keep list of minmers, sequence# , their position within seq , here while parsing sequence 
       * Note : position is local within each contig
       * Hashes saved here are non-unique, ordered as they appear in the reference
       */

      //Frequency histogram of minmers
      //[... ,x -> y, ...] implies y number of minmers occur x times
      std::map<uint64_t, uint64_t> minmerFreqHistogram;

      // Instance of the SequenceIdManager
      SequenceIdManager& idManager;

      public:

      /**
       * @brief   constructor
       *          also builds, indexes the minmer table
       */
      Sketch(skch::Parameters p,
             SequenceIdManager& idMgr,
             const std::vector<std::string>& targets = {},
             const std::string& indexFilename = "")
        : param(std::move(p)),
          idManager(idMgr)
      {
        if (!indexFilename.empty()) {
          loadIndex(indexFilename);
        } else {
          initialize(targets);
        }
      }

      void loadIndex(const std::string& indexFilename) {
        std::ifstream inStream(indexFilename, std::ios::binary);
        if (!inStream) {
          std::cerr << "Error: Unable to open index file: " << indexFilename << std::endl;
          exit(1);
        }
        readParameters(inStream);
        readSketchBinary(inStream);
        readPosListBinary(inStream);
        readFreqKmersBinary(inStream);
        inStream.close();
        isInitialized = true;
        std::cerr << "[mashmap::skch::Sketch] Sketch loaded from index file: " << indexFilename << std::endl;
      }

    public:
      void initialize(const std::vector<std::string>& targets = {}) {
        std::cerr << "[mashmap::skch::Sketch] Initializing Sketch..." << std::endl;
        
        // Calculate total sequence length
        /*
        for (const auto& fileName : param.refSequences) {
            std::cerr << "targets are " << targets.size() << " ";
            for (const auto& target : targets) {
                std::cerr << target << " ";
            }
            std::cerr << std::endl;
            seqiter::for_each_seq_in_file(
                fileName,
                targets,
                [&](const std::string& seq_name, const std::string& seq) {
                    total_seq_length += seq.length();
                });
        }
        */
        
        this->build(true, targets);
        this->computeFreqHist();
        this->computeFreqSeedSet();
        this->dropFreqSeedSet();
        this->hashFreq.clear();
        if (!param.indexFilename.empty())
        {
            this->writeIndex();
        }

        std::cerr << "[mashmap::skch::Sketch] Unique minmer hashes after pruning = " << (minmerPosLookupIndex.size() - this->frequentSeeds.size()) << std::endl;
        std::cerr << "[mashmap::skch::Sketch] Total minmer windows after pruning = " << minmerIndex.size() << std::endl;
        std::cerr << "[mashmap::skch::Sketch] Number of sequences = " << idManager.size() << std::endl;
        isInitialized = true;
        std::cerr << "[mashmap::skch::Sketch] Sketch initialization complete." << std::endl;
      }

      private:
      void reader_thread(const std::vector<std::string>& targets, std::atomic<bool>& reader_done) {
          for (const auto& fileName : param.refSequences) {
              seqiter::for_each_seq_in_file(
                  fileName,
                  targets,
                  [&](const std::string& seq_name, const std::string& seq) {
                      if (seq.length() >= param.segLength) {
                          seqno_t seqId = idManager.getSequenceId(seq_name);
                          auto record = new InputSeqContainer(seq, seq_name, seqId);
                          input_queue.push(record);
                      }
                      // We don't update progress here anymore
                  });
          }
          reader_done.store(true);
      }

      void worker_thread(std::atomic<bool>& reader_done, progress_meter::ProgressMeter& progress) {
          while (true) {
              InputSeqContainer* record = nullptr;
              if (input_queue.try_pop(record)) {
                  auto minmers = new MI_Type();
                  CommonFunc::addMinmers(*minmers, &(record->seq[0]), record->len, 
                                         param.kmerSize, param.segLength, param.alphabetSize, 
                                         param.sketchSize, record->seqId);
                  auto output_pair = new std::pair<uint64_t,MI_Type*>(record->len, minmers);
                  output_queue.push(output_pair);
                  delete record;
              } else if (reader_done.load() && input_queue.was_empty()) {
                  break;
              } else {
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }
      }

      void writer_thread(std::atomic<bool>& workers_done, progress_meter::ProgressMeter& progress) {
          while (true) {
              std::pair<uint64_t, MI_Type*>* output = nullptr;
              if (output_queue.try_pop(output)) {
                  uint64_t seq_length = output->first;
                  MI_Type* minmers = output->second;
                  for (const auto& mi : *minmers) {
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
                  //this->minmerIndex.insert(this->minmerIndex.end(), minmers->begin(), minmers->end());
                  this->minmerIndex.insert(
                      this->minmerIndex.end(), 
                      std::make_move_iterator(minmers->begin()), 
                      std::make_move_iterator(minmers->end()));

                  // Update progress meter
                  progress.increment(seq_length);
                  
                  delete output->second;
                  delete output;
              } else if (workers_done.load() && output_queue.was_empty()) {
                  break;
              } else {
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
              }
          }

          // Finalize progress meter
          progress.finish();
      }

      /**
       * @brief     Get sequence metadata and optionally build the sketch table
       *
       * @details   Iterate through ref sequences to get metadata and
       *            optionally compute and save minmers from the reference sequence(s)
       *            assuming a fixed window size
       * @param     compute_seeds   Whether to compute seeds or just collect metadata
       * @param     target_ids      Set of target sequence IDs to sketch over
       */
      void build(bool compute_seeds, const std::vector<std::string>& target_names = {})
      {
        std::chrono::time_point<std::chrono::system_clock> t0 = skch::Time::now();

        if (compute_seeds) {
            std::cerr << "creating seeds" << std::endl;

          //Create the thread pool 
          ThreadPool<InputSeqContainer, MI_Type> threadPool([this](InputSeqContainer* e) { return buildHelper(e); }, param.threads);

          size_t totalSeqProcessed = 0;
          size_t totalSeqSkipped = 0;
          size_t shortestSeqLength = std::numeric_limits<size_t>::max();
          for (const auto& fileName : param.refSequences) {
            std::cerr << "[mashmap::skch::Sketch::build] Processing file: " << fileName << std::endl;

            seqiter::for_each_seq_in_file(
              fileName,
              target_names,
              [&](const std::string& seq_name, const std::string& seq) {
                  std::cerr << "on sequence " << seq_name << std::endl;
                  if (seq.length() >= param.segLength) {
                      seqno_t seqId = idManager.getSequenceId(seq_name);
                      threadPool.runWhenThreadAvailable(new InputSeqContainer(seq, seq_name, seqId));
                      totalSeqProcessed++;
                      shortestSeqLength = std::min(shortestSeqLength, seq.length());
                      std::cerr << "DEBUG: Processing sequence: " << seq_name << " (length: " << seq.length() << ")" << std::endl;

                      //Collect output if available
                      while (threadPool.outputAvailable()) {
                          this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
                      }
                    
                      // Update metadata
                      // Metadata is now handled by idManager, no need to push_back here
                  } else {
                      totalSeqSkipped++;
                      std::cerr << "WARNING, skch::Sketch::build, skipping short sequence: " << seq_name << " (length: " << seq.length() << ")" << std::endl;
                  }
              });
          }
          
          // Update sequencesByFileInfo
          // Removed as sequencesByFileInfo is no longer used
          std::cerr << "[mashmap::skch::Sketch::build] Shortest sequence length: " << shortestSeqLength << std::endl;

          //Collect remaining output objects
          while (threadPool.running())
            this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());

          std::cerr << "[mashmap::skch::Sketch::build] Total sequences processed: " << totalSeqProcessed << std::endl;
          std::cerr << "[mashmap::skch::Sketch::build] Total sequences skipped: " << totalSeqSkipped << std::endl;
          std::cerr << "[mashmap::skch::Sketch::build] Unique minmer hashes before pruning = " << minmerPosLookupIndex.size() << std::endl;
          std::cerr << "[mashmap::skch::Sketch::build] Total minmer windows before pruning = " << minmerIndex.size() << std::endl;
        }

        std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
        std::cerr << "[mashmap::skch::Sketch::build] time spent computing the reference index: " << timeRefSketch.count() << " sec" << std::endl;

        if (this->minmerIndex.size() == 0)
        {
          std::cerr << "[mashmap::skch::Sketch::build] ERROR, reference sketch is empty. Reference sequences shorter than the kmer size are not indexed" << std::endl;
          exit(1);
        }
      }

      public:

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
                input->seqId);

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
       * @brief Write parameters 
       */
      void writeParameters(std::ofstream& outStream)
      {
        // Write segment length, sketch size, and kmer size
        outStream.write((char*) &param.segLength, sizeof(param.segLength));
        outStream.write((char*) &param.sketchSize, sizeof(param.sketchSize));
        outStream.write((char*) &param.kmerSize, sizeof(param.kmerSize));
      }


      /**
       * @brief  Write all index data structures to disk
       */
      void writeIndex(const std::string& filename = "") 
      {
        fs::path freqListFilename = filename.empty() ? fs::path(param.indexFilename) : fs::path(filename);
        std::ofstream outStream;
        outStream.open(freqListFilename, std::ios::binary);

        writeParameters(outStream);
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
       * @brief  Read parameters and compare to CLI params
       */
      void readParameters(std::ifstream& inStream)
      {
        // Read segment length, sketch size, and kmer size
        decltype(param.segLength) index_segLength;
        decltype(param.sketchSize) index_sketchSize;
        decltype(param.kmerSize) index_kmerSize;

        inStream.read((char*) &index_segLength, sizeof(index_segLength));
        inStream.read((char*) &index_sketchSize, sizeof(index_sketchSize));
        inStream.read((char*) &index_kmerSize, sizeof(index_kmerSize));

        if (param.segLength != index_segLength 
            || param.sketchSize != index_sketchSize
            || param.kmerSize != index_kmerSize)
        {
          std::cerr << "[mashmap::skch::Sketch::build] ERROR: Parameters of indexed sketch differ from CLI parameters" << std::endl;
          std::cerr << "[mashmap::skch::Sketch::build] ERROR: Index --> segLength=" << index_segLength
            << " sketchSize=" << index_sketchSize << " kmerSize=" << index_kmerSize << std::endl;
          std::cerr << "[mashmap::skch::Sketch::build] ERROR: CLI   --> segLength=" << param.segLength
            << " sketchSize=" << param.sketchSize << " kmerSize=" << param.kmerSize << std::endl;
          exit(1);
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
        readParameters(inStream);
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

              for (auto& e : this->minmerPosLookupIndex)
                  this->minmerFreqHistogram[e.second.size()]++;

              std::cerr << "[mashmap::skch::Sketch::computeFreqHist] Frequency histogram of minmer interval points = "
                        << *this->minmerFreqHistogram.begin() << " ... " << *this->minmerFreqHistogram.rbegin()
                        << std::endl;

              //2. Compute frequency threshold to ignore most frequent minmers

              int64_t totalUniqueMinmers = this->minmerPosLookupIndex.size();
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
                            << "\%, ignore minmers with more than >= " << this->freqThreshold << " interval points during mapping."
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
        for(auto &e : this->minmerPosLookupIndex) {
          if (e.second.size() >= this->freqThreshold) {
            this->frequentSeeds.insert(e.first);
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
      }

      bool isFreqSeed(hash_t h) const
      {
        return frequentSeeds.find(h) != frequentSeeds.end();
      }

      void clear()
      {
        hashFreq.clear();
        minmerPosLookupIndex.clear();
        minmerIndex.clear();
        minmerFreqHistogram.clear();
        frequentSeeds.clear();
        freqThreshold = std::numeric_limits<uint64_t>::max();
      }

    }; //End of class Sketch
} //End of namespace skch

#endif
