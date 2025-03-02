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

      //Make the default constructor protected, non-accessible
      protected:
      Sketch(SequenceIdManager& idMgr) : idManager(idMgr) {}

      public:

      //Flag to indicate if the Sketch is fully initialized
      bool isInitialized = false;

      using MI_Type = std::vector< MinmerInfo >;
      using MIIter_t = MI_Type::const_iterator;
      using HF_Map_t = ankerl::unordered_dense::map<hash_t, uint64_t>;

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

      double hgNumerator;

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
             std::ifstream* indexStream = nullptr)
        : param(std::move(p)),
          idManager(idMgr)
      {
        if (indexStream) {
          readIndex(*indexStream, targets);
        } else {
          initialize(targets);
        }
      }

    public:
      void initialize(const std::vector<std::string>& targets = {}) {
        this->build(true, targets);
        this->hgNumerator = param.hgNumerator;
        isInitialized = true;
      }

    private:

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
          // Calculate total sequence length from id manager
          uint64_t total_seq_length = 0;
          for (const auto& seqName : target_names) {
              seqno_t seqId = idManager.getSequenceId(seqName);
              total_seq_length += idManager.getSequenceLength(seqId);
          }

          // First progress meter for sketch computation
          progress_meter::ProgressMeter sketch_progress(
              total_seq_length,
              "[wfmash::mashmap] computing sketch");

          // Create the thread pool 
          ThreadPool<InputSeqContainer, MI_Type> threadPool(
              [this, &sketch_progress](InputSeqContainer* e) { 
                  return buildHelper(e, &sketch_progress); 
              }, 
              param.threads);

          size_t totalSeqProcessed = 0;
          size_t totalSeqSkipped = 0;
          size_t shortestSeqLength = std::numeric_limits<size_t>::max();
          
          // Vector to store all thread outputs
          std::vector<MI_Type*> threadOutputs;

          for (const auto& fileName : param.refSequences) {
              seqiter::for_each_seq_in_file(
                  fileName,
                  target_names,
                  [&](const std::string& seq_name, const std::string& seq) {
                      if (seq.length() >= param.segLength) {
                          seqno_t seqId = idManager.getSequenceId(seq_name);
                          threadPool.runWhenThreadAvailable(new InputSeqContainer(seq, seq_name, seqId));
                          totalSeqProcessed++;
                          shortestSeqLength = std::min(shortestSeqLength, seq.length());

                          while (threadPool.outputAvailable()) {
                              auto output = threadPool.popOutputWhenAvailable();
                              threadOutputs.push_back(output);
                          }
                      } else {
                          totalSeqSkipped++;
                          std::cerr << "WARNING, skch::Sketch::build, skipping short sequence: " << seq_name 
                                   << " (length: " << seq.length() << ")" << std::endl;
                      }
                  });
          }

          while (threadPool.running()) {
              auto output = threadPool.popOutputWhenAvailable();
              threadOutputs.push_back(output);
          }

          // Make sure to finish first progress meter before starting the next
          sketch_progress.finish();

          // Calculate total windows for index building progress
          uint64_t total_windows = 0;
          for (const auto& output : threadOutputs) {
              total_windows += output->size();
          }

          // Second progress meter for index building
          progress_meter::ProgressMeter index_progress(
              total_windows,
              "[wfmash::mashmap] building index");

          // Parallel k-mer frequency counting
          std::vector<HF_Map_t> thread_kmer_freqs(param.threads);
          std::vector<std::thread> freq_threads;
          
          // Split outputs into chunks for parallel processing
          size_t chunk_size = (threadOutputs.size() + param.threads - 1) / param.threads;
          
          for (size_t t = 0; t < param.threads; ++t) {
              freq_threads.emplace_back([&, t]() {
                  size_t start = t * chunk_size;
                  size_t end = std::min(start + chunk_size, threadOutputs.size());
                  
                  for (size_t i = start; i < end; ++i) {
                      for (const MinmerInfo& mi : *threadOutputs[i]) {
                          thread_kmer_freqs[t][mi.hash]++;
                      }
                  }
              });
          }
          
          for (auto& thread : freq_threads) {
              thread.join();
          }

          // Merge frequency maps
          HF_Map_t kmer_freqs;
          for (const auto& thread_freq : thread_kmer_freqs) {
              for (const auto& [hash, freq] : thread_freq) {
                  kmer_freqs[hash] += freq;
              }
          }

          // Parallel index building
          std::vector<MI_Map_t> thread_pos_indexes(param.threads);
          std::vector<MI_Type> thread_minmer_indexes(param.threads);
          std::vector<uint64_t> thread_total_kmers(param.threads, 0);
          std::vector<uint64_t> thread_filtered_kmers(param.threads, 0);
          std::vector<std::thread> index_threads;

          for (size_t t = 0; t < param.threads; ++t) {
              index_threads.emplace_back([&, t]() {
                  size_t start = t * chunk_size;
                  size_t end = std::min(start + chunk_size, threadOutputs.size());

                  for (size_t i = start; i < end; ++i) {
                      for (const MinmerInfo& mi : *threadOutputs[i]) {
                          thread_total_kmers[t]++;
                          
                          auto freq_it = kmer_freqs.find(mi.hash);
                          if (freq_it == kmer_freqs.end()) {
                              continue;  // Should never happen
                          }

                          uint64_t freq = freq_it->second;
                          uint64_t min_occ = 10;
                          uint64_t max_occ = std::numeric_limits<uint64_t>::max();
                          uint64_t count_threshold;
                          
                          if (param.max_kmer_freq <= 1.0) {
                              count_threshold = std::min(max_occ, 
                                                       std::max(min_occ, 
                                                              (uint64_t)(total_windows * param.max_kmer_freq)));
                          } else {
                              count_threshold = std::min(max_occ,
                                                       std::max(min_occ,
                                                              (uint64_t)param.max_kmer_freq));
                          }

                          if (freq > count_threshold && freq > min_occ) {
                              thread_filtered_kmers[t]++;
                              continue;
                          }

                          auto& pos_list = thread_pos_indexes[t][mi.hash];
                          if (pos_list.size() == 0 
                                  || pos_list.back().hash != mi.hash 
                                  || pos_list.back().pos != mi.wpos) {
                              pos_list.push_back(IntervalPoint {mi.wpos, mi.hash, mi.seqId, side::OPEN});
                              pos_list.push_back(IntervalPoint {mi.wpos_end, mi.hash, mi.seqId, side::CLOSE});
                          } else {
                              pos_list.back().pos = mi.wpos_end;
                          }

                          thread_minmer_indexes[t].push_back(mi);
                          index_progress.increment(1);
                      }
                      delete threadOutputs[i];
                  }
              });
          }

          for (auto& thread : index_threads) {
              thread.join();
          }

          // Merge results
          uint64_t total_kmers = std::accumulate(thread_total_kmers.begin(), thread_total_kmers.end(), 0ULL);
          uint64_t filtered_kmers = std::accumulate(thread_filtered_kmers.begin(), thread_filtered_kmers.end(), 0ULL);

          // Clear and resize main indexes
          minmerPosLookupIndex.clear();
          minmerIndex.clear();
          
          // Reserve approximate space
          size_t total_minmers = 0;
          for (const auto& thread_index : thread_minmer_indexes) {
              total_minmers += thread_index.size();
          }
          minmerIndex.reserve(total_minmers);

          // Merge position lookup indexes
          for (auto& thread_pos_index : thread_pos_indexes) {
              for (auto& [hash, pos_list] : thread_pos_index) {
                  auto& main_pos_list = minmerPosLookupIndex[hash];
                  main_pos_list.insert(main_pos_list.end(), pos_list.begin(), pos_list.end());
              }
          }

          // Merge minmer indexes
          for (auto& thread_index : thread_minmer_indexes) {
              minmerIndex.insert(minmerIndex.end(), 
                               std::make_move_iterator(thread_index.begin()),
                               std::make_move_iterator(thread_index.end()));
          }
          
          // Finish second progress meter
          index_progress.finish();

          uint64_t freq_cutoff;
          if (param.max_kmer_freq <= 1.0) {
              freq_cutoff = std::max(1UL, (uint64_t)(total_windows * param.max_kmer_freq));
          } else {
              freq_cutoff = (uint64_t)param.max_kmer_freq;
          }
          std::cerr << "[wfmash::mashmap] Processed " << totalSeqProcessed << " sequences (" << totalSeqSkipped << " skipped, " << total_seq_length << " total bp), " 
                    << minmerPosLookupIndex.size() << " unique hashes, " << minmerIndex.size() << " windows" << std::endl
                    << "[wfmash::mashmap] Filtered " << filtered_kmers << "/" << total_kmers 
                    << " k-mers occurring > " << freq_cutoff << " times"
                    << " (target: " << (param.max_kmer_freq <= 1.0 ? 
                                      ([&]() { 
                                          std::stringstream ss;
                                          ss << std::fixed << std::setprecision(2) << (param.max_kmer_freq * 100);
                                          return ss.str();
                                      })() + "%" :
                                      ">" + std::to_string((int)param.max_kmer_freq) + " occurrences") 
                    << ")" << std::endl;
        }

        std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
        std::cerr << "[wfmash::mashmap] reference index computed in " << timeRefSketch.count() << "s" << std::endl;

        if (this->minmerIndex.size() == 0)
        {
          std::cerr << "[wfmash::mashmap] ERROR, reference sketch is empty. "
                    << "Reference sequences shorter than the kmer size are not indexed" << std::endl;
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
      MI_Type* buildHelper(InputSeqContainer *input, progress_meter::ProgressMeter* progress)
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
                input->seqId,
                progress);

        return thread_output;
      }

      /**
       * @brief                 routine to handle thread's local minmer index
       * @param[in] output      thread local minmer output
       */

      void buildHandleThreadOutput(MI_Type* contigMinmerIndex)
      {
          // Count k-mer frequencies first
          HF_Map_t kmer_freqs;
          for (const auto& mi : *contigMinmerIndex) {
              kmer_freqs[mi.hash]++;
          }

          // This function is kept for compatibility but should not be used
          // when parallel index building is enabled
          for (MinmerInfo& mi : *contigMinmerIndex) {
              // Skip high-frequency k-mers
              auto freq_it = kmer_freqs.find(mi.hash);
              if (freq_it != kmer_freqs.end() && freq_it->second > param.max_kmer_freq) {
                  continue;
              }

              if (minmerPosLookupIndex[mi.hash].size() == 0 
                      || minmerPosLookupIndex[mi.hash].back().hash != mi.hash 
                      || minmerPosLookupIndex[mi.hash].back().pos != mi.wpos) {
                  minmerPosLookupIndex[mi.hash].push_back(IntervalPoint {mi.wpos, mi.hash, mi.seqId, side::OPEN});
                  minmerPosLookupIndex[mi.hash].push_back(IntervalPoint {mi.wpos_end, mi.hash, mi.seqId, side::CLOSE});
              } else {
                  minmerPosLookupIndex[mi.hash].back().pos = mi.wpos_end;
              }
          }

          // Only add k-mers that aren't too frequent
          MI_Type filtered_minmers;
          for (const auto& mi : *contigMinmerIndex) {
              auto freq_it = kmer_freqs.find(mi.hash);
              if (freq_it == kmer_freqs.end() || freq_it->second <= param.max_kmer_freq) {
                  filtered_minmers.push_back(mi);
              }
          }

          this->minmerIndex.insert(
              this->minmerIndex.end(), 
              std::make_move_iterator(filtered_minmers.begin()), 
              std::make_move_iterator(filtered_minmers.end()));

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
      // Removed writeFreqKmersBinary function


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
      void writeIndex(const std::vector<std::string>& target_subset, const std::string& filename = "", bool append = false, size_t batch_idx = 0, size_t total_batches = 1) 
      {
        fs::path indexFilename = filename.empty() ? fs::path(param.indexFilename) : fs::path(filename);
        std::ofstream outStream;
        if (append) {
            outStream.open(indexFilename, std::ios::binary | std::ios::app);
        } else {
            outStream.open(indexFilename, std::ios::binary);
        }
        if (!outStream) {
            std::cerr << "Error: Unable to open index file for writing: " << indexFilename << std::endl;
            exit(1);
        }
        writeSubIndexHeader(outStream, target_subset, batch_idx, total_batches);
        writeParameters(outStream);
        writeSketchBinary(outStream);
        writePosListBinary(outStream);
        // Removed writeFreqKmersBinary call
        outStream.close();
      }

      void writeSubIndexHeader(std::ofstream& outStream, const std::vector<std::string>& target_subset, size_t batch_idx = 0, size_t total_batches = 1) 
      {
        const uint64_t magic_number = 0xDEADBEEFCAFEBABE;
        outStream.write(reinterpret_cast<const char*>(&magic_number), sizeof(magic_number));
  
        // Write batch information
        outStream.write(reinterpret_cast<const char*>(&batch_idx), sizeof(batch_idx));
        outStream.write(reinterpret_cast<const char*>(&total_batches), sizeof(total_batches));
  
        uint64_t num_sequences = target_subset.size();
        outStream.write(reinterpret_cast<const char*>(&num_sequences), sizeof(num_sequences));
        for (const auto& seqName : target_subset) {
            uint64_t name_length = seqName.size();
            outStream.write(reinterpret_cast<const char*>(&name_length), sizeof(name_length));
            outStream.write(seqName.c_str(), name_length);
        }
        
        // Write sequence ID mappings to ensure consistency
        idManager.exportIdMapping(outStream);
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
       * @brief  Read parameters and compare to CLI params
       */
      void readParameters(std::ifstream& inStream)
      {
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
          std::cerr << "[wfmash::mashmap] ERROR: Parameters of indexed sketch differ from current parameters" << std::endl;
          std::cerr << "[wfmash::mashmap] Index --> segLength=" << index_segLength
                    << " sketchSize=" << index_sketchSize << " kmerSize=" << index_kmerSize << std::endl;
          std::cerr << "[wfmash::mashmap] Current --> segLength=" << param.segLength
                    << " sketchSize=" << param.sketchSize << " kmerSize=" << param.kmerSize << std::endl;
          exit(1);
        }
      }


      /**
       * @brief  Read all index data structures from file
       */
      void readIndex(std::ifstream& inStream, const std::vector<std::string>& targetSequenceNames) 
      {
        std::cerr << "[wfmash::mashmap] Reading index" << std::endl;
        size_t batch_idx, total_batches;
        if (!readSubIndexHeader(inStream, targetSequenceNames, batch_idx, total_batches)) {
            std::cerr << "Error: Sequences in the index do not match the expected target sequences." << std::endl;
            exit(1);
        }
        readParameters(inStream);
        readSketchBinary(inStream);
        readPosListBinary(inStream);
        // Removed readFreqKmersBinary call
      }

      bool readSubIndexHeader(std::ifstream& inStream, const std::vector<std::string>& targetSequenceNames, size_t& batch_idx, size_t& total_batches) 
      {
        uint64_t magic_number = 0;
        inStream.read(reinterpret_cast<char*>(&magic_number), sizeof(magic_number));
        if (magic_number != 0xDEADBEEFCAFEBABE) {
            std::cerr << "Error: Invalid magic number in index file." << std::endl;
            exit(1);
        }
  
        // Read batch information
        inStream.read(reinterpret_cast<char*>(&batch_idx), sizeof(batch_idx));
        inStream.read(reinterpret_cast<char*>(&total_batches), sizeof(total_batches));
  
        uint64_t num_sequences = 0;
        inStream.read(reinterpret_cast<char*>(&num_sequences), sizeof(num_sequences));
        
        if (num_sequences > 1000000) { // Sanity check
            std::cerr << "Error: Invalid number of sequences in index file: " << num_sequences << std::endl;
            exit(1);
        }
        
        std::vector<std::string> sequenceNames;
        for (uint64_t i = 0; i < num_sequences; ++i) {
            uint64_t name_length = 0;
            inStream.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            
            if (name_length > 10000) { // Sanity check
                std::cerr << "Error: Invalid sequence name length in index file: " << name_length << std::endl;
                exit(1);
            }
            
            std::string seqName(name_length, '\0');
            inStream.read(&seqName[0], name_length);
            sequenceNames.push_back(seqName);
        }
        
        // Read and restore sequence ID mappings from index
        idManager.importIdMapping(inStream);
        
        // Debug output of sequence ID mappings
        std::cerr << "[wfmash::mashmap] Index contains " << sequenceNames.size() << " sequences." << std::endl;
        
        // Check for sequence name matches
        bool all_found = true;
        std::vector<std::string> missing;
        for (const auto& seqName : targetSequenceNames) {
            try {
                idManager.getSequenceId(seqName);
            } catch (const std::runtime_error&) {
                missing.push_back(seqName);
                all_found = false;
            }
        }
        
        if (!missing.empty()) {
            std::cerr << "Warning: " << missing.size() << " sequence(s) not found in index:" << std::endl;
            for (size_t i = 0; i < std::min(missing.size(), size_t(5)); ++i) {
                std::cerr << "  - '" << missing[i] << "'" << std::endl;
            }
            if (missing.size() > 5) {
                std::cerr << "  - ... and " << (missing.size() - 5) << " more" << std::endl;
            }
            
            if (targetSequenceNames.size() == missing.size()) {
                std::cerr << "ERROR: None of the target sequences found in index!" << std::endl;
                
                // Print first few index sequences to help debugging
                std::cerr << "Index contains these sequences:" << std::endl;
                for (size_t i = 0; i < std::min(sequenceNames.size(), size_t(5)); ++i) {
                    std::cerr << "  - '" << sequenceNames[i] << "'" << std::endl;
                }
                if (sequenceNames.size() > 5) {
                    std::cerr << "  - ... and " << (sequenceNames.size() - 5) << " more" << std::endl;
                }
            } else {
                std::cerr << "These sequences will be skipped during mapping." << std::endl;
            }
        }
        
        // Allow proceeding even with missing sequences
        return true;
      }


      public:

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

      void clear()
      {
        minmerPosLookupIndex.clear();
        minmerIndex.clear();
        minmerFreqHistogram.clear();
      }

    }; //End of class Sketch
} //End of namespace skch

#endif
