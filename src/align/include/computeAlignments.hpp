/**
 * @file    computeAlignments.hpp
 * @brief   logic for generating alignments when given mashmap 
 *          mappings as input
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMPUTE_ALIGNMENTS_HPP 
#define COMPUTE_ALIGNMENTS_HPP

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <cassert>
#include <thread>
#include <memory>
#include <htslib/faidx.h>

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/wflign/src/wflign.hpp"
#include "common/atomic_queue/atomic_queue.h"
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "common/utils.hpp"

namespace align
{

long double float2phred(long double prob) {
    if (prob == 1)
        return 255;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
}

struct seq_record_t {
    MappingBoundaryRow currentRecord;
    std::string mappingRecordLine;
    std::string refSequence;  
    std::string querySequence;
    uint64_t refStartPos;
    uint64_t refLen;
    uint64_t refTotalLength;
    uint64_t queryStartPos;
    uint64_t queryLen;
    uint64_t queryTotalLength;
    seq_record_t(const MappingBoundaryRow& c, const std::string& r, 
                 const std::string& ref, uint64_t refStart, uint64_t refLength, uint64_t refTotalLength,
                 const std::string& query, uint64_t queryStart, uint64_t queryLength, uint64_t queryTotalLength)
        : currentRecord(c)
        , mappingRecordLine(r)
        , refSequence(ref)
        , querySequence(query)
        , refStartPos(refStart)
        , refLen(refLength)
        , refTotalLength(refTotalLength)
        , queryStartPos(queryStart)
        , queryLen(queryLength)
        , queryTotalLength(queryTotalLength)
        { }
};

/**
 * @brief A single-producer, multi-consumer (SPMC) atomic queue for storing pointers to seq_record_t objects.
 *
 * This queue is designed for a setup where there is a single producer and multiple consumers.
 * The producer enqueues pointers to seq_record_t objects, which represent sequences to be processed.
 * Multiple consumers dequeue these pointers and perform long-running alignment processes on the sequences.
 *
 * The queue has the following characteristics:
 * - Capacity: 1024 elements
 * - Default value for empty elements: nullptr
 * - MINIMIZE_CONTENTION: true (minimizes contention among consumers)
 * - MAXIMIZE_THROUGHPUT: true (optimized for high throughput)
 * - TOTAL_ORDER: false (relaxed memory ordering for better performance)
 * - SPSC: false (single-producer, multi-consumer mode)
 */
typedef atomic_queue::AtomicQueue<seq_record_t*, 1024, nullptr, true, true, false, false> seq_atomic_queue_t;

/**
 * @brief A multi-producer, single-consumer (MPSC) atomic queue for storing pointers to std::string objects.
 *
 * This queue is designed for a setup where there are multiple producers and a single consumer.
 * Multiple producers enqueue pointers to std::string objects, which represent PAF (Pairwise Alignment Format) strings.
 * The single consumer dequeues these pointers and writes out the PAF strings.
 *
 * The queue has the following characteristics:
 * - Capacity: 1024 elements
 * - Default value for empty elements: nullptr
 * - MINIMIZE_CONTENTION: true (minimizes contention among producers)
 * - MAXIMIZE_THROUGHPUT: true (optimized for high throughput)
 * - TOTAL_ORDER: false (relaxed memory ordering for better performance)
 * - SPSC: false (multi-producer, single-consumer mode)
 */
typedef atomic_queue::AtomicQueue<std::string*, 1024, nullptr, true, true, false, false> paf_atomic_queue_t;


  /**
   * @class     align::Aligner
   * @brief     compute alignments and generate sam output
   *            from mashmap mappings
   */
  class Aligner
  {
    private:

      //algorithm parameters
      const align::Parameters &param;

      faidx_t* ref_faidx;
      faidx_t* query_faidx;

    public:

      explicit Aligner(const align::Parameters &p) : param(p) {
          assert(param.refSequences.size() == 1);
          assert(param.querySequences.size() == 1);
          ref_faidx = fai_load(param.refSequences.front().c_str());
          query_faidx = fai_load(param.querySequences.front().c_str());
      }

      ~Aligner() {
          fai_destroy(ref_faidx);
          fai_destroy(query_faidx);  
      }
      
      /**
       * @brief                 compute alignments
       */
      void compute()
      {
        this->computeAlignments();
      }

      /**
       * @brief       parse mashmap row sequence
       * @param[in]   mappingRecordLine
       * @param[out]  currentRecord
       */
      inline static void parseMashmapRow(const std::string &mappingRecordLine, MappingBoundaryRow &currentRecord) {
          std::stringstream ss(mappingRecordLine); // Insert the string into a stream
          std::string word; // Have a buffer string

          vector<std::string> tokens; // Create vector to hold our words

          while (ss >> word) {
              tokens.push_back(word);
          }

          // Check if the number of tokens is at least 13
          if (tokens.size() < 13) {
              throw std::runtime_error("[wfmash::align::parseMashmapRow] Error! Invalid mashmap mapping record: " + mappingRecordLine);
          }

          // Extract the mashmap identity from the string
          const vector<string> mm_id_vec = skch::CommonFunc::split(tokens[12], ':');
          // if the estimated identity is missing, avoid assuming too low values
          const float mm_id = wfmash::is_a_number(mm_id_vec.back()) ? std::stof(mm_id_vec.back()) : skch::fixed::percentage_identity;

          // Save words into currentRecord
          {
              currentRecord.qId = tokens[0];
              currentRecord.qStartPos = std::stoi(tokens[2]);
              currentRecord.qEndPos = std::stoi(tokens[3]);
              currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
              currentRecord.refId = tokens[5];
              currentRecord.rStartPos = std::stoi(tokens[7]);
              currentRecord.rEndPos = std::stoi(tokens[8]);
              currentRecord.mashmap_estimated_identity = mm_id;
          }
      }

  private:

      /**
       * @brief                 parse query sequences and mashmap mappings
       *                        to compute sequence alignments
       */
      void computeAlignments()
      {
          uint64_t total_seqs = 0;

          // Count the number of mapped bases to align
          uint64_t total_alignment_length = 0;
          {
              std::ifstream mappingListStream(param.mashmapPafFile);
              std::string mappingRecordLine;
              MappingBoundaryRow currentRecord;
              uint64_t total_records = 0;

              while(!mappingListStream.eof()) {
                  std::getline(mappingListStream, mappingRecordLine);
                  if (!mappingRecordLine.empty()) {
                      parseMashmapRow(mappingRecordLine, currentRecord);
                      total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
                      ++total_records;
                  }
              }
          }

          progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align::computeAlignments] aligned");

          // input atomic queue
          seq_atomic_queue_t seq_queue;
          // output atomic queues
          paf_atomic_queue_t paf_queue, tsv_queue, patching_tsv_queue;
          // flag when we're done reading
          std::atomic<bool> reader_done;
          reader_done.store(false);

          auto& nthreads = param.threads;
          //for (

          // atomics to record if we're working or not
          std::vector<std::atomic<bool>> working(nthreads);
          for (auto& w : working) {
              w.store(true);
          }

          size_t total_alignments_queued = 0;
          auto reader_thread = [&]() {
              std::ifstream mappingListStream(param.mashmapPafFile);
              if (!mappingListStream.is_open()) {
                  throw std::runtime_error("[wfmash::align::computeAlignments] Error! Failed to open input mapping file: " + param.mashmapPafFile);
              }

              std::string mappingRecordLine;
              MappingBoundaryRow currentRecord;

              while (!mappingListStream.eof()) {
                  std::getline(mappingListStream, mappingRecordLine);
                  if (!mappingRecordLine.empty()) {
                      parseMashmapRow(mappingRecordLine, currentRecord);
            
                      // Get the reference sequence length
                      const int64_t ref_size = faidx_seq_len(ref_faidx, currentRecord.refId.c_str());
                      // Get the query sequence length
                      const int64_t query_size = faidx_seq_len(query_faidx, currentRecord.qId.c_str());

                      // Compute padding
                      const uint64_t head_padding = currentRecord.rStartPos >= param.wflign_max_len_minor
                          ? param.wflign_max_len_minor : currentRecord.rStartPos;
                      const uint64_t tail_padding = ref_size - currentRecord.rEndPos >= param.wflign_max_len_minor
                          ? param.wflign_max_len_minor : ref_size - currentRecord.rEndPos;

                      // Extract reference sequence
                      int64_t ref_len;
                      char* ref_seq = faidx_fetch_seq64(ref_faidx, currentRecord.refId.c_str(),
                                                        currentRecord.rStartPos - head_padding, currentRecord.rEndPos + tail_padding, &ref_len);

                      // Extract query sequence
                      int64_t query_len;
                      char* query_seq = faidx_fetch_seq64(query_faidx, currentRecord.qId.c_str(),
                                                          currentRecord.qStartPos, currentRecord.qEndPos, &query_len);

                      // Create a new seq_record_t object for the alignment using std::move
                      seq_record_t* rec = new seq_record_t(currentRecord, mappingRecordLine,
                                                           std::string(ref_seq, ref_len), currentRecord.rStartPos - head_padding, ref_len, ref_size,
                                                           std::string(query_seq, query_len), currentRecord.qStartPos, query_len, query_size);

                      // Clean up
                      free(ref_seq);
                      free(query_seq);

                      ++total_alignments_queued;
                      seq_queue.push(rec);
                  }
              }

              mappingListStream.close();
              reader_done.store(true);
          };

          // helper to check if we're still aligning
          auto still_working =
              [&](const std::vector<std::atomic<bool>>& working) {
                  bool ongoing = false;
                  for (auto& w : working) {
                      ongoing = ongoing || w.load();
                  }
                  return ongoing;
              };

          // writer, picks output from queue and writes it to our output stream
          std::ofstream outstrm(param.pafOutputFile, ios::app);

          size_t total_alignments_written = 0;
          auto writer_thread =
              [&]() {
                  while (true) {
                      std::string* paf_lines = nullptr;
                      if (!paf_queue.try_pop(paf_lines)
                          && !still_working(working)) {
                          break;
                      } else if (paf_lines != nullptr) {
                          ++total_alignments_written;
                          outstrm << *paf_lines;
                          delete paf_lines;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
              };

#ifdef WFA_PNG_TSV_TIMING
          auto writer_thread_tsv =
                  [&]() {
              if (!param.tsvOutputPrefix.empty()) {
                  uint64_t num_alignments_completed = 0;

                  while (true) {
                      std::string* tsv_lines = nullptr;
                      if (!tsv_queue.try_pop(tsv_lines)
                      && !still_working(working)) {
                          break;
                      } else if (tsv_lines != nullptr) {
                          std::ofstream ofstream_tsv(param.tsvOutputPrefix + std::to_string(num_alignments_completed++) + ".tsv");
                          ofstream_tsv << *tsv_lines;
                          ofstream_tsv.close();

                          delete tsv_lines;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
              }
          };

          std::ofstream ofstream_patching_tsv(param.path_patching_info_in_tsv);
          auto writer_thread_patching_tsv =
                  [&]() {
                      if (!param.path_patching_info_in_tsv.empty()) {
                          while (true) {
                              std::string* tsv_lines = nullptr;
                              if (!patching_tsv_queue.try_pop(tsv_lines)
                                  && !still_working(working)) {
                                  break;
                              } else if (tsv_lines != nullptr) {
                                  ofstream_patching_tsv << *tsv_lines;

                                  delete tsv_lines;
                              } else {
                                  std::this_thread::sleep_for(100ns);
                              }
                          }
                      }
                  };
#endif

          // worker, takes candidate alignments and runs wfa alignment on them
          auto worker_thread = 
              [&](uint64_t tid,
                  std::atomic<bool>& is_working) {
                  is_working.store(true);
                  while (true) {
                      seq_record_t* rec = nullptr;
                      if (!seq_queue.try_pop(rec)
                          && reader_done.load()) {
                          break;
                      } else if (rec != nullptr) {
                          std::stringstream output;
#ifdef WFA_PNG_TSV_TIMING
                          std::stringstream output_tsv;
                          std::stringstream patching_output_tsv;
#endif
                          doAlignment(
                              output,
#ifdef WFA_PNG_TSV_TIMING
                              output_tsv,
                              patching_output_tsv,
#endif
                              rec,
                              tid);
                          progress.increment(rec->currentRecord.qEndPos - rec->currentRecord.qStartPos);

                          auto* paf_rec = new std::string(output.str());
                          if (!paf_rec->empty()) {
                              paf_queue.push(paf_rec);
                          } else {
                              delete paf_rec;
                          }

#ifdef WFA_PNG_TSV_TIMING
                          auto* tsv_rec = new std::string(output_tsv.str());
                          if (!tsv_rec->empty()) {
                              tsv_queue.push(tsv_rec);
                          } else {
                              delete tsv_rec;
                          }

                          auto* patching_tsv_rec = new std::string(patching_output_tsv.str());
                          if (!patching_tsv_rec->empty()) {
                              patching_tsv_queue.push(patching_tsv_rec);
                          } else {
                              delete patching_tsv_rec;
                          }
#endif

                          delete rec;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
                  is_working.store(false);
              };

          // launch reader
          std::thread reader(reader_thread);
          // launch PAF/SAM writer
          std::thread writer(writer_thread);
#ifdef WFA_PNG_TSV_TIMING
          // launch TSV writer
          std::thread writer_tsv(writer_thread_tsv);
          std::thread writer_patching_tsv(writer_thread_patching_tsv);
#endif
          // launch workers
          std::vector<std::thread> workers; workers.reserve(nthreads);
          for (uint64_t t = 0; t < nthreads; ++t) {
              workers.emplace_back(worker_thread,
                                   t,
                                   std::ref(working[t]));
          }

          // wait for reader and workers to complete
          reader.join();
          for (auto& worker : workers) {
              worker.join();
          }
          // and finally the writer
          writer.join();
#ifdef WFA_PNG_TSV_TIMING
          writer_tsv.join();
          writer_patching_tsv.join();
          ofstream_patching_tsv.close();
#endif

          progress.finish();
          std::cerr << "[wfmash::align::computeAlignments] "
                    << "count of mapped reads = " << total_seqs
                    << ", total aligned bp = " << total_alignment_length << std::endl;
      }

      // core alignment computation function
      void doAlignment(
          std::stringstream& output,
#ifdef WFA_PNG_TSV_TIMING
          std::stringstream& output_tsv,
          std::stringstream& patching_output_tsv,
#endif
          seq_record_t* rec,
          uint64_t tid) {
#ifdef DEBUG
          std::cerr << "INFO, align::Aligner::doAlignment, aligning mashmap record: " << rec->mappingRecordLine << std::endl;
#endif

          std::string& ref_seq = rec->refSequence;
          std::string& query_seq = rec->querySequence;

          skch::CommonFunc::makeUpperCaseAndValidDNA(ref_seq.data(), ref_seq.length());
          skch::CommonFunc::makeUpperCaseAndValidDNA(query_seq.data(), query_seq.length());

          // Adjust the reference sequence to start from the original start position
          char* ref_seq_ptr = &ref_seq[rec->currentRecord.rStartPos - rec->refStartPos];

          char* queryRegionStrand = new char[query_seq.size() + 1];

          if(rec->currentRecord.strand == skch::strnd::FWD) {
              strncpy(queryRegionStrand, query_seq.data(), query_seq.size());
          } else {
              skch::CommonFunc::reverseComplement(query_seq.data(), queryRegionStrand, query_seq.size());
          }

          // To distinguish split alignment in SAM output format (currentRecord.rankMapping == 0 to avoid the suffix there is just one alignment for the query)
          //const std::string query_name_suffix = param.split && param.sam_format ? "_" + std::to_string(rec->currentRecord.rankMapping) : "";

          wflign::wavefront::WFlign* wflign = new wflign::wavefront::WFlign(
              param.wflambda_segment_length,
              param.min_identity,
              param.force_biwfa_alignment,
              param.wfa_mismatch_score,
              param.wfa_gap_opening_score,
              param.wfa_gap_extension_score,
              param.wfa_patching_mismatch_score,
              param.wfa_patching_gap_opening_score1,
              param.wfa_patching_gap_extension_score1,
              param.wfa_patching_gap_opening_score2,
              param.wfa_patching_gap_extension_score2,
              rec->currentRecord.mashmap_estimated_identity,
              param.wflign_mismatch_score,
              param.wflign_gap_opening_score,
              param.wflign_gap_extension_score,
              param.wflign_max_mash_dist,
              param.wflign_min_wavefront_length,
              param.wflign_max_distance_threshold,
              param.wflign_max_len_major,
              param.wflign_max_len_minor,
              param.wflign_erode_k,
              param.chain_gap,
              param.wflign_min_inv_patch_len,
              param.wflign_max_patching_score);
          wflign->set_output(
              &output,
#ifdef WFA_PNG_TSV_TIMING
              !param.tsvOutputPrefix.empty(),
              &output_tsv,
              param.prefix_wavefront_plot_in_png,
              param.wfplot_max_size,
              !param.path_patching_info_in_tsv.empty(),
              &patching_output_tsv,
#endif
              true, // merge alignments
              param.emit_md_tag,
              !param.sam_format,
              param.no_seq_in_sam);
          wflign->wflign_affine_wavefront(
              rec->currentRecord.qId,// + query_name_suffix,
              queryRegionStrand,
              rec->queryTotalLength,
              rec->queryStartPos,
              rec->queryLen,
              rec->currentRecord.strand != skch::strnd::FWD,
              rec->currentRecord.refId,
              ref_seq_ptr,
              rec->refTotalLength,
              rec->currentRecord.rStartPos,
              rec->currentRecord.rEndPos - rec->currentRecord.rStartPos);
          delete wflign;

          delete[] queryRegionStrand;
          // n.b. rec is deleted in calling context
      }
  };
}


#endif
