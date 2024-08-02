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

seq_record_t* createSeqRecord(const MappingBoundaryRow& currentRecord, 
                              const std::string& mappingRecordLine,
                              faidx_t* ref_faidx,
                              faidx_t* query_faidx) {
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
                                      currentRecord.rStartPos - head_padding, 
                                      currentRecord.rEndPos + tail_padding, &ref_len);

    // Extract query sequence
    int64_t query_len;
    char* query_seq = faidx_fetch_seq64(query_faidx, currentRecord.qId.c_str(),
                                        currentRecord.qStartPos, currentRecord.qEndPos, &query_len);

    // Create a new seq_record_t object for the alignment
    seq_record_t* rec = new seq_record_t(currentRecord, mappingRecordLine,
                                         std::string(ref_seq, ref_len), 
                                         currentRecord.rStartPos - head_padding, ref_len, ref_size,
                                         std::string(query_seq, query_len), 
                                         currentRecord.qStartPos, query_len, query_size);

    // Clean up
    free(ref_seq);
    free(query_seq);

    return rec;
}

std::string processAlignment(seq_record_t* rec) {
    std::string& ref_seq = rec->refSequence;
    std::string& query_seq = rec->querySequence;

    skch::CommonFunc::makeUpperCaseAndValidDNA(ref_seq.data(), ref_seq.length());
    skch::CommonFunc::makeUpperCaseAndValidDNA(query_seq.data(), query_seq.length());

    // Adjust the reference sequence to start from the original start position
    char* ref_seq_ptr = &ref_seq[rec->currentRecord.rStartPos - rec->refStartPos];

    std::vector<char> queryRegionStrand(query_seq.size() + 1);

    if(rec->currentRecord.strand == skch::strnd::FWD) {
        std::copy(query_seq.begin(), query_seq.end(), queryRegionStrand.begin());
    } else {
        skch::CommonFunc::reverseComplement(query_seq.data(), queryRegionStrand.data(), query_seq.size());
    }

    wflign::wavefront::WFlign wflign(
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

    std::stringstream output;
    wflign.set_output(
        &output,
#ifdef WFA_PNG_TSV_TIMING
        !param.tsvOutputPrefix.empty(),
        nullptr,
        param.prefix_wavefront_plot_in_png,
        param.wfplot_max_size,
        !param.path_patching_info_in_tsv.empty(),
        nullptr,
#endif
        true, // merge alignments
        param.emit_md_tag,
        !param.sam_format,
        param.no_seq_in_sam);

    wflign.wflign_affine_wavefront(
        rec->currentRecord.qId,
        queryRegionStrand.data(),
        rec->queryTotalLength,
        rec->queryStartPos,
        rec->queryLen,
        rec->currentRecord.strand != skch::strnd::FWD,
        rec->currentRecord.refId,
        ref_seq_ptr,
        rec->refTotalLength,
        rec->currentRecord.rStartPos,
        rec->currentRecord.rEndPos - rec->currentRecord.rStartPos);

    return output.str();
}

void single_reader_thread(const std::string& input_file,
                          atomic_queue::AtomicQueue<std::string*, 1024>& line_queue,
                          std::atomic<bool>& reader_done) {
    std::ifstream mappingListStream(input_file);
    if (!mappingListStream.is_open()) {
        throw std::runtime_error("[wfmash::align::computeAlignments] Error! Failed to open input mapping file: " + input_file);
    }

    std::string line;
    while (std::getline(mappingListStream, line)) {
        if (!line.empty()) {
            std::string* line_ptr = new std::string(std::move(line));
            line_queue.push(line_ptr);
        }
    }

    mappingListStream.close();
    reader_done.store(true);
}

void processor_thread(std::atomic<size_t>& total_alignments_queued,
                      std::atomic<bool>& reader_done,
                      atomic_queue::AtomicQueue<std::string*, 1024>& line_queue,
                      seq_atomic_queue_t& seq_queue,
                      std::atomic<bool>& thread_should_exit) {
    faidx_t* local_ref_faidx = fai_load(param.refSequences.front().c_str());
    faidx_t* local_query_faidx = fai_load(param.querySequences.front().c_str());

    while (!thread_should_exit.load()) {
        std::string* line_ptr = nullptr;
        //std::cerr << "size of line queue " << line_queue.was_size() << std::endl;
        if (line_queue.try_pop(line_ptr)) {
            MappingBoundaryRow currentRecord;
            parseMashmapRow(*line_ptr, currentRecord);
            
            // Process the record and create seq_record_t
            seq_record_t* rec = createSeqRecord(currentRecord, *line_ptr, local_ref_faidx, local_query_faidx);
            //std::cerr << "size of seq_queue " << seq_queue.was_size() << std::endl;

            while (!seq_queue.try_push(rec)) {
                if (thread_should_exit.load()) {
                    delete rec;
                    delete line_ptr;
                    goto cleanup;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }

            ++total_alignments_queued;
            delete line_ptr;
        } else if (reader_done.load() && line_queue.was_empty()) {
            break;
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }

cleanup:
    fai_destroy(local_ref_faidx);
    fai_destroy(local_query_faidx);
}

void processor_manager(seq_atomic_queue_t& seq_queue,
                       atomic_queue::AtomicQueue<std::string*, 1024>& line_queue,
                       std::atomic<size_t>& total_alignments_queued,
                       std::atomic<bool>& reader_done,
                       std::atomic<bool>& processor_done,
                       size_t max_processors) {
    std::vector<std::thread> processor_threads;
    std::vector<std::atomic<bool>> thread_should_exit(max_processors);

    const size_t queue_capacity = seq_queue.capacity();
    const size_t low_threshold = queue_capacity * 0.2;
    const size_t high_threshold = queue_capacity * 0.8;

    auto spawn_processor = [&](size_t id) {
        //std::cerr << "spawn_processor: " << id << std::endl;
        thread_should_exit[id].store(false);
        processor_threads.emplace_back([this, &total_alignments_queued, &reader_done, &line_queue, &seq_queue, &thread_should_exit, id]() {
            this->processor_thread(total_alignments_queued, reader_done, line_queue, seq_queue, thread_should_exit[id]);
        });
    };

    // Start with one processor
    spawn_processor(0);
    size_t current_processors = 1;

    while (!reader_done.load() || !line_queue.was_empty() || !seq_queue.was_empty()) {
        size_t queue_size = seq_queue.was_size();

        //std::cerr << "queue_size: " << queue_size << std::endl;

        if (queue_size < low_threshold && current_processors < max_processors) {
            //std::cerr << "spawn_processor: " << current_processors << std::endl;
            spawn_processor(current_processors++);
        } else if (queue_size > high_threshold && current_processors > 1) {
            //std::cerr << "kill_processor: " << current_processors << std::endl;
            thread_should_exit[--current_processors].store(true);
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    // Signal all remaining threads to exit
    for (size_t i = 0; i < current_processors; ++i) {
        thread_should_exit[i].store(true);
    }

    // Wait for all processor threads to finish
    for (auto& thread : processor_threads) {
        thread.join();
    }
    processor_done.store(true);
}

void worker_thread(uint64_t tid,
                   std::atomic<bool>& is_working,
                   seq_atomic_queue_t& seq_queue,
                   paf_atomic_queue_t& paf_queue,
                   std::atomic<bool>& reader_done,
                   progress_meter::ProgressMeter& progress,
                   std::atomic<uint64_t>& processed_alignment_length) {
    is_working.store(true);
    while (true) {
        seq_record_t* rec = nullptr;
        if (seq_queue.try_pop(rec)) {
            is_working.store(true);
            std::string alignment_output = processAlignment(rec);
            
            // Push the alignment output to the paf_queue
            paf_queue.push(new std::string(std::move(alignment_output)));
            
            // Update progress meter and processed alignment length
            uint64_t alignment_length = rec->currentRecord.qEndPos - rec->currentRecord.qStartPos;
            progress.increment(alignment_length);
            processed_alignment_length.fetch_add(alignment_length, std::memory_order_relaxed);
            
            delete rec;
        } else if (reader_done.load() && seq_queue.was_empty()) {
            break;
        } else {
            is_working.store(false);
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }
    is_working.store(false);
}

void writer_thread(const std::string& output_file,
                   paf_atomic_queue_t& paf_queue,
                   std::atomic<bool>& reader_done,
                   std::atomic<bool>& processor_done,
                   const std::vector<std::atomic<bool>>& worker_working) {
    std::ofstream outstream(output_file);
    if (!outstream.is_open()) {
        throw std::runtime_error("[wfmash::align::computeAlignments] Error! Failed to open output file: " + output_file);
    }

    auto all_workers_done = [&]() {
        return std::all_of(worker_working.begin(), worker_working.end(),
                           [](const std::atomic<bool>& w) { return !w.load(); });
    };

    while (true) {
        std::string* paf_output = nullptr;
        if (paf_queue.try_pop(paf_output)) {
            outstream << *paf_output;
            outstream.flush();
            delete paf_output;
        } else if (reader_done.load() && processor_done.load() && paf_queue.was_empty() && all_workers_done()) {
            break;
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }

    outstream.close();
}

void computeAlignments() {
    std::atomic<size_t> total_alignments_queued(0);
    std::atomic<bool> reader_done(false);
    std::atomic<bool> processor_done(false);

    // Create queues
    atomic_queue::AtomicQueue<std::string*, 1024> line_queue;
    seq_atomic_queue_t seq_queue;
    paf_atomic_queue_t paf_queue;  // Add this line

    // Calculate max_processors based on the number of worker threads
    size_t max_processors = std::max(1UL, static_cast<unsigned long>(param.threads));

    // Calculate total alignment length
    uint64_t total_alignment_length = 0;
    {
        std::ifstream mappingListStream(param.mashmapPafFile);
        std::string mappingRecordLine;
        MappingBoundaryRow currentRecord;

        while(std::getline(mappingListStream, mappingRecordLine)) {
            if (!mappingRecordLine.empty()) {
                parseMashmapRow(mappingRecordLine, currentRecord);
                total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
            }
        }
    }

    // Create progress meter
    progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align::computeAlignments] aligned");

    // Create atomic counter for processed alignment length
    std::atomic<uint64_t> processed_alignment_length(0);

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Launch single reader thread
    std::thread single_reader([this, &line_queue, &reader_done]() {
        this->single_reader_thread(param.mashmapPafFile, line_queue, reader_done);
    });

    // Launch processor manager
    std::thread processor_manager_thread([this, &seq_queue, &line_queue, &total_alignments_queued, &reader_done, &processor_done, max_processors]() {
        this->processor_manager(seq_queue, line_queue, total_alignments_queued, reader_done, processor_done, max_processors);
    });

    // Launch worker threads
    std::vector<std::thread> workers;
    std::vector<std::atomic<bool>> worker_working(param.threads);
    for (uint64_t t = 0; t < param.threads; ++t) {
        workers.emplace_back([this, t, &worker_working, &seq_queue, &paf_queue, &reader_done, &progress, &processed_alignment_length]() {
            this->worker_thread(t, worker_working[t], seq_queue, paf_queue, reader_done, progress, processed_alignment_length);
        });
    }

    // Launch writer thread
    std::thread writer([this, &paf_queue, &reader_done, &processor_done, &worker_working]() {
        this->writer_thread(param.pafOutputFile, paf_queue, reader_done, processor_done, worker_working);
    });

    // Wait for all threads to complete
    single_reader.join();
    processor_manager_thread.join();
    for (auto& worker : workers) {
        worker.join();
    }
    writer.join();

    // Stop timing
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    // Finish progress meter
    progress.finish();

    std::cerr << "[wfmash::align::computeAlignments] "
              << "total aligned records = " << total_alignments_queued.load() 
              << ", total aligned bp = " << processed_alignment_length.load()
              << ", time taken = " << duration.count() << " seconds" << std::endl;
}
      
  };
}


#endif
