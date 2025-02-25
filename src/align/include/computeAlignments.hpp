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
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>
#include <htslib/hfile.h>

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
#include <any>
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/pipeline.hpp>

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
      hts_tpool* thread_pool;

    public:

      explicit Aligner(const align::Parameters &p) : param(p) {
          assert(param.refSequences.size() == 1);
          assert(param.querySequences.size() == 1);
          
          // Initialize thread pool
          thread_pool = hts_tpool_init(param.threads);
          if (!thread_pool) {
              throw std::runtime_error("[wfmash::align] Error! Failed to create thread pool");
          }
          
          // Load faidx indices
          ref_faidx = fai_load(param.refSequences.front().c_str());
          if (!ref_faidx) {
              hts_tpool_destroy(thread_pool);
              throw std::runtime_error("[wfmash::align] Error! Failed to load reference FASTA index: " + param.refSequences.front());
          }
          
          query_faidx = fai_load(param.querySequences.front().c_str());
          if (!query_faidx) {
              fai_destroy(ref_faidx);
              hts_tpool_destroy(thread_pool);
              throw std::runtime_error("[wfmash::align] Error! Failed to load query FASTA index: " + param.querySequences.front());
          }
          
          // Attach thread pool to faidx objects using the proper API
          // Calculate optimal queue size based on thread count (1/4 of threads)
          int bgzf_queue_size = std::max(1, param.threads / 4);
          
          if (fai_thread_pool(ref_faidx, thread_pool, bgzf_queue_size) != 0) {
              std::cerr << "[wfmash::align] Warning: Failed to attach thread pool to reference FASTA index" << std::endl;
          }
          
          if (fai_thread_pool(query_faidx, thread_pool, bgzf_queue_size) != 0) {
              std::cerr << "[wfmash::align] Warning: Failed to attach thread pool to query FASTA index" << std::endl;
          }
      }

      ~Aligner() {
          fai_destroy(ref_faidx);
          fai_destroy(query_faidx);
          hts_tpool_destroy(thread_pool);
      }
      
      /**
       * @brief                 compute alignments
       */

      void compute()
      {
        this->computeAlignmentsTaskflow();
      }

      /**
       * @brief       parse mashmap row sequence
       * @param[in]   mappingRecordLine
       * @param[out]  currentRecord
       */
      inline static void parseMashmapRow(const std::string &mappingRecordLine, MappingBoundaryRow &currentRecord, const uint64_t target_padding) {
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

          // Parse chain info if present (expecting format "chain:i:id.pos.len" in tokens[14])
          int32_t chain_id = -1;
          int32_t chain_length = 1;
          int32_t chain_pos = 1;
          if (tokens.size() > 14) {
              const vector<string> chain_vec = skch::CommonFunc::split(tokens[14], ':');
              if (chain_vec.size() == 3 && chain_vec[0] == "chain" && chain_vec[1] == "i") {
                  // Split the id.pos.len format
                  const vector<string> chain_parts = skch::CommonFunc::split(chain_vec[2], '.');
                  if (chain_parts.size() == 3) {
                      chain_id = std::stoi(chain_parts[0]);
                      chain_pos = std::stoi(chain_parts[1]); 
                      chain_length = std::stoi(chain_parts[2]);
                  }
              }
          }

          // Save words into currentRecord
          {
              currentRecord.qId = tokens[0];
              currentRecord.qStartPos = std::stoi(tokens[2]);
              currentRecord.qEndPos = std::stoi(tokens[3]);
              currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
              currentRecord.refId = tokens[5];
              const uint64_t ref_len = std::stoi(tokens[6]);
              currentRecord.chain_id = chain_id;
              currentRecord.chain_length = chain_length;
              currentRecord.chain_pos = chain_pos;
              
              // Apply target padding while ensuring we don't go below 0 or above reference length
              uint64_t rStartPos = std::stoi(tokens[7]);
              uint64_t rEndPos = std::stoi(tokens[8]);
              
              // Always apply target padding
              if (target_padding > 0) {
                  if (rStartPos >= target_padding) {
                      rStartPos -= target_padding;
                  } else {
                      rStartPos = 0;
                  }
                  if (rEndPos + target_padding <= ref_len) {
                      rEndPos += target_padding;
                  } else {
                      rEndPos = ref_len;
                  }
              }

              // Validate coordinates against reference length
              if (rStartPos >= ref_len || rEndPos > ref_len) {
                  std::cerr << "[parse-debug] ERROR: Coordinates exceed reference length!" << std::endl;
                  throw std::runtime_error("[wfmash::align::parseMashmapRow] Error! Coordinates exceed reference length: " 
                                         + std::to_string(rStartPos) + "-" + std::to_string(rEndPos) 
                                         + " (ref_len=" + std::to_string(ref_len) + ")");
              }
              
              currentRecord.rStartPos = rStartPos;
              currentRecord.rEndPos = rEndPos;
              currentRecord.mashmap_estimated_identity = mm_id;
          }
      }

  private:

void computeAlignmentsTaskflow() {
    // Prepare output file
    std::ofstream outstream(param.pafOutputFile);
    if (!outstream.is_open()) {
        throw std::runtime_error("[wfmash::align] Error! Failed to open output file: " + param.pafOutputFile);
    }
    
    // Write SAM header if needed
    if (param.sam_format) {
        write_sam_header(outstream);
    }
    
    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Create concurrent queues for communication between stages
    const size_t queue_capacity = 65536; // Increased 8x from 8192 to buffer more records (2^16)
    atomic_queue::AtomicQueue<std::string*, queue_capacity, nullptr> line_queue;
    std::atomic<bool> reader_done(false);
    std::atomic<uint64_t> total_alignments_processed(0);
    std::atomic<uint64_t> processed_alignment_length(0);
    
    // Create a taskflow to manage the reading task
    tf::Taskflow reader_taskflow;
    tf::Executor reader_executor(1); // Single thread for reading
    
    // Calculate total alignment length for progress tracking
    uint64_t total_alignment_length = 0;
    {
        std::ifstream mappingListStream(param.mashmapPafFile);
        std::string mappingRecordLine;
        MappingBoundaryRow currentRecord;

        while(std::getline(mappingListStream, mappingRecordLine)) {
            if (!mappingRecordLine.empty()) {
                parseMashmapRow(mappingRecordLine, currentRecord, param.target_padding);
                total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
            }
        }
    }
    
    // Progress meter
    progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align] aligned");
    
    // Create a task to read input in a separate thread
    auto read_task = reader_taskflow.emplace([&]() {
        std::ifstream mappingStream(param.mashmapPafFile);
        if (!mappingStream.is_open()) {
            throw std::runtime_error("[wfmash::align] Error! Failed to open input mapping file: " + param.mashmapPafFile);
        }
        
        // Use batch reading to reduce overhead
        const size_t batch_size = 100;
        std::vector<std::string> batch_lines;
        batch_lines.reserve(batch_size);
        
        std::string line;
        while(true) {
            // Read a batch of lines
            batch_lines.clear();
            for (size_t i = 0; i < batch_size && std::getline(mappingStream, line); ++i) {
                if (!line.empty()) {
                    batch_lines.push_back(std::move(line));
                }
                line.clear(); // Ensure line is empty for next read
            }
            
            // If no lines were read, we're done
            if (batch_lines.empty()) {
                break;
            }
            
            // Process the batch
            std::vector<std::string*> line_ptrs;
            line_ptrs.reserve(batch_lines.size());
            
            // Create heap strings
            for (auto& bline : batch_lines) {
                line_ptrs.push_back(new std::string(std::move(bline)));
            }
            
            // Push batch to queue with adaptive backoff
            size_t pushed = 0;
            size_t backoff_time = 1; // Start with 1ms
            
            while (pushed < line_ptrs.size()) {
                bool made_progress = false;
                
                // Try to push as many as possible
                while (pushed < line_ptrs.size() && line_queue.try_push(line_ptrs[pushed])) {
                    ++pushed;
                    made_progress = true;
                    backoff_time = 1; // Reset backoff on progress
                }
                
                // If we couldn't push everything, wait with adaptive backoff
                if (pushed < line_ptrs.size()) {
                    if (!made_progress) {
                        std::this_thread::sleep_for(std::chrono::milliseconds(backoff_time));
                        backoff_time = std::min(backoff_time * 2, size_t(100)); // Exponential backoff, max 100ms
                    }
                }
            }
        }
        
        reader_done.store(true);
    }).name("read_input_file");
    
    // Launch reader task asynchronously
    reader_executor.run(reader_taskflow);
    
    // Create storage for pipeline data
    const size_t num_pipeline_lines = std::max<size_t>(256, param.threads * 64);
    std::vector<seq_record_t*> records(num_pipeline_lines, nullptr);
    std::vector<std::string> alignment_outputs(num_pipeline_lines);
    
    // Create a taskflow for the pipeline
    tf::Taskflow taskflow;
    tf::Executor executor(param.threads);
    
    // Define our pipeline
    tf::Pipeline pipeline(num_pipeline_lines,
        // Stage 1: Parse mapping record (SERIAL)
        tf::Pipe{tf::PipeType::SERIAL, [&, this](tf::Pipeflow& pf) {
            // Try to get a line from the queue with retry logic
            std::string* line_ptr = nullptr;
            
            // Try a few times with short intervals before yielding
            for (int retry = 0; retry < 3; ++retry) {
                if (line_queue.try_pop(line_ptr)) {
                    break;
                }
                // Very short pause between immediate retries
                std::this_thread::sleep_for(std::chrono::microseconds(10));
            }
            
            if (line_ptr) {
                // Parse the current record
                MappingBoundaryRow currentRecord;
                parseMashmapRow(*line_ptr, currentRecord, param.target_padding);
                
                // Create seq_record_t using the thread-pool enabled faidx
                records[pf.line()] = createSeqRecord(currentRecord, *line_ptr, 
                                                    this->ref_faidx, this->query_faidx);
                
                // Clean up the line
                delete line_ptr;
            } else if (reader_done.load() && line_queue.was_empty()) {
                // No more input to process
                pf.stop();
                return;
            } else {
                // Queue is temporarily empty but reader is still working
                // Use a progressive backoff strategy rather than just yielding
                static thread_local int empty_count = 0;
                if (++empty_count > 10) {
                    // If we've had many empty attempts, sleep for a bit longer
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    empty_count = 0;
                } else {
                    // Otherwise just yield to give other tasks a chance
                    std::this_thread::yield();
                }
            }
        }},
        
        // Stage 2: Perform alignment (PARALLEL)
        tf::Pipe{tf::PipeType::PARALLEL, [&, this](tf::Pipeflow& pf) {
            // Process alignment for this record
            seq_record_t* rec = records[pf.line()];
            if (rec == nullptr) return;
            
            alignment_outputs[pf.line()] = processAlignment(rec);
            
            // Update progress
            uint64_t alignment_length = rec->currentRecord.qEndPos - rec->currentRecord.qStartPos;
            progress.increment(alignment_length);
            processed_alignment_length.fetch_add(alignment_length, std::memory_order_relaxed);
            total_alignments_processed.fetch_add(1, std::memory_order_relaxed);
        }},
        
        // Stage 3: Write results (SERIAL)
        tf::Pipe{tf::PipeType::SERIAL, [&](tf::Pipeflow& pf) {
            // Write the output for this record
            if (records[pf.line()] == nullptr) return;
            if (alignment_outputs[pf.line()].empty()) return;
            
            std::stringstream ss(alignment_outputs[pf.line()]);
            std::string line;
            while (std::getline(ss, line)) {
                if (line.empty()) continue;
                
                std::vector<std::string> fields;
                std::stringstream field_ss(line);
                std::string field;
                while (field_ss >> field) {
                    fields.push_back(field);
                }

                // Find the CIGAR string field (should be after cg:Z:)
                auto cigar_it = std::find_if(fields.begin(), fields.end(),
                    [](const std::string& s) { 
                        return s.size() >= 5 && s.substr(0, 5) == "cg:Z:"; 
                    });
                
                if (cigar_it != fields.end()) {
                    // Reconstruct the line
                    std::string new_line;
                    for (const auto& f : fields) {
                        if (!new_line.empty()) new_line += '\t';
                        new_line += f;
                    }
                    new_line += '\n';
                    
                    outstream << new_line;
                } else {
                    // If no CIGAR string found, output the line unchanged
                    outstream << line << '\n';
                }
            }
            
            outstream.flush();
            
            // Clean up
            delete records[pf.line()];
            records[pf.line()] = nullptr;
            alignment_outputs[pf.line()].clear();
        }}
    );
    
    // Build the pipeline task
    tf::Task task = taskflow.composed_of(pipeline).name("alignment_pipeline");
    
    // Run the pipeline and wait for completion
    executor.run(taskflow).wait();
    
    // Wait for reader to finish (should already be done)
    reader_executor.wait_for_all();
    
    // Clean up any remaining items in the line queue
    std::string* remaining_line = nullptr;
    while (line_queue.try_pop(remaining_line)) {
        delete remaining_line;
    }
    
    // Close output stream
    outstream.close();
    
    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    // Finish progress meter
    progress.finish();
    
    std::cerr << "[wfmash::align] "
              << "total aligned records = " << total_alignments_processed.load()
              << ", total aligned bp = " << processed_alignment_length.load()
              << ", time taken = " << duration.count() << " seconds" << std::endl;
}

seq_record_t* createSeqRecord(const MappingBoundaryRow& currentRecord, 
                              const std::string& mappingRecordLine,
                              faidx_t* ref_faidx,
                              faidx_t* query_faidx) {
    // Get the reference sequence length
    const int64_t ref_size = faidx_seq_len(ref_faidx, currentRecord.refId.c_str());
    // Get the query sequence length
    const int64_t query_size = faidx_seq_len(query_faidx, currentRecord.qId.c_str());

    // Compute padding for sequence extraction
    const uint64_t head_padding = currentRecord.rStartPos >= param.wflign_max_len_minor
        ? param.wflign_max_len_minor : currentRecord.rStartPos;
    const uint64_t tail_padding = ref_size - currentRecord.rEndPos >= param.wflign_max_len_minor
        ? param.wflign_max_len_minor : ref_size - currentRecord.rEndPos;

    // Extract reference sequence
    int64_t ref_len;
    char* ref_seq = faidx_fetch_seq64(ref_faidx, currentRecord.refId.c_str(),
                                      currentRecord.rStartPos - head_padding, 
                                      currentRecord.rEndPos - 1 + tail_padding, &ref_len);

    // Extract query sequence
    int64_t query_len;
    char* query_seq = faidx_fetch_seq64(query_faidx, currentRecord.qId.c_str(),
                                        currentRecord.qStartPos, currentRecord.qEndPos - 1, &query_len);

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

    // Set up penalties for biWFA
    wflign_penalties_t wfa_penalties;
    wfa_penalties.match = 0;
    wfa_penalties.mismatch = param.wfa_patching_mismatch_score;
    wfa_penalties.gap_opening1 = param.wfa_patching_gap_opening_score1;
    wfa_penalties.gap_extension1 = param.wfa_patching_gap_extension_score1;
    wfa_penalties.gap_opening2 = param.wfa_patching_gap_opening_score2;
    wfa_penalties.gap_extension2 = param.wfa_patching_gap_extension_score2;

    std::stringstream output;

    // Do direct biWFA alignment
    wflign::wavefront::do_biwfa_alignment(
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
        rec->currentRecord.rEndPos - rec->currentRecord.rStartPos,
        output,
        wfa_penalties,
        param.emit_md_tag,
        !param.sam_format,
        param.no_seq_in_sam,
        param.min_identity,
        param.wflign_max_len_minor,
        rec->currentRecord.mashmap_estimated_identity,
        rec->currentRecord.chain_id,
        rec->currentRecord.chain_length,
        rec->currentRecord.chain_pos);

    return output.str();
}

void single_reader_thread(const std::string& input_file,
                          atomic_queue::AtomicQueue<std::string*, 1024>& line_queue,
                          std::atomic<bool>& reader_done) {
    std::ifstream mappingListStream(input_file);
    if (!mappingListStream.is_open()) {
        throw std::runtime_error("[wfmash::align] Error! Failed to open input mapping file: " + input_file);
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
    // Create a local thread pool for this thread
    hts_tpool* local_thread_pool = hts_tpool_init(1); // Single thread pool for this thread
    if (!local_thread_pool) {
        std::cerr << "[wfmash::align] Warning: Failed to create local thread pool in processor thread" << std::endl;
    }
    
    // Load FASTA indices
    faidx_t* local_ref_faidx = fai_load(param.refSequences.front().c_str());
    faidx_t* local_query_faidx = fai_load(param.querySequences.front().c_str());
    
    // Attach thread pool to faidx objects if available
    if (local_thread_pool) {
        // Single-threaded pool, so use queue size of 1
        if (local_ref_faidx && fai_thread_pool(local_ref_faidx, local_thread_pool, 1) != 0) {
            std::cerr << "[wfmash::align] Warning: Failed to attach local thread pool to reference FASTA index" << std::endl;
        }
        
        if (local_query_faidx && fai_thread_pool(local_query_faidx, local_thread_pool, 1) != 0) {
            std::cerr << "[wfmash::align] Warning: Failed to attach local thread pool to query FASTA index" << std::endl;
        }
    }

    while (!thread_should_exit.load()) {
        std::string* line_ptr = nullptr;
        if (line_queue.try_pop(line_ptr)) {
            MappingBoundaryRow currentRecord;
            parseMashmapRow(*line_ptr, currentRecord, param.target_padding);
            
            // Process the record and create seq_record_t
            seq_record_t* rec = createSeqRecord(currentRecord, *line_ptr, local_ref_faidx, local_query_faidx);

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
    hts_tpool_destroy(local_thread_pool);
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
    const size_t low_threshold = 1;
    const size_t high_threshold = queue_capacity * 0.8;

    auto spawn_processor = [&](size_t id) {
        thread_should_exit[id].store(false);
        processor_threads.emplace_back([this, &total_alignments_queued, &reader_done, &line_queue, &seq_queue, &thread_should_exit, id]() {
            this->processor_thread(total_alignments_queued, reader_done, line_queue, seq_queue, thread_should_exit[id]);
        });
    };

    // Start with one processor
    spawn_processor(0);
    size_t current_processors = 1;
    uint64_t exhausted = 0;

    while (!reader_done.load() || !line_queue.was_empty() || !seq_queue.was_empty()) {
        size_t queue_size = seq_queue.was_size();

        if (param.multithread_fasta_input) {
            if (queue_size < low_threshold && current_processors < max_processors) {
                ++exhausted;
            } else if (queue_size > high_threshold && current_processors > 1) {
                thread_should_exit[--current_processors].store(true);
            }

            if (exhausted > 20 && queue_size < low_threshold) {
                spawn_processor(current_processors++);
                exhausted = 0;
            }
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
                   std::atomic<bool>& processor_done,
                   progress_meter::ProgressMeter& progress,
                   std::atomic<uint64_t>& processed_alignment_length) {
    is_working.store(true);
    while (true) {
        seq_record_t* rec = nullptr;
        if (seq_queue.try_pop(rec)) {
            is_working.store(true);
            std::string alignment_output = processAlignment(rec);

            // Parse the alignment output to find CIGAR string and coordinates
            std::stringstream ss(alignment_output);
            std::string line;
            while (std::getline(ss, line)) {
                if (line.empty()) continue;
                
                std::vector<std::string> fields;
                std::stringstream field_ss(line);
                std::string field;
                while (field_ss >> field) {
                    fields.push_back(field);
                }

                // Find the CIGAR string field (should be after cg:Z:)
                auto cigar_it = std::find_if(fields.begin(), fields.end(),
                    [](const std::string& s) { return s.substr(0, 5) == "cg:Z:"; });
                
                if (cigar_it != fields.end()) {
                    std::string cigar = cigar_it->substr(5); // Remove cg:Z: prefix
                    uint64_t ref_start = std::stoull(fields[7]);
                    uint64_t ref_end = std::stoull(fields[8]);


                    // Just pass through the CIGAR string and coordinates unchanged
                    // The trimming is now handled in wflign namespace

                    // Reconstruct the line
                    std::string new_line;
                    for (const auto& f : fields) {
                        if (!new_line.empty()) new_line += '\t';
                        new_line += f;
                    }
                    new_line += '\n';
                    
                    // Push the modified alignment output to the paf_queue
                    paf_queue.push(new std::string(std::move(new_line)));
                } else {
                    // If no CIGAR string found, output the line unchanged
                    paf_queue.push(new std::string(line + '\n'));
                }
            }
            
            // Update progress meter and processed alignment length
            uint64_t alignment_length = rec->currentRecord.qEndPos - rec->currentRecord.qStartPos;
            progress.increment(alignment_length);
            processed_alignment_length.fetch_add(alignment_length, std::memory_order_relaxed);
            
            delete rec;
        } else if (reader_done.load() && processor_done.load() && seq_queue.was_empty()) {
            break;
        } else {
            is_working.store(false);
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }
    is_working.store(false);
}

void write_sam_header(std::ofstream& outstream) {
    for(const auto &fileName : param.refSequences) {
        // check if there is a .fai
        std::string fai_name = fileName + ".fai";
        if (fs::exists(fai_name)) {
            // if so, process the .fai to determine our sequence length
            std::string line;
            std::ifstream in(fai_name.c_str());
            while (std::getline(in, line)) {
                auto line_split = skch::CommonFunc::split(line, '\t');
                const std::string seq_name = line_split[0];
                const uint64_t seq_len = std::stoull(line_split[1]);
                outstream << "@SQ\tSN:" << seq_name << "\tLN:" << seq_len << "\n";
            }
        } else {
            // if not, warn that this is expensive
            std::cerr << "[wfmash::align] WARNING, no .fai index found for " << fileName << ", reading the file to prepare SAM header (slow)" << std::endl;
            seqiter::for_each_seq_in_file(
                fileName, {}, "",
                [&](const std::string& seq_name,
                    const std::string& seq) {
                    outstream << "@SQ\tSN:" << seq_name << "\tLN:" << seq.length() << "\n";
                });
        }
    }
    outstream << "@PG\tID:wfmash\tPN:wfmash\tVN:" << WFMASH_GIT_VERSION << "\tCL:wfmash\n";
}

void writer_thread(const std::string& output_file,
                   paf_atomic_queue_t& paf_queue,
                   std::atomic<bool>& reader_done,
                   std::atomic<bool>& processor_done,
                   const std::vector<std::atomic<bool>>& worker_working) {
    std::ofstream outstream(output_file);
    // if the output file is SAM, we write the header
    if (param.sam_format) {
        write_sam_header(outstream);
    }

    if (!outstream.is_open()) {
        throw std::runtime_error("[wfmash::align] Error! Failed to open output file: " + output_file);
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
                parseMashmapRow(mappingRecordLine, currentRecord, param.target_padding);
                total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
            }
        }
    }

    // Create progress meter
    progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align] aligned");

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
        workers.emplace_back([this, t, &worker_working, &seq_queue, &paf_queue, &reader_done, &processor_done, &progress, &processed_alignment_length]() {
            this->worker_thread(t, worker_working[t], seq_queue, paf_queue, reader_done, processor_done, progress, processed_alignment_length);
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

    std::cerr << "[wfmash::align] "
              << "total aligned records = " << total_alignments_queued.load() 
              << ", total aligned bp = " << processed_alignment_length.load()
              << ", time taken = " << duration.count() << " seconds" << std::endl;
}
      
  };
}


#endif
