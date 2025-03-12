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
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>
#include <htslib/hfile.h>
#include "common/thread_safe_faidx.hpp"

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/wflign/src/wflign.hpp"
#include "common/seqiter.hpp"
#include "common/progress.hpp"
#include "common/utils.hpp"
#include <any>
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/pipeline.hpp>

namespace align
{

// String_view based version of split
inline std::vector<std::string_view> split_view(std::string_view s, char delim) {
    std::vector<std::string_view> result;
    size_t pos = 0;
    size_t found;
    
    while ((found = s.find(delim, pos)) != std::string_view::npos) {
        result.emplace_back(s.substr(pos, found - pos));
        pos = found + 1;
    }
    
    // Add the last part
    if (pos <= s.size()) {
        result.emplace_back(s.substr(pos));
    }
    
    return result;
}

// String_view based version of whitespace tokenization
inline std::vector<std::string_view> tokenize_view(std::string_view s) {
    std::vector<std::string_view> tokens;
    size_t pos = 0;
    size_t start;
    
    while (pos < s.size()) {
        // Skip whitespace
        while (pos < s.size() && std::isspace(s[pos])) pos++;
        if (pos >= s.size()) break;
        
        // Find token end
        start = pos;
        while (pos < s.size() && !std::isspace(s[pos])) pos++;
        
        // Create string_view (zero copy)
        tokens.emplace_back(s.substr(start, pos - start));
    }
    
    return tokens;
}

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
    
    // Raw pointers that own the memory
    char* raw_ref_sequence = nullptr;
    char* raw_query_sequence = nullptr;
    
    // Non-owning views of the sequences
    std::string_view refSequence;
    std::string_view querySequence;
    
    uint64_t refStartPos;
    uint64_t refLen;
    uint64_t refTotalLength;
    uint64_t queryStartPos;
    uint64_t queryLen;
    uint64_t queryTotalLength;

    seq_record_t(const MappingBoundaryRow& c, const std::string& r, 
                 char* ref_seq, uint64_t refStart, uint64_t refLength, uint64_t refTotalLength,
                 char* query_seq, uint64_t queryStart, uint64_t queryLength, uint64_t queryTotalLength)
        : currentRecord(c)
        , mappingRecordLine(r)
        , raw_ref_sequence(ref_seq)
        , raw_query_sequence(query_seq)
        , refSequence(ref_seq, refLength)
        , querySequence(query_seq, queryLength)
        , refStartPos(refStart)
        , refLen(refLength)
        , refTotalLength(refTotalLength)
        , queryStartPos(queryStart)
        , queryLen(queryLength)
        , queryTotalLength(queryTotalLength)
        { }
    
    ~seq_record_t() {
        free(raw_ref_sequence);
        free(raw_query_sequence);
    }
};



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

      ts_faidx::FastaReader* ref_faidx;
      ts_faidx::FastaReader* query_faidx;

    public:

      explicit Aligner(const align::Parameters &p) : param(p) {
          assert(param.refSequences.size() == 1);
          assert(param.querySequences.size() == 1);
          
          // Load thread-safe FASTA readers
          try {
              ref_faidx = new ts_faidx::FastaReader(param.refSequences.front());
          } catch (const std::exception& e) {
              throw std::runtime_error("[wfmash::align] Error! Failed to load reference FASTA index: " 
                                      + param.refSequences.front() + " - " + e.what());
          }
          
          try {
              query_faidx = new ts_faidx::FastaReader(param.querySequences.front());
          } catch (const std::exception& e) {
              delete ref_faidx;
              throw std::runtime_error("[wfmash::align] Error! Failed to load query FASTA index: " 
                                      + param.querySequences.front() + " - " + e.what());
          }
          
          std::cerr << "[wfmash::align] Successfully loaded thread-safe FASTA indices" << std::endl;
      }

      ~Aligner() {
          delete ref_faidx;
          delete query_faidx;
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
          auto tokens = tokenize_view(mappingRecordLine);
          
          // Check if the number of tokens is at least 13
          if (tokens.size() < 13) {
              throw std::runtime_error("[wfmash::align::parseMashmapRow] Error! Invalid mashmap mapping record: " + mappingRecordLine);
          }

          // Extract the mashmap identity from the string
          const auto mm_id_vec = split_view(std::string_view(tokens[12]), ':');
          // if the estimated identity is missing, avoid assuming too low values
          const float mm_id = !mm_id_vec.empty() && wfmash::is_a_number(std::string(mm_id_vec.back())) 
                              ? std::stof(std::string(mm_id_vec.back())) 
                              : skch::fixed::percentage_identity;

          // Parse chain info if present (expecting format "chain:i:id.pos.len" in tokens[14])
          int32_t chain_id = -1;
          int32_t chain_length = 1;
          int32_t chain_pos = 1;
          if (tokens.size() > 14) {
              const auto chain_vec = split_view(std::string_view(tokens[14]), ':');
              if (chain_vec.size() == 3 && chain_vec[0] == "chain" && chain_vec[1] == "i") {
                  // Split the id.pos.len format
                  const auto chain_parts = split_view(std::string_view(chain_vec[2]), '.');
                  if (chain_parts.size() == 3) {
                      chain_id = std::stoi(std::string(chain_parts[0]));
                      chain_pos = std::stoi(std::string(chain_parts[1])); 
                      chain_length = std::stoi(std::string(chain_parts[2]));
                  }
              }
          }

          // Save values into currentRecord
          {
              currentRecord.qId = std::string(tokens[0]);  // Need to copy ID strings
              currentRecord.qStartPos = std::stoi(std::string(tokens[2]));
              currentRecord.qEndPos = std::stoi(std::string(tokens[3]));
              currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
              currentRecord.refId = std::string(tokens[5]);  // Need to copy ID strings
              const uint64_t ref_len = std::stoull(std::string(tokens[6]));
              currentRecord.chain_id = chain_id;
              currentRecord.chain_length = chain_length;
              currentRecord.chain_pos = chain_pos;
              
              // Apply target padding while ensuring we don't go below 0 or above reference length
              uint64_t rStartPos = std::stoi(std::string(tokens[7]));
              uint64_t rEndPos = std::stoi(std::string(tokens[8]));
              
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

// Structure to hold alignment results and metadata for transmission between tasks
struct alignment_result_t {
    std::string output;
    uint64_t alignment_length;
    bool success = false;

    alignment_result_t() = default;
    alignment_result_t(std::string&& out, uint64_t len) 
        : output(std::move(out)), alignment_length(len), success(true) {}
};

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
    
    // Tracking variables for statistics
    std::atomic<uint64_t> total_alignments_processed(0);
    std::atomic<uint64_t> processed_alignment_length(0);
    
    // Calculate total alignment length for progress tracking
    uint64_t total_alignment_length = 0;
    size_t total_records = 0;
    {
        std::ifstream mappingListStream(param.mashmapPafFile);
        std::string mappingRecordLine;
        MappingBoundaryRow currentRecord;

        while(std::getline(mappingListStream, mappingRecordLine)) {
            if (!mappingRecordLine.empty()) {
                parseMashmapRow(mappingRecordLine, currentRecord, param.target_padding);
                total_alignment_length += currentRecord.qEndPos - currentRecord.qStartPos;
                total_records++;
            }
        }
    }
    
    // Progress meter
    progress_meter::ProgressMeter progress(total_alignment_length, "[wfmash::align] aligned");
    
    // Create three taskflows with separate executors
    tf::Taskflow reader_taskflow;
    tf::Taskflow worker_taskflow;
    tf::Taskflow writer_taskflow;
    
    // Executors with appropriate thread counts
    tf::Executor reader_executor(1);  // Single reader thread
    tf::Executor worker_executor(param.threads);  // Worker threads for alignment
    tf::Executor writer_executor(1);  // Single writer thread
    
    // Create thread-safe queues for communication between tasks
    using MappingQueue = tf::BoundedTaskQueue<std::string>;
    using ResultQueue = tf::BoundedTaskQueue<alignment_result_t>;
    
    // Create queues with reasonable capacity for our workload
    constexpr size_t queue_capacity = 1024;
    auto mapping_queue = std::make_shared<MappingQueue>(queue_capacity);
    auto result_queue = std::make_shared<ResultQueue>(queue_capacity);
    
    // Signal flags for completion
    std::atomic<bool> reader_done(false);
    std::atomic<bool> all_workers_done(false);
    
    // Create reader task - reads mapping records and pushes to mapping_queue
    auto reader_task = reader_taskflow.emplace([&]() {
        std::ifstream mappingStream(param.mashmapPafFile);
        if (!mappingStream.is_open()) {
            throw std::runtime_error("[wfmash::align] Error! Failed to open input mapping file: " + param.mashmapPafFile);
        }
        
        std::string line;
        size_t count = 0;
        while(std::getline(mappingStream, line)) {
            if (!line.empty()) {
                mapping_queue->push(std::move(line));
                count++;
            }
        }
        
        std::cerr << "[wfmash::align] Reader task completed, pushed " << count << " records" << std::endl;
        reader_done.store(true);
    });
    
    // Create worker tasks - take records from mapping_queue, process, push results to result_queue
    const int num_workers = param.threads;
    std::vector<tf::Task> worker_tasks;
    
    for (int i = 0; i < num_workers; i++) {
        worker_tasks.push_back(worker_taskflow.emplace([&, this, i]() {
            std::string mapping_record;
            size_t processed = 0;
            
            while (true) {
                // Try to get work from the queue
                bool got_work = mapping_queue->pop(mapping_record);
                
                if (!got_work) {
                    // If reader is done and queue is empty, we're done
                    if (reader_done.load() && mapping_queue->empty()) {
                        break;
                    }
                    // Otherwise, yield and try again
                    std::this_thread::yield();
                    continue;
                }
                
                try {
                    // Process the mapping record
                    MappingBoundaryRow currentRecord;
                    parseMashmapRow(mapping_record, currentRecord, param.target_padding);
                    
                    // Create sequence record by fetching from thread-safe readers
                    seq_record_t* rec = createSeqRecord(currentRecord, mapping_record, 
                                                   this->ref_faidx, this->query_faidx);
                    
                    // Process the alignment
                    std::string alignment_output = processAlignment(rec);
                    
                    // Push result to output queue
                    uint64_t alignment_length = currentRecord.qEndPos - currentRecord.qStartPos;
                    if (!alignment_output.empty()) {
                        result_queue->push(alignment_result_t(std::move(alignment_output), alignment_length));
                    }
                    
                    // Update progress and stats
                    progress.increment(alignment_length);
                    processed_alignment_length.fetch_add(alignment_length, std::memory_order_relaxed);
                    total_alignments_processed.fetch_add(1, std::memory_order_relaxed);
                    
                    // Clean up
                    delete rec;
                    processed++;
                } catch (const std::exception& e) {
                    std::cerr << "[wfmash::align] Worker " << i << " error: " << e.what() << std::endl;
                }
            }
            
            std::cerr << "[wfmash::align] Worker " << i << " processed " << processed << " records" << std::endl;
        }));
    }
    
    // Create writer task - takes results from result_queue and writes to output file
    auto writer_task = writer_taskflow.emplace([&]() {
        alignment_result_t result;
        size_t written = 0;
        
        while (true) {
            // Try to get a result from the queue
            bool got_result = result_queue->pop(result);
            
            if (!got_result) {
                // If all workers are done and queue is empty, we're done
                if (all_workers_done.load() && result_queue->empty()) {
                    break;
                }
                // Otherwise, yield and try again
                std::this_thread::yield();
                continue;
            }
            
            if (result.success) {
                // Process the alignment output
                std::string_view output(result.output);
                size_t start = 0;
                size_t end = output.find('\n');
                
                while (end != std::string_view::npos) {
                    // Process one line at a time without copying
                    std::string_view line = output.substr(start, end - start);
                    if (!line.empty()) {
                        // Parse fields using string_view
                        auto fields = tokenize_view(line);
                        
                        // Find the CIGAR string field (should be after cg:Z:)
                        auto cigar_it = std::find_if(fields.begin(), fields.end(),
                            [](std::string_view s) { 
                                return s.size() >= 5 && s.substr(0, 5) == "cg:Z:"; 
                            });
                        
                        if (cigar_it != fields.end()) {
                            // Reconstruct the line efficiently with a single allocation
                            std::string new_line;
                            new_line.reserve(line.size() + 1); // +1 for newline
                            
                            for (const auto& f : fields) {
                                if (!new_line.empty()) new_line += '\t';
                                new_line.append(f.data(), f.size());
                            }
                            new_line += '\n';
                            
                            outstream << new_line;
                        } else {
                            // If no CIGAR string found, output the line unchanged
                            outstream << line << '\n';
                        }
                    }
                    
                    start = end + 1;
                    end = output.find('\n', start);
                }
                
                outstream.flush();
                written++;
            }
        }
        
        std::cerr << "[wfmash::align] Writer task completed, wrote " << written << " records" << std::endl;
    });
    
    // Start the reader task
    reader_executor.run(reader_taskflow);
    
    // Start the worker tasks
    worker_executor.run(worker_taskflow);
    
    // Start the writer task
    writer_executor.run(writer_taskflow);
    
    // Wait for reader and worker tasks to complete
    reader_executor.wait_for_all();
    worker_executor.wait_for_all();
    
    // Signal that all workers are done
    all_workers_done.store(true);
    
    // Wait for writer task to complete
    writer_executor.wait_for_all();
    
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

// Creates a sequence record using the thread-safe FASTA readers
seq_record_t* createSeqRecord(const MappingBoundaryRow& currentRecord, 
                              const std::string& mappingRecordLine,
                              ts_faidx::FastaReader* ref_faidx,
                              ts_faidx::FastaReader* query_faidx) {
    try {
        // Get the reference sequence length
        const int64_t ref_size = ref_faidx->get_sequence_length(currentRecord.refId);
        if (ref_size < 0) {
            throw std::runtime_error("Reference sequence not found: " + currentRecord.refId);
        }
        
        // Get the query sequence length
        const int64_t query_size = query_faidx->get_sequence_length(currentRecord.qId);
        if (query_size < 0) {
            throw std::runtime_error("Query sequence not found: " + currentRecord.qId);
        }

        // Compute padding for sequence extraction
        const uint64_t head_padding = currentRecord.rStartPos >= param.wflign_max_len_minor
            ? param.wflign_max_len_minor : currentRecord.rStartPos;
        const uint64_t tail_padding = ref_size - currentRecord.rEndPos >= param.wflign_max_len_minor
            ? param.wflign_max_len_minor : ref_size - currentRecord.rEndPos;

        // Extract reference sequence (thread-safe)
        std::string ref_seq_str = ref_faidx->fetch_sequence(
            currentRecord.refId,
            currentRecord.rStartPos - head_padding,
            currentRecord.rEndPos + tail_padding);
        int64_t ref_len = ref_seq_str.length();
        char* ref_seq = strdup(ref_seq_str.c_str());

        // Extract query sequence (thread-safe)
        std::string query_seq_str = query_faidx->fetch_sequence(
            currentRecord.qId,
            currentRecord.qStartPos,
            currentRecord.qEndPos);
        int64_t query_len = query_seq_str.length();
        char* query_seq = strdup(query_seq_str.c_str());

        // Create a new seq_record_t object that takes ownership of the sequences
        seq_record_t* rec = new seq_record_t(currentRecord, mappingRecordLine,
                                           ref_seq, // Transfer ownership of ref_seq
                                           currentRecord.rStartPos - head_padding, ref_len, ref_size,
                                           query_seq, // Transfer ownership of query_seq
                                           currentRecord.qStartPos, query_len, query_size);

        return rec;
    } catch (const std::exception& e) {
        std::cerr << "[wfmash::align] Error extracting sequence: " << e.what() << std::endl;
        throw;
    }
}

std::string processAlignment(seq_record_t* rec) {
    // Thread-local buffer for query strand
    thread_local std::vector<char> queryRegionStrand;
    
    // Resize buffer only if needed
    if (queryRegionStrand.size() < rec->querySequence.size() + 1) {
        queryRegionStrand.resize(rec->querySequence.size() + 1);
    }
    
    // Make sequences uppercase and valid DNA
    skch::CommonFunc::makeUpperCaseAndValidDNA(rec->raw_ref_sequence, rec->refLen);
    skch::CommonFunc::makeUpperCaseAndValidDNA(rec->raw_query_sequence, rec->queryLen);

    // Adjust the reference sequence to start from the original start position
    char* ref_seq_ptr = &rec->raw_ref_sequence[rec->currentRecord.rStartPos - rec->refStartPos];

    if(rec->currentRecord.strand == skch::strnd::FWD) {
        std::copy(rec->querySequence.begin(), rec->querySequence.end(), queryRegionStrand.begin());
    } else {
        skch::CommonFunc::reverseComplement(rec->raw_query_sequence, queryRegionStrand.data(), rec->queryLen);
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

void write_sam_header(std::ofstream& outstream) {
    // Use our thread-safe FastaReader to get sequence names and lengths
    for (const auto& seq_name : ref_faidx->get_sequence_names()) {
        int64_t seq_len = ref_faidx->get_sequence_length(seq_name);
        outstream << "@SQ\tSN:" << seq_name << "\tLN:" << seq_len << "\n";
    }
    outstream << "@PG\tID:wfmash\tPN:wfmash\tVN:" << WFMASH_GIT_VERSION << "\tCL:wfmash\n";
}
      
  };
}


#endif
