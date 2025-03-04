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
                // (using main class faidx with optimized BGZF queue size of threads/4)
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
            
            std::string_view output(alignment_outputs[pf.line()]);
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

// Creates a sequence record using the provided faidx objects
// When called from taskflow pipeline, this uses the optimized thread pool with BGZF queue size of threads/4
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

    // Create a new seq_record_t object that takes ownership of the sequences
    seq_record_t* rec = new seq_record_t(currentRecord, mappingRecordLine,
                                         ref_seq, // Transfer ownership of ref_seq
                                         currentRecord.rStartPos - head_padding, ref_len, ref_size,
                                         query_seq, // Transfer ownership of query_seq
                                         currentRecord.qStartPos, query_len, query_size);

    // No free() calls here - the seq_record_t destructor will handle cleanup

    return rec;
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
      
  };
}


#endif
