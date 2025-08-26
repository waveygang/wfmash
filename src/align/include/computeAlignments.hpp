/**
 * @file    computeAlignments.hpp
 * @brief   logic for generating alignments when given mashmap
 *          mappings as input
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMPUTE_ALIGNMENTS_HPP
#define COMPUTE_ALIGNMENTS_HPP

// #include <vector>
// #include <algorithm>
// #include <fstream>
// #include <sstream>
// #include <cassert>
// #include <thread>
// #include <memory>

// Include the reentrant FASTA/BGZF index implementation
#define REENTRANT_FAIDX_IMPLEMENTATION
#include "common/faigz.h"

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/wflign/src/wflign.hpp"
#include "common/seqiter.hpp"
// #include "common/progress.hpp"
#include "common/utils.hpp"
#include <any>
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/algorithm/partitioner.hpp>

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

      // Shared FASTA index metadata (thread-safe)
      faidx_meta_t* ref_meta;
      faidx_meta_t* query_meta;

    public:

      explicit Aligner(const align::Parameters &p) : param(p) {
          assert(param.refSequences.size() == 1);
          assert(param.querySequences.size() == 1);

          // Load index metadata (shared across threads)
          ref_meta = faidx_meta_load(param.refSequences.front().c_str(), FAI_FASTA, FAI_CREATE);
          if (!ref_meta) {
              throw std::runtime_error("[wfmash::align] Error! Failed to load reference FASTA index: "
                                      + param.refSequences.front());
          }

          query_meta = faidx_meta_load(param.querySequences.front().c_str(), FAI_FASTA, FAI_CREATE);
          if (!query_meta) {
              faidx_meta_destroy(ref_meta);
              throw std::runtime_error("[wfmash::align] Error! Failed to load query FASTA index: "
                                      + param.querySequences.front());
          }

          std::cerr << "[wfmash::align] Successfully loaded thread-safe FASTA indices" << std::endl;
      }

      ~Aligner() {
          if (ref_meta) faidx_meta_destroy(ref_meta);
          if (query_meta) faidx_meta_destroy(query_meta);
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
      inline static void parseMashmapRow(const std::string &mappingRecordLine, MappingBoundaryRow &currentRecord, const uint64_t target_padding, const uint64_t query_padding = 0) {
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

          // Parse chain info if present (expecting format "ch:Z:id.pos.len" in tokens[14])
          int64_t chain_id = -1;
          int64_t chain_length = 1;
          int64_t chain_pos = 1;
          if (tokens.size() > 14) {
              const auto chain_vec = split_view(std::string_view(tokens[14]), ':');
              if (chain_vec.size() == 3 && chain_vec[0] == "ch" && chain_vec[1] == "Z") {
                  // Split the id.pos.len format
                  const auto chain_parts = split_view(std::string_view(chain_vec[2]), '.');
                  if (chain_parts.size() == 3) {
			  try{
				  chain_id = std::stoll(std::string(chain_parts[0]));
				  chain_pos = std::stoll(std::string(chain_parts[1]));
				  chain_length = std::stoll(std::string(chain_parts[2]));
			  } catch (const std::exception& e) {
				  throw std::runtime_error("Failed to parse chain info '" +
						  std::string(tokens[14]) + "': " + e.what());
			  }
		  }
              }
          }

          // Save values into currentRecord
          {
              currentRecord.qId = std::string(tokens[0]);  // Need to copy ID strings
              currentRecord.qStartPos = std::stoll(std::string(tokens[2]));
              currentRecord.qEndPos = std::stoll(std::string(tokens[3]));
              currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
              currentRecord.refId = std::string(tokens[5]);  // Need to copy ID strings
              const uint64_t ref_len = std::stoull(std::string(tokens[6]));
              currentRecord.chain_id = chain_id;
              currentRecord.chain_length = chain_length;
              currentRecord.chain_pos = chain_pos;

              // Parse position values
              uint64_t rStartPos = std::stoll(std::string(tokens[7]));
              uint64_t rEndPos = std::stoll(std::string(tokens[8]));
              uint64_t qStartPos = currentRecord.qStartPos;
              uint64_t qEndPos = currentRecord.qEndPos;
              const uint64_t query_len = std::stoull(std::string(tokens[1]));  // Query sequence length


              // Apply target padding while ensuring we don't go below 0 or above reference length
              if (target_padding > 0) {
                // Always applied to reduce target holes
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

              // Apply query padding while ensuring we don't go below 0 or above query length
              if (query_padding > 0) {
                // Apply query padding only at the ends (left first piece and right last piece)
                // Do not pad the query in the middle to avoid overlaps between consecutive pieces of the same chain
                if (chain_pos == 1) {
                  if (qStartPos >= query_padding) {
                      qStartPos -= query_padding;
                  } else {
                      qStartPos = 0;
                  }
                }
                if (chain_pos == chain_length) {
                  if (qEndPos + query_padding <= query_len) {
                      qEndPos += query_padding;
                  } else {
                      qEndPos = query_len;
                  }

                  // Update the query positions
                  currentRecord.qStartPos = qStartPos;
                  currentRecord.qEndPos = qEndPos;
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

    // Store file offsets instead of full records
    struct MappingOffset {
        std::streampos offset;
        uint64_t query_length;
        uint64_t target_length;
    };
    std::vector<MappingOffset> mapping_offsets;
    uint64_t total_query_length = 0;
    uint64_t total_target_length = 0;
    
    {
        std::ifstream mappingListStream(param.mashmapPafFile);
        if (!mappingListStream.is_open()) {
            throw std::runtime_error("[wfmash::align] Error! Failed to open input mapping file: "
                                    + param.mashmapPafFile);
        }

        std::string mappingRecordLine;
        MappingBoundaryRow currentRecord;
        std::streampos current_pos = mappingListStream.tellg();

        while(std::getline(mappingListStream, mappingRecordLine)) {
            if (!mappingRecordLine.empty()) {
                try {
                    parseMashmapRow(mappingRecordLine, currentRecord, param.target_padding, param.query_padding);
                    uint64_t qlen = currentRecord.qEndPos - currentRecord.qStartPos;
                    uint64_t tlen = currentRecord.rEndPos - currentRecord.rStartPos;
                    
                    total_query_length += qlen;
                    total_target_length += tlen;
                    
                    // Store offset and lengths
                    mapping_offsets.push_back({current_pos, qlen, tlen});
                } catch (const std::exception& e) {
                    std::cerr << "[wfmash::align] Warning: Skipping invalid record: " << e.what() << std::endl;
                }
            }
            current_pos = mappingListStream.tellg();
        }
    }

    std::cerr << "[wfmash::align] Found " << mapping_offsets.size()
              << " mapping records for alignment ("
              << total_query_length << " query bp, "
              << total_target_length << " target bp)" << std::endl;

    // Progress meter
    auto progress = std::make_shared<progress_meter::ProgressMeter>(
        total_query_length,
        "[wfmash::align] aligning",
        param.use_progress_bar
        );

    // Create taskflow executor with thread count
    tf::Executor executor(param.threads);
    tf::Taskflow taskflow;

    // Mutex for synchronized output writing
    std::mutex output_mutex;

    // Using for_each with iterators and DynamicPartitioner for better load balancing
    taskflow.for_each(
        mapping_offsets.begin(),
        mapping_offsets.end(),
        [&](const MappingOffset& offset) {
            // Thread-local file handle for efficient reading
            thread_local std::ifstream tl_file;
            thread_local bool tl_file_opened = false;
            
            if (!tl_file_opened) {
                tl_file.open(param.mashmapPafFile);
                if (!tl_file.is_open()) {
                    std::cerr << "[wfmash::align] Error: Could not open mapping file for reading" << std::endl;
                    return;
                }
                tl_file_opened = true;
            }
            
            // Read the record from file
            std::string record;
            tl_file.seekg(offset.offset);
            std::getline(tl_file, record);
            
            if (!record.empty()) {
                processMappingRecord(
                    record,
                    ref_meta,
                    query_meta,
                    param,
                    total_alignments_processed,
                    processed_alignment_length,
                    progress,
                    output_mutex,
                    outstream
                );
            }
        },
        tf::DynamicPartitioner()
    );

    // Run the taskflow
    executor.run(taskflow).wait();

    // Close output stream
    outstream.close();

    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    // First finish progress meter to ensure its thread stops updating
    progress->finish();

    // Print summary after progress meter is fully stopped
    std::cerr << "[wfmash::align] "
              << "total aligned records = " << total_alignments_processed.load()
              << ", total aligned bp = " << processed_alignment_length.load()
              << ", completed in " << duration.count() << " seconds" << std::endl;
}

// Process a single mapping record (extracted to avoid lambda issues with for_each)
void processMappingRecord(
    const std::string& record,
    faidx_meta_t* ref_meta,
    faidx_meta_t* query_meta,
    const align::Parameters& param,
    std::atomic<uint64_t>& total_alignments_processed,
    std::atomic<uint64_t>& processed_alignment_length,
    std::shared_ptr<progress_meter::ProgressMeter>& progress,
    std::mutex& output_mutex,
    std::ofstream& outstream) {

    try {
        // Parse the mapping record
        MappingBoundaryRow currentRecord;
        parseMashmapRow(record, currentRecord, param.target_padding, param.query_padding);

        // Create sequence record
        std::unique_ptr<seq_record_t> seq_rec(
            createSeqRecord(currentRecord, record, ref_meta, query_meta)
        );

        // Process alignment
        std::string alignment_output = processAlignment(seq_rec.get());
        uint64_t alignment_length = currentRecord.qEndPos - currentRecord.qStartPos;

        // Format the output
        std::stringstream formatted_buffer;
        std::string_view output(alignment_output);
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

                    formatted_buffer << new_line;
                } else {
                    // If no CIGAR string found, output the line unchanged
                    formatted_buffer << line << '\n';
                }
            }

            start = end + 1;
            end = output.find('\n', start);
        }

        // Capture the formatted output
        std::string formatted_output = formatted_buffer.str();

        // Update statistics
        processed_alignment_length.fetch_add(alignment_length, std::memory_order_relaxed);
        total_alignments_processed.fetch_add(1, std::memory_order_relaxed);

        // Update progress
        progress->increment(alignment_length);

        // Write to output with minimal critical section
        if (!formatted_output.empty()) {
            std::lock_guard<std::mutex> lock(output_mutex);
            outstream << formatted_output;
            // Only flush occasionally to reduce I/O overhead
            if (total_alignments_processed.load(std::memory_order_relaxed) % 1000 == 0) {
                outstream.flush();
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "[wfmash::align] Error processing record: " << e.what() << std::endl;
    }
}

// Thread-local readers that persist for the lifetime of the thread
struct ThreadLocalReaders {
    faidx_reader_t* ref_reader;
    faidx_reader_t* query_reader;

    ThreadLocalReaders(faidx_meta_t* ref_meta, faidx_meta_t* query_meta)
        : ref_reader(nullptr), query_reader(nullptr) {
        // Create thread-local readers
        ref_reader = faidx_reader_create(ref_meta);
        if (!ref_reader) {
            throw std::runtime_error("Failed to create reference reader");
        }

        query_reader = faidx_reader_create(query_meta);
        if (!query_reader) {
            faidx_reader_destroy(ref_reader);
            throw std::runtime_error("Failed to create query reader");
        }
    }

    ~ThreadLocalReaders() {
        if (ref_reader) faidx_reader_destroy(ref_reader);
        if (query_reader) faidx_reader_destroy(query_reader);
    }
};

// Get thread-local readers from shared metadata
static ThreadLocalReaders& getThreadLocalReaders(faidx_meta_t* ref_meta, faidx_meta_t* query_meta) {
    thread_local ThreadLocalReaders readers(ref_meta, query_meta);
    return readers;
}

// Creates a sequence record using thread-local FASTA readers from the shared metadata
seq_record_t* createSeqRecord(const MappingBoundaryRow& currentRecord,
                              const std::string& mappingRecordLine,
                              faidx_meta_t* ref_meta,
                              faidx_meta_t* query_meta) {
    try {
        // Get thread-local readers (created once per thread)
        ThreadLocalReaders& readers = getThreadLocalReaders(ref_meta, query_meta);
        faidx_reader_t* ref_reader = readers.ref_reader;
        faidx_reader_t* query_reader = readers.query_reader;

        // Get the reference sequence length
        const hts_pos_t ref_size = faidx_meta_seq_len(ref_meta, currentRecord.refId.c_str());
        if (ref_size < 0) {
            faidx_reader_destroy(ref_reader);
            faidx_reader_destroy(query_reader);
            throw std::runtime_error("Reference sequence not found: " + currentRecord.refId);
        }

        // Get the query sequence length
        const hts_pos_t query_size = faidx_meta_seq_len(query_meta, currentRecord.qId.c_str());
        if (query_size < 0) {
            faidx_reader_destroy(ref_reader);
            faidx_reader_destroy(query_reader);
            throw std::runtime_error("Query sequence not found: " + currentRecord.qId);
        }

        // Compute padding for sequence extraction
        const uint64_t head_padding = currentRecord.rStartPos >= param.wflign_max_len_minor
            ? param.wflign_max_len_minor : currentRecord.rStartPos;
        const uint64_t tail_padding = ref_size - currentRecord.rEndPos >= param.wflign_max_len_minor
            ? param.wflign_max_len_minor : ref_size - currentRecord.rEndPos;

        // Extract reference sequence
        hts_pos_t ref_len;
        char* ref_seq = faidx_reader_fetch_seq(
            ref_reader,
            currentRecord.refId.c_str(),
            currentRecord.rStartPos - head_padding,
            currentRecord.rEndPos + tail_padding - 1, // faigz uses inclusive end
            &ref_len);

        if (!ref_seq || ref_len <= 0) {
            faidx_reader_destroy(ref_reader);
            faidx_reader_destroy(query_reader);
            throw std::runtime_error("Failed to fetch reference sequence");
        }

        // Extract query sequence
        hts_pos_t query_len;
        char* query_seq = faidx_reader_fetch_seq(
            query_reader,
            currentRecord.qId.c_str(),
            currentRecord.qStartPos,
            currentRecord.qEndPos - 1, // faigz uses inclusive end
            &query_len);

        if (!query_seq || query_len <= 0) {
            free(ref_seq);
            faidx_reader_destroy(ref_reader);
            faidx_reader_destroy(query_reader);
            throw std::runtime_error("Failed to fetch query sequence");
        }

        // Create a new seq_record_t object that takes ownership of the sequences
        seq_record_t* rec = new seq_record_t(currentRecord, mappingRecordLine,
                                          ref_seq, // Transfer ownership
                                          currentRecord.rStartPos - head_padding, ref_len, ref_size,
                                          query_seq, // Transfer ownership
                                          currentRecord.qStartPos, query_len, query_size);

        // Note: We no longer destroy readers here as they are thread-local and will be cleaned up automatically

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
        param.disable_chain_patching,
        param.min_identity,
        param.min_alignment_length,
        param.min_block_identity,
        param.wflign_max_len_minor,
        rec->currentRecord.mashmap_estimated_identity,
        rec->currentRecord.chain_id,
        rec->currentRecord.chain_length,
        rec->currentRecord.chain_pos);

    return output.str();
}

void write_sam_header(std::ofstream& outstream) {
    // Use the FASTA metadata to get sequence names and lengths
    int num_seqs = faidx_meta_nseq(ref_meta);
    for (int i = 0; i < num_seqs; i++) {
        const char* seq_name = faidx_meta_iseq(ref_meta, i);
        if (seq_name) {
            hts_pos_t seq_len = faidx_meta_seq_len(ref_meta, seq_name);
            outstream << "@SQ\tSN:" << seq_name << "\tLN:" << seq_len << "\n";
        }
    }
    outstream << "@PG\tID:wfmash\tPN:wfmash\tVN:" << WFMASH_GIT_VERSION << "\tCL:wfmash\n";
}

  };
}


#endif
