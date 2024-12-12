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
#include <iomanip>
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

bool left_align_leading_deletion(
    const std::string& reference,
    const std::string& query,
    int& match_len,
    int del_len,
    int max_shifts)
{
    if (match_len == 0) {
        return true;
    }

    int shifts = 0;
    int total_possible_shifts = std::min({match_len, max_shifts, 
                                        static_cast<int>(query.size()), 
                                        static_cast<int>(reference.size() - del_len)});
    
    while (shifts < total_possible_shifts && 
           (match_len - shifts - 1) >= 0 &&
           (match_len + del_len - shifts - 1) < reference.size()) {
        char query_base = query[match_len - shifts - 1];
        char ref_base = reference[match_len + del_len - shifts - 1];

        if (query_base != ref_base) {
            break;
        }
        shifts++;
    }

    if (shifts == 0) {
        return false;
    }

    match_len -= shifts;
    return true;
}

bool right_align_trailing_deletion(
    const std::string& reference,
    const std::string& query,
    int& del_pos,
    int del_len,
    int& match_len,
    int max_shifts)
{
    if (match_len == 0) {
        return true;
    }

    int shifts = 0;
    int total_possible_shifts = std::min({match_len, max_shifts,
                                        static_cast<int>(query.size() - del_pos),
                                        static_cast<int>(reference.size() - del_pos - del_len)});

    while (shifts < total_possible_shifts &&
           (del_pos + shifts) < query.size() &&
           (del_pos + del_len + shifts) < reference.size()) {
        char query_base = query[del_pos + shifts];
        char ref_base = reference[del_pos + del_len + shifts];

        if (query_base != ref_base) {
            break;
        }
        shifts++;
    }

    if (shifts == 0) {
        return false;
    }

    del_pos += shifts;
    match_len -= shifts;
    return true;
}


void recompute_identity_metrics(const std::string& cigar,
                              int& matches,
                              int& mismatches,
                              int& insertions,
                              int& insertion_events,
                              int& deletions,
                              int& deletion_events,
                              double& gap_compressed_identity,
                              double& blast_identity) {
    matches = mismatches = insertions = insertion_events = deletions = deletion_events = 0;

    size_t i = 0;
    char last_op = '\0';
    while (i < cigar.size()) {
        size_t j = i;
        while (j < cigar.size() && isdigit(cigar[j])) j++;
        int count = std::stoi(cigar.substr(i, j - i));
        char op = cigar[j];
        i = j + 1;

        switch (op) {
            case '=':
                matches += count;
                break;
            case 'X':
                mismatches += count;
                break;
            case 'I':
                insertions += count;
                if (last_op != 'I') {
                    insertion_events++;
                }
                break;
            case 'D':
                deletions += count;
                if (last_op != 'D') {
                    deletion_events++;
                }
                break;
        }
        last_op = op;
    }

    int total_columns = matches + mismatches + insertions + deletions;
    blast_identity = total_columns > 0 ? (double)matches / total_columns : 0.0;

    int total_differences = mismatches + insertion_events + deletion_events;
    int total_bases = matches + total_differences;
    gap_compressed_identity = total_bases > 0 ? (double)matches / total_bases : 0.0;
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

std::string merge_cigar_operations(const std::string& cigar) {
    if (cigar.empty()) return cigar;
    
    std::string result;
    size_t i = 0;
    
    // Parse first operation
    size_t j = i;
    while (j < cigar.size() && isdigit(cigar[j])) j++;
    int current_count = std::stoi(cigar.substr(i, j - i));
    char current_op = cigar[j];
    i = j + 1;
    
    // Process remaining operations
    while (i < cigar.size()) {
        j = i;
        while (j < cigar.size() && isdigit(cigar[j])) j++;
        int next_count = std::stoi(cigar.substr(i, j - i));
        char next_op = cigar[j];
        i = j + 1;
        
        // If operations are the same, merge them
        if (next_op == current_op) {
            current_count += next_count;
        } else {
            // Write out current operation and start new one
            result += std::to_string(current_count) + current_op;
            current_count = next_count;
            current_op = next_op;
        }
    }
    
    // Write final operation
    result += std::to_string(current_count) + current_op;
    
    return result;
}

void verify_cigar_alignment(const std::string& cigar,
                           const char* query_seq,
                           const char* target_seq,
                           int64_t query_start,
                           int64_t target_start,
                           int64_t query_length,
                           int64_t target_length) {
    size_t q_pos = 0;  // position in query sequence
    size_t t_pos = 0;  // position in target sequence

    size_t i = 0;
    while (i < cigar.size()) {
        size_t j = i;
        // Get the length of the operation
        while (j < cigar.size() && isdigit(cigar[j])) j++;
        int op_len = std::stoi(cigar.substr(i, j - i));
        char op = cigar[j];
        i = j + 1;

        switch (op) {
            case '=':  // match
                // Verify that query_seq[q_pos .. q_pos + op_len] matches target_seq[t_pos .. t_pos + op_len]
                for (int k = 0; k < op_len; ++k) {
                    char q_char = query_seq[q_pos + k];
                    char t_char = target_seq[t_pos + k];
                    // Check bounds before comparing
                    if (q_pos + k >= query_length || t_pos + k >= target_length) {
                        std::cerr << "[wfmash::align] Error: Position out of bounds during alignment verification "
                                  << "at query pos " << query_start + q_pos + k
                                  << " vs target pos " << target_start + t_pos + k << "\n";
                        exit(1);
                    }
                    if (q_char != t_char) {
                        // Print error message
                        std::cerr << "[wfmash::align] Error: Mismatch at position "
                                  << "query pos " << query_start + q_pos + k
                                  << " vs target pos " << target_start + t_pos + k
                                  << ": query char '" << q_char
                                  << "' vs target char '" << t_char << "' in '=' operation of CIGAR.\n";
                        
                        // Calculate context range for query sequence
                        size_t query_context_start = (q_pos + k >= 5) ? q_pos + k - 5 : 0;
                        size_t query_context_end = std::min(q_pos + k + 5, static_cast<size_t>(query_length - 1));
                
                        // Calculate context range for target sequence
                        size_t target_context_start = (t_pos + k >= 5) ? t_pos + k - 5 : 0;
                        size_t target_context_end = std::min(t_pos + k + 5, static_cast<size_t>(target_length - 1));
                        
                        // Extract context sequences
                        std::string query_context(query_seq + query_context_start, query_seq + query_context_end + 1);
                        std::string target_context(target_seq + target_context_start, target_seq + target_context_end + 1);
                        
                        // Print the alignment state
                        std::cerr << "[DEBUG] Alignment state at mismatch:\n";
                        std::cerr << "CIGAR string: " << cigar << "\n";
                        std::cerr << "Query sequence around mismatch (positions " << query_start + query_context_start
                                 << "-" << query_start + query_context_end << "):\n"
                                 << query_context << "\n";
                        std::cerr << "Target sequence around mismatch (positions " << target_start + target_context_start
                                 << "-" << target_start + target_context_end << "):\n"
                                 << target_context << "\n";
                        std::cerr << "q_pos: " << q_pos << ", t_pos: " << t_pos << ", k: " << k << "\n";
                        std::cerr << "Query sequence length: " << query_length << ", Target sequence length: " << target_length << "\n";
                        std::cerr << "Total query start: " << query_start << ", Total target start: " << target_start << "\n";
                        
                        exit(1);
                    }
                }
                q_pos += op_len;
                t_pos += op_len;
                break;
            case 'X':  // mismatch
                q_pos += op_len;
                t_pos += op_len;
                break;
            case 'I':  // insertion in target; nucleotides present in query
                q_pos += op_len;
                break;
            case 'D':  // deletion in target; nucleotides absent in query
                t_pos += op_len;
                break;
            default:
                std::cerr << "[wfmash::align] Error: Unsupported CIGAR operation '" << op << "'.\n";
                exit(1);
        }
    }
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

std::string adjust_cigar_string(const std::string& cigar,
                               const std::string& query_seq,
                               const std::string& target_seq,
                               int64_t query_start,
                               int64_t target_start,
                               uint64_t max_shift) {
    std::cerr << "[DEBUG] Original CIGAR string: " << cigar << "\n";
    
    // Find the first two operations
    size_t first_op_end = 0;
    while (first_op_end < cigar.size() && isdigit(cigar[first_op_end])) first_op_end++;
    if (first_op_end >= cigar.size()) return cigar;
    
    size_t second_op_start = first_op_end + 1;
    size_t second_op_end = second_op_start;
    while (second_op_end < cigar.size() && isdigit(cigar[second_op_end])) second_op_end++;
    if (second_op_end >= cigar.size()) return cigar;
    
    // Extract the operations
    int first_count = std::stoi(cigar.substr(0, first_op_end));
    char first_op = cigar[first_op_end];
    int second_count = std::stoi(cigar.substr(second_op_start, second_op_end - second_op_start));
    char second_op = cigar[second_op_end];
    
    // Check if we need to swap leading operations
    if (first_op == '=' && second_op == 'D') {
        // Check if swapping is valid by verifying sequence matches after the deletion
        bool can_swap = true;
        for (int k = 0; k < first_count && 
             (query_start + k) < query_seq.size() && 
             (target_start + second_count + k) < target_seq.size(); ++k) {
            if (query_seq[query_start + k] != 
                target_seq[target_start + second_count + k]) {
                can_swap = false;
                break;
            }
        }
        
        if (can_swap) {
            // Directly construct the swapped string
            return std::to_string(second_count) + 'D' + 
                   std::to_string(first_count) + '=' +
                   cigar.substr(second_op_end + 1);
        }
    }
    
    // Find and check the last two operations
    size_t last_op_end = cigar.size() - 1;
    size_t last_op_start = last_op_end;
    while (last_op_start > 0 && isdigit(cigar[last_op_start - 1])) last_op_start--;
    
    size_t second_last_op_end = last_op_start - 1;
    if (second_last_op_end <= 0) return cigar;
    size_t second_last_op_start = second_last_op_end;
    while (second_last_op_start > 0 && isdigit(cigar[second_last_op_start - 1])) second_last_op_start--;
    
    if (second_last_op_start == 0) return cigar;
    
    // Extract the last operations
    int last_count = std::stoi(cigar.substr(last_op_start, last_op_end - last_op_start));
    char last_op = cigar[last_op_end];
    int second_last_count = std::stoi(cigar.substr(second_last_op_start, second_last_op_end - second_last_op_start));
    char second_last_op = cigar[second_last_op_end];
    
    // Check if we need to swap trailing operations
    if (second_last_op == 'D' && last_op == '=') {
        int query_pos = query_start + query_seq.size() - last_count - second_last_count;
        int target_pos = target_start + target_seq.size() - last_count - second_last_count;
        
        bool can_swap = true;
        for (int k = 0; k < last_count && can_swap; ++k) {
            int64_t q_idx = query_pos + k;
            int64_t t_idx = target_pos + k;
            
            if (q_idx >= query_seq.size() || t_idx >= target_seq.size() ||
                query_seq[q_idx] != target_seq[t_idx]) {
                can_swap = false;
                break;
            }
        }
        
        if (can_swap) {
            // Directly construct the swapped string
            return cigar.substr(0, second_last_op_start) +
                   std::to_string(last_count) + '=' +
                   std::to_string(second_last_count) + 'D';
        }
    }
    
    return cigar;
}


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

          // Save words into currentRecord
          {
              currentRecord.qId = tokens[0];
              currentRecord.qStartPos = std::stoi(tokens[2]);
              currentRecord.qEndPos = std::stoi(tokens[3]);
              currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
              currentRecord.refId = tokens[5];
              const uint64_t ref_len = std::stoi(tokens[6]);
              // Apply target padding while ensuring we don't go below 0 or above reference length
              uint64_t rStartPos = std::stoi(tokens[7]);
              uint64_t rEndPos = std::stoi(tokens[8]);
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
              currentRecord.rStartPos = rStartPos;
              currentRecord.rEndPos = rEndPos;
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

    // Calculate padded positions ensuring they stay within bounds
    int64_t ref_fetch_start = static_cast<int64_t>(currentRecord.rStartPos) - 
                             static_cast<int64_t>(param.wflign_max_len_minor);
    uint64_t head_padding = param.wflign_max_len_minor;
    if (ref_fetch_start < 0) {
        head_padding += ref_fetch_start;  // Reduce padding
        ref_fetch_start = 0;
    }

    int64_t ref_fetch_end = static_cast<int64_t>(currentRecord.rEndPos - 1) + 
                           static_cast<int64_t>(param.wflign_max_len_minor);
    uint64_t tail_padding = param.wflign_max_len_minor;
    if (ref_fetch_end >= ref_size) {
        tail_padding -= (ref_fetch_end - (ref_size - 1));
        ref_fetch_end = ref_size - 1;
    }

    // Extract reference sequence with validated bounds
    int64_t ref_len;
    char* ref_seq = faidx_fetch_seq64(ref_faidx, currentRecord.refId.c_str(),
                                      ref_fetch_start,
                                      ref_fetch_end, &ref_len);

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
    queryRegionStrand[query_seq.size()] = '\0';

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
        rec->currentRecord.mashmap_estimated_identity);

    // Get the alignment output as a string
    std::string alignment_output = output.str();

    // Extract and adjust the CIGAR string
    size_t cg_pos = alignment_output.find("cg:Z:");
    if (cg_pos != std::string::npos) {
        size_t cigar_start = cg_pos + 5;
        size_t cigar_end = alignment_output.find('\t', cigar_start);
        if (cigar_end == std::string::npos) cigar_end = alignment_output.length();
        std::string original_cigar = alignment_output.substr(cigar_start, cigar_end - cigar_start);

        // Add debugging print statements
        std::cerr << "[DEBUG] Alignment output before modification:\n" << alignment_output << "\n";
        std::cerr << "[DEBUG] Original CIGAR string: " << original_cigar << "\n";
        std::cerr << "[DEBUG] Target positions: start=" << rec->currentRecord.rStartPos 
                  << ", end=" << rec->currentRecord.rEndPos << "\n";
        std::cerr << "[DEBUG] Query positions: start=" << rec->currentRecord.qStartPos 
                  << ", end=" << rec->currentRecord.qEndPos << "\n";

        // Adjust the CIGAR string
        std::string adjusted_cigar = adjust_cigar_string(original_cigar,
                                                       queryRegionStrand.data(),
                                                       ref_seq_ptr,
                                                       rec->currentRecord.qStartPos,
                                                       rec->currentRecord.rStartPos,
                                                       param.target_padding);

        // Merge any equivalent successive operations
        adjusted_cigar = merge_cigar_operations(adjusted_cigar);

        // Skip empty alignments
        if (adjusted_cigar.empty()) {
            return "";
        }

        // Use original sequence pointers
        char* adjusted_ref_seq_ptr = ref_seq_ptr;
        char* adjusted_query_seq_ptr = queryRegionStrand.data();

        // Debug output for adjusted sequences
        std::cerr << "[DEBUG] Adjusted positions and sequences:\n"
                  << "Adjusted target_start: " << rec->currentRecord.rStartPos << "\n"
                  << "Adjusted target_end: " << rec->currentRecord.rEndPos << "\n"
                  << "Adjusted query_start: " << rec->currentRecord.qStartPos << "\n"
                  << "Adjusted query_end: " << rec->currentRecord.qEndPos << "\n";

        // Print reference sequence context
        int ref_context_size = 10;
        int64_t ref_seq_offset = rec->currentRecord.rStartPos - rec->refStartPos;
        std::string ref_context = ref_seq.substr(
            std::max<int64_t>(0, ref_seq_offset - ref_context_size),
            std::min<size_t>(ref_seq.size() - ref_seq_offset + ref_context_size, 2 * ref_context_size)
        );
        std::cerr << "Reference sequence around adjusted target_start (positions "
                  << rec->currentRecord.rStartPos - ref_context_size << " to " << rec->currentRecord.rStartPos + ref_context_size << "):\n"
                  << ref_context << "\n";

        // Print query sequence context
        int query_context_size = 10;
        int64_t query_seq_offset = rec->currentRecord.qStartPos - rec->currentRecord.qStartPos;
        std::string query_context = std::string(queryRegionStrand.data()).substr(
            std::max<int64_t>(0, query_seq_offset - query_context_size),
            std::min<size_t>(queryRegionStrand.size() - query_seq_offset + query_context_size, 2 * query_context_size)
        );
        std::cerr << "Query sequence around adjusted query_start (positions "
                  << rec->currentRecord.qStartPos - query_context_size << " to " << rec->currentRecord.qStartPos + query_context_size << "):\n"
                  << query_context << "\n";

        // Verify the alignment matches in '=' operations using adjusted pointers
        verify_cigar_alignment(adjusted_cigar,
                             adjusted_query_seq_ptr,
                             adjusted_ref_seq_ptr,
                             rec->currentRecord.qStartPos,
                             rec->currentRecord.rStartPos,
                             rec->queryLen,
                             rec->refLen);

        // Recompute identity metrics
        int matches, mismatches, insertions, insertion_events, deletions, deletion_events;
        double gap_compressed_identity, blast_identity;
        recompute_identity_metrics(adjusted_cigar, matches, mismatches, insertions, insertion_events,
                                   deletions, deletion_events, gap_compressed_identity, blast_identity);

        // Update the alignment output with new positions and metrics
        std::string updated_output;
        std::istringstream iss(alignment_output);
        std::string field;
        std::vector<std::string> fields;
        
        while (std::getline(iss, field, '\t')) {
            fields.push_back(field);
        }

        // Update positions based on format
        if (!param.sam_format) {  // PAF format
            // Query positions (0-based)
            fields[2] = std::to_string(rec->currentRecord.qStartPos);
            fields[3] = std::to_string(rec->currentRecord.qEndPos);
            // Target positions (0-based)
            fields[7] = std::to_string(rec->currentRecord.rStartPos);
            fields[8] = std::to_string(rec->currentRecord.rEndPos);
        } else {  // SAM format
            // Target position (1-based)
            fields[3] = std::to_string(rec->currentRecord.rStartPos + 1);
            // If necessary, adjust the query positions stored in optional fields
        }

        // Replace the original CIGAR with the adjusted one
        if (!param.sam_format) {
            // Replace the 'cg:Z:' tag
            for (size_t i = 12; i < fields.size(); ++i) {
                if (fields[i].substr(0, 5) == "cg:Z:") {
                    fields[i] = "cg:Z:" + adjusted_cigar;
                    break;
                }
            }
        } else {
            // Replace the CIGAR field
            fields[5] = adjusted_cigar;
        }

        // Update or replace 'gi:f:' and 'bi:f:' tags
        bool gi_found = false, bi_found = false;
        std::ostringstream gi_stream, bi_stream;
        gi_stream << std::fixed << std::setprecision(6) << "gi:f:" << gap_compressed_identity;
        bi_stream << std::fixed << std::setprecision(6) << "bi:f:" << blast_identity;

        for (size_t i = 12; i < fields.size(); ++i) {
            if (fields[i].substr(0, 5) == "gi:f:") {
                fields[i] = gi_stream.str();
                gi_found = true;
            }
            if (fields[i].substr(0, 5) == "bi:f:") {
                fields[i] = bi_stream.str();
                bi_found = true;
            }
        }
        if (!gi_found) {
            fields.push_back(gi_stream.str());
        }
        if (!bi_found) {
            fields.push_back(bi_stream.str());
        }

        // Reconstruct the updated alignment output
        updated_output = fields[0];
        for (size_t i = 1; i < fields.size(); ++i) {
            updated_output += "\t" + fields[i];
        }
        // Add a newline only if one is not already present
        if (!updated_output.empty() && updated_output.back() != '\n') {
            updated_output += "\n";
        }

        // Verify the alignment matches in '=' operations
        verify_cigar_alignment(adjusted_cigar,
                             queryRegionStrand.data(),
                             ref_seq_ptr,
                             rec->queryStartPos,
                             rec->refStartPos,
                             rec->queryLen,
                             rec->refLen);

        std::cerr << "[DEBUG] Alignment output after modification:\n" << updated_output << "\n";
        return updated_output;
    }

    // If 'cg:Z:' tag is not found, return the original alignment output
    return alignment_output;
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
    faidx_t* local_ref_faidx = fai_load(param.refSequences.front().c_str());
    faidx_t* local_query_faidx = fai_load(param.querySequences.front().c_str());

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
            
            // Push the alignment output to the paf_queue
            paf_queue.push(new std::string(std::move(alignment_output)));
            
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
