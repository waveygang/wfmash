#include <cassert>
#include <chrono>
#include <string>

#include "wflign.hpp"
#include "wflign_patch.hpp"
#include "wflign_swizzle.hpp"


// Namespaces
namespace wflign {
namespace wavefront {

/*
* Configuration
*/
#define MIN_WF_LENGTH            256

std::string erode_short_matches_in_cigar(const std::string& cigar, int max_match_length = 3, bool is_head_cigar = true) {
    if (cigar.length() < 6) return cigar; // Too short for indel-match-indel pattern
    
    // Parse CIGAR into operations
    std::vector<std::pair<int, char>> ops;
    ops.reserve(cigar.length() / 2); // Pre-allocate to avoid reallocations
    
    size_t pos = 0;
    while (pos < cigar.length()) {
        size_t start = pos;
        while (pos < cigar.length() && std::isdigit(cigar[pos])) pos++;
        if (pos < cigar.length()) {
            int count = std::stoi(cigar.substr(start, pos - start));
            char type = cigar[pos++];
            ops.push_back({count, type});
        }
    }
    
    // Check if we have enough operations
    if (ops.size() < 3) return cigar;
    
    // Process operations, marking some for removal
    bool modified = false;
    
    // Define the range of operations to check based on whether this is head or tail CIGAR
    size_t start_idx = 1;
    size_t end_idx = ops.size() - 1;
    
    if (is_head_cigar) {
        // Only check the first 3 operations for head CIGAR
        end_idx = std::min(end_idx, (size_t)3);
    } else {
        // Only check the last 3 operations for tail CIGAR
        start_idx = std::max(start_idx, ops.size() - 3);
    }
    
    for (size_t i = start_idx; i < end_idx; i++) {
        // Check for indel-match-indel pattern with short match
        bool is_match = (ops[i].second == 'M' || ops[i].second == '=' || ops[i].second == 'X');
        bool prev_is_indel = (ops[i-1].second == 'I' || ops[i-1].second == 'D');
        bool next_is_indel = (ops[i+1].second == 'I' || ops[i+1].second == 'D');
        
        if (is_match && ops[i].first <= max_match_length && prev_is_indel && next_is_indel) {
            // Only erode if:
            // 1. The indels are of different types (I and D or D and I)
            // 2. Both indels are longer than the match
            if (((ops[i-1].second == 'I' && ops[i+1].second == 'D') || 
                 (ops[i-1].second == 'D' && ops[i+1].second == 'I')) &&
                (ops[i-1].first > ops[i].first && ops[i+1].first > ops[i].first)) {
                // Increase both surrounding indels
                ops[i-1].first += ops[i].first;
                ops[i+1].first += ops[i].first;
                
                // Mark for removal by setting count to 0
                ops[i].first = 0;
                modified = true;
            }
        }
    }
    
    // Fast path: if nothing changed, return original
    if (!modified) return cigar;
    
    // Build new operations vector, excluding those with count 0 and merging consecutive same-type ops
    std::vector<std::pair<int, char>> merged_ops;
    merged_ops.reserve(ops.size());
    
    for (const auto& op : ops) {
        if (op.first > 0) {
            // If we have a previous operation of the same type, merge them
            if (!merged_ops.empty() && merged_ops.back().second == op.second) {
                merged_ops.back().first += op.first;
            } else {
                merged_ops.push_back(op);
            }
        }
    }
    
    // Build result string from merged operations
    std::string result;
    result.reserve(cigar.length());
    
    for (const auto& op : merged_ops) {
        result += std::to_string(op.first) + op.second;
    }
    
    return result;
}

void do_biwfa_alignment(
    const std::string& query_name,
    char* const query,
    const uint64_t query_total_length,
    const uint64_t query_offset,
    const uint64_t query_length,
    const bool query_is_rev,
    const std::string& target_name,
    char* const target,
    const uint64_t target_total_length,
    const uint64_t target_offset,
    const uint64_t target_length,
    std::ostream& out,
    const wflign_penalties_t& penalties,
    const bool emit_md_tag,
    const bool paf_format_else_sam,
    const bool no_seq_in_sam,
    const bool disable_chain_patching,
    const float min_identity,
    const uint64_t wflign_max_len_minor,
    const float mashmap_estimated_identity,
    const int32_t chain_id,
    const int32_t chain_length,
    const int32_t chain_pos) {
    
    // Create WFA aligner for the main alignment
    wfa::WFAlignerGapAffine2Pieces wf_aligner(
        0,  // match
        penalties.mismatch,
        penalties.gap_opening1,
        penalties.gap_extension1,
        penalties.gap_opening2,
        penalties.gap_extension2,
        wfa::WFAligner::Alignment,
        wfa::WFAligner::MemoryUltralow);
    wf_aligner.setHeuristicNone();
    
    // Perform the main end-to-end alignment first
    const int status = wf_aligner.alignEnd2End(target, (int)target_length, query, (int)query_length);
    
    if (status != 0) {
        return; // Alignment failed
    }

    // Create alignment record on stack
    alignment_t aln;
    aln.ok = true;
    aln.j = 0;
    aln.i = 0;
    aln.query_length = query_length;
    aln.target_length = target_length;
    aln.is_rev = false;
    
    // Copy alignment CIGAR from the main alignment
    wflign_edit_cigar_copy(wf_aligner, &aln.edit_cigar);
    std::string main_cigar = wfa_edit_cigar_to_string(aln.edit_cigar);
    
    if (!disable_chain_patching) {
        // Set up constants for patching
        const int MIN_PATCH_LENGTH = 128;       // Minimum length to expose for patching
        const int MAX_ERODE_LENGTH = 4096;      // Maximum erosion before stopping
        const int MIN_CONSECUTIVE_MATCHES = 11; // Minimum consecutive matches to stop erosion

        // Helper function to parse a single CIGAR operation (e.g., "10M")
        auto parse_cigar_op = [](const std::string& cigar, size_t& pos) -> std::pair<int, char> {
            size_t start = pos;
            while (pos < cigar.length() && isdigit(cigar[pos])) pos++;
            int count = std::stoi(cigar.substr(start, pos - start));
            char op = cigar[pos++];
            return {count, op};
        };
        
        // Helper function to convert CIGAR string from long form to short form (run-length encoded)
        auto compress_cigar = [](const std::string& long_cigar) -> std::string {
            if (long_cigar.empty()) return "";
            
            std::string short_cigar;
            char prev_op = long_cigar[0];
            int count = 1;
            
            for (size_t i = 1; i < long_cigar.length(); i++) {
                char op = long_cigar[i];
                if (op == prev_op) {
                    count++;
                } else {
                    // Convert M to = for consistency with wfa_edit_cigar_to_string
                    char out_op = (prev_op == 'M') ? '=' : prev_op;
                    short_cigar += std::to_string(count) + out_op;
                    prev_op = op;
                    count = 1;
                }
            }
            
            // Add the last operation
            char out_op = (prev_op == 'M') ? '=' : prev_op;
            short_cigar += std::to_string(count) + out_op;
            
            return short_cigar;
        };

        // Helper function to merge adjacent CIGAR operations at concatenation point
        auto merge_adjacent_ops = [](const std::string& cigar1, const std::string& cigar2) -> std::string {
            if (cigar1.empty()) return cigar2;
            if (cigar2.empty()) return cigar1;
            
            // Find the last operation in cigar1
            size_t pos1 = cigar1.length() - 1;
            while (pos1 > 0 && !std::isdigit(cigar1[pos1])) pos1--;
            char op1 = cigar1.back();
            size_t start1 = cigar1.rfind(std::isdigit(cigar1[pos1]) ? cigar1[pos1] : ' ');
            while (start1 > 0 && std::isdigit(cigar1[start1-1])) start1--;
            int count1 = std::stoi(cigar1.substr(start1, pos1 - start1 + 1));
            
            // Find the first operation in cigar2
            size_t pos2 = 0;
            while (pos2 < cigar2.length() && std::isdigit(cigar2[pos2])) pos2++;
            if (pos2 >= cigar2.length()) return cigar1 + cigar2;
            char op2 = cigar2[pos2];
            int count2 = std::stoi(cigar2.substr(0, pos2));
            
            // If operations match, merge them
            if (op1 == op2) {
                return cigar1.substr(0, start1) + 
                    std::to_string(count1 + count2) + op1 + 
                    cigar2.substr(pos2 + 1);
            }
            
            return cigar1 + cigar2;
        };

        // Perform head patching
        {
            u_int64_t query_eroded = 0;
            u_int64_t target_eroded = 0;
            size_t cigar_pos = 0;
            size_t erode_end_pos = 0;
            bool found_consecutive_matches = false;

            // Continue eroding until stopping conditions are met
            while (cigar_pos < main_cigar.length()) {
                auto [count, op] = parse_cigar_op(main_cigar, cigar_pos);

                if (op == '=' && count >= MIN_CONSECUTIVE_MATCHES) {
                    found_consecutive_matches = true;
                }

                // If we've found MIN_CONSECUTIVE_MATCHES and satisfied MIN_PATCH_LENGTH, stop
                if (found_consecutive_matches && 
                    query_eroded >= MIN_PATCH_LENGTH && target_eroded >= MIN_PATCH_LENGTH) {
                    break;
                }
                // Stop if we've reached MAX_ERODE_LENGTH
                if (query_eroded >= MAX_ERODE_LENGTH || target_eroded >= MAX_ERODE_LENGTH) {
                    break;
                }

                // Update counts based on operation
                if (op == 'M' || op == 'X' || op == '=') {
                    query_eroded += count;
                    target_eroded += count;
                } else if (op == 'I') {
                    query_eroded += count;
                } else if (op == 'D') {
                    target_eroded += count;
                }
                erode_end_pos = cigar_pos;
            }

            if (query_eroded > 3 || target_eroded > 3) {
                // Create a dedicated aligner for head patching
                wfa::WFAlignerGapAffine2Pieces head_aligner(
                    0,  // match
                    penalties.mismatch,
                    penalties.gap_opening1,
                    penalties.gap_extension1,
                    penalties.gap_opening2,
                    penalties.gap_extension2,
                    wfa::WFAligner::Alignment,
                    wfa::WFAligner::MemoryMed);
                head_aligner.setHeuristicNone();
                
                // Extract sequences for head patching
                int head_query_length = query_eroded;
                int head_target_length = target_eroded;
                
                std::string head_query_str(query, head_query_length);
                std::string head_target_str(target, head_target_length);
                
                // Do semi-global alignment for head patching
                // Allow free gaps at the beginning of both sequences
                const int head_status = head_aligner.alignEndsFree(
                    head_target_str,
                    target_eroded, 0,
                    head_query_str,
                    query_eroded, 0
                );

                if (head_status == 0) {
                    // Get the head CIGAR in long form
                    std::string head_cigar_long = head_aligner.getAlignment();

                    // Convert to short form using our helper function
                    std::string head_cigar_short = compress_cigar(head_cigar_long);
                    
                    head_cigar_short = erode_short_matches_in_cigar(head_cigar_short, 3);
                    
                    // Remove the eroded part from the beginning of main_cigar
                    main_cigar = merge_adjacent_ops(head_cigar_short, main_cigar.substr(erode_end_pos));
                }
            }
        }
        
        // Perform tail patching
        {
            // For the tail, we need to parse the entire CIGAR first to know where to start
            std::vector<std::pair<int, char>> cigar_ops;
            size_t pos = 0;
            while (pos < main_cigar.length()) {
                cigar_ops.push_back(parse_cigar_op(main_cigar, pos));
            }
            
            u_int64_t query_eroded = 0;
            u_int64_t target_eroded = 0;
            size_t erode_start_idx = cigar_ops.size();
            bool found_consecutive_matches = false;

            // Work backwards from the end of the CIGAR
            for (int i = cigar_ops.size() - 1; i >= 0; i--) {
                auto [count, op] = cigar_ops[i];
                                
                if (op == '=' && count >= MIN_CONSECUTIVE_MATCHES) {
                    found_consecutive_matches = true;
                }

                // If we've found MIN_CONSECUTIVE_MATCHES and satisfied MIN_PATCH_LENGTH, stop
                if (found_consecutive_matches && 
                    query_eroded >= MIN_PATCH_LENGTH && target_eroded >= MIN_PATCH_LENGTH) {
                    break;
                }
                // Stop if we've reached MAX_ERODE_LENGTH
                if (query_eroded >= MAX_ERODE_LENGTH || target_eroded >= MAX_ERODE_LENGTH) {
                    break;
                }

                // Update counts based on operation
                if (op == 'M' || op == 'X' || op == '=') {
                    query_eroded += count;
                    target_eroded += count;
                } else if (op == 'I') {
                    query_eroded += count;
                } else if (op == 'D') {
                    target_eroded += count;
                }
                erode_start_idx = i;
            }

            if (query_eroded > 3 || target_eroded > 3) {
                // Create a dedicated aligner for tail patching
                wfa::WFAlignerGapAffine2Pieces tail_aligner(
                    0,  // match
                    penalties.mismatch,
                    penalties.gap_opening1,
                    penalties.gap_extension1,
                    penalties.gap_opening2,
                    penalties.gap_extension2,
                    wfa::WFAligner::Alignment,
                    wfa::WFAligner::MemoryMed);
                tail_aligner.setHeuristicNone();
                
                // Extract sequences for tail patching
                int tail_query_length = query_eroded;
                int tail_target_length = target_eroded;
                
                // Get the starting positions for the tail patching
                char* query_tail = query + query_length - tail_query_length;
                char* target_tail = target + target_length - tail_target_length;
                
                std::string tail_query_str(query_tail, tail_query_length);
                std::string tail_target_str(target_tail, tail_target_length);
                
                // Do semi-global alignment for tail patching
                // Allow free gaps at the end of both sequences
                const int tail_status = tail_aligner.alignEndsFree(
                    tail_target_str,
                    0, tail_target_length,  // textBeginFree, textEndFree
                    tail_query_str,
                    0, tail_query_length   // patternBeginFree, patternEndFree
                );
                
                if (tail_status == 0) {
                    // Get the tail CIGAR in long form
                    std::string tail_cigar_long = tail_aligner.getAlignment();
                    
                    // Convert to short form using our helper function
                    std::string tail_cigar_short = compress_cigar(tail_cigar_long);
                    
                    tail_cigar_short = erode_short_matches_in_cigar(tail_cigar_short, 3, false);
                    
                    // Rebuild the CIGAR string up to the erode_start_idx
                    std::string truncated_cigar;
                    for (size_t i = 0; i < erode_start_idx; i++) {
                        truncated_cigar += std::to_string(cigar_ops[i].first) + cigar_ops[i].second;
                    }
                    
                    // Combine the truncated main CIGAR with the tail CIGAR
                    main_cigar = merge_adjacent_ops(truncated_cigar, tail_cigar_short);
                }
            }
        }

    }

    // Try swizzling the CIGAR at both ends
    std::string swizzled = try_swap_start_pattern(main_cigar, query, target, 0, 0);
    if (swizzled != main_cigar) {
        main_cigar = swizzled;
    }

    swizzled = try_swap_end_pattern(main_cigar, query, target, 0, 0);
    if (swizzled != main_cigar) {
        main_cigar = swizzled;
    }
    
    // Write alignment
    if (paf_format_else_sam) {
        write_alignment_paf(
            out,
            aln,
            main_cigar,
            query_name,
            query_total_length,
            query_offset,
            query_length,
            query_is_rev,
            target_name,
            target_total_length,
            target_offset,
            target_length,
            min_identity,
            mashmap_estimated_identity,
            chain_id,
            chain_length,
            chain_pos);
    } else {
        // Write SAM output
        write_alignment_sam(
            out,
            aln,
            main_cigar,
            query_name,
            query_total_length,
            query_offset,
            query_length,
            query_is_rev,
            target_name,
            target_total_length,
            target_offset,
            target_length,
            min_identity,
            mashmap_estimated_identity,
            no_seq_in_sam,
            emit_md_tag,
            query,
            target,
            0,
            chain_id,
            chain_length,
            chain_pos); // No target pointer shift for biwfa
    }
}

/*
* Configuration
*/
#define MAX_LEN_FOR_STANDARD_WFA 1000

/*
* Utils
*/
inline uint64_t encode_pair(
    int v,
    int h) {
    return ((uint64_t)v << 32) | (uint64_t)h;
}
inline void decode_pair(uint64_t pair, int *v, int *h) {
    *v = (int)(pair >> 32);
    *h = (int)(pair & 0x00000000FFFFFFFF);
}
void clean_up_sketches(std::vector<std::vector<rkmh::hash_t>*> &sketches) {
    // The C++ language guarantees that `delete p` will do nothing if p is equal to NULL
    for (auto &s : sketches) {
        delete s;
        s = nullptr;
    }
}

/*
* Setup
*/
WFlign::WFlign(
    const uint16_t segment_length,
    const float min_identity,
    const bool force_wflign_,
    const int wfa_mismatch_score,
    const int wfa_gap_opening_score,
    const int wfa_gap_extension_score,
    const int wfa_patching_mismatch_score,
    const int wfa_patching_gap_opening_score1,
    const int wfa_patching_gap_extension_score1,
    const int wfa_patching_gap_opening_score2,
    const int wfa_patching_gap_extension_score2,
    const float mashmap_estimated_identity,
    const int wflign_mismatch_score,
    const int wflign_gap_opening_score,
    const int wflign_gap_extension_score,
    const float wflign_max_mash_dist,
    const int wflign_min_wavefront_length,
    const int wflign_max_distance_threshold,
    const uint64_t wflign_max_len_major,
    const uint64_t wflign_max_len_minor,
    const int erode_k,
    const int64_t chain_gap,
    const int min_inversion_length,
    const int max_patching_score) {
    // Parameters
    this->segment_length = segment_length;
    this->min_identity = min_identity;

    this->force_wflign = force_wflign_;

    this->wfa_mismatch_score = wfa_mismatch_score;
    this->wfa_gap_opening_score = wfa_gap_opening_score;
    this->wfa_gap_extension_score = wfa_gap_extension_score;

    this->wfa_patching_mismatch_score = wfa_patching_mismatch_score;
    this->wfa_patching_gap_opening_score1 = wfa_patching_gap_opening_score1;
    this->wfa_patching_gap_extension_score1 = wfa_patching_gap_extension_score1;
    this->wfa_patching_gap_opening_score2 = wfa_patching_gap_opening_score2;
    this->wfa_patching_gap_extension_score2 = wfa_patching_gap_extension_score2;

    this->mashmap_estimated_identity = mashmap_estimated_identity;
    this->wflign_mismatch_score = wflign_mismatch_score;
    this->wflign_gap_opening_score = wflign_gap_opening_score;
    this->wflign_gap_extension_score = wflign_gap_extension_score;
    this->wflign_max_mash_dist = wflign_max_mash_dist;
    this->wflign_min_wavefront_length = wflign_min_wavefront_length;
    this->wflign_max_distance_threshold = wflign_max_distance_threshold;
    this->wflign_max_len_major = wflign_max_len_major;
    this->wflign_max_len_minor = wflign_max_len_minor;
    this->erode_k = erode_k;
    this->chain_gap = chain_gap;
    this->max_patching_score = max_patching_score;
    this->min_inversion_length = min_inversion_length;
    // Query
    this->query_name = nullptr;
    this->query = nullptr;
    this->query_total_length = 0;
    this->query_offset = 0;
    this->query_length = 0;
    this->query_is_rev = false;
    // Target
    this->target_name = nullptr;
    this->target = nullptr;
    this->target_total_length = 0;
    this->target_offset = 0;
    this->target_length = 0;
    // Output
    this->out = nullptr;
#ifdef WFA_PNG_TSV_TIMING
    this->emit_tsv = false;
    this->out_tsv = nullptr;
    this->prefix_wavefront_plot_in_png = nullptr;
    this->wfplot_max_size = 0;
    this->emit_patching_tsv = false;
    this->out_patching_tsv = nullptr;
#endif
    this->merge_alignments = false;
    this->emit_md_tag = false;
    this->paf_format_else_sam = false;
    this->no_seq_in_sam = false;
    this->disable_chain_patching = false;
}
/*
* Output configuration
*/
void WFlign::set_output(
    std::ostream* const out,
#ifdef WFA_PNG_TSV_TIMING
    const bool emit_tsv,
    std::ostream* const out_tsv,
    const std::string &wfplot_filepath,
    const uint64_t wfplot_max_size,
    const bool emit_patching_tsv,
    std::ostream* const out_patching_tsv,
#endif
    const bool merge_alignments,
    const bool emit_md_tag,
    const bool paf_format_else_sam,
    const bool no_seq_in_sam) {
    this->out = out;
#ifdef WFA_PNG_TSV_TIMING
    this->emit_tsv = emit_tsv;
    this->out_tsv = out_tsv;
    this->prefix_wavefront_plot_in_png = &wfplot_filepath;
    this->wfplot_max_size = wfplot_max_size;
    this->emit_patching_tsv = emit_patching_tsv;
    this->out_patching_tsv = out_patching_tsv;
#endif
    this->merge_alignments = merge_alignments;
    this->emit_md_tag = emit_md_tag;
    this->paf_format_else_sam = paf_format_else_sam;
    this->no_seq_in_sam = no_seq_in_sam;
}
/*
* WFlambda
*/
int wflambda_extend_match(
    const int v,
    const int h,
    void* arguments) {
    // Extract arguments
    wflign_extend_data_t* extend_data = (wflign_extend_data_t*)arguments;
    // Expand arguments
    const WFlign& wflign = *(extend_data->wflign);
    const uint64_t query_length = wflign.query_length;
    const uint64_t target_length = wflign.target_length;
    const int step_size = extend_data->step_size;
    const int segment_length_to_use = extend_data->segment_length_to_use;
    const int pattern_length = extend_data->pattern_length;
    const int text_length = extend_data->text_length;
    robin_hood::unordered_flat_map<uint64_t,alignment_t*>& alignments = *(extend_data->alignments);
    std::vector<std::vector<rkmh::hash_t>*>& query_sketches = *(extend_data->query_sketches);
    std::vector<std::vector<rkmh::hash_t>*>& target_sketches = *(extend_data->target_sketches);
#ifdef WFA_PNG_TSV_TIMING
    // wfplots
    const bool emit_png = extend_data->emit_png;
    robin_hood::unordered_set<uint64_t>& high_order_dp_matrix_mismatch = *(extend_data->high_order_dp_matrix_mismatch);
#endif
    // Check match
    bool is_a_match = false;
    if (v >= 0 && h >= 0 && v < pattern_length && h < text_length) {
        const uint64_t k = encode_pair(v, h);
        const auto f = alignments.find(k); // high-level of WF-inception
        if (f != alignments.end()) {
            is_a_match = (alignments[k] != nullptr);
        } else {
            const int64_t query_begin = v * step_size;
            const int64_t target_begin = h * step_size;

            // The last fragment can be longer than segment_length_to_use (max 2*segment_length_to_use - 1)
            const uint16_t segment_length_to_use_q =
                    (v == pattern_length - 1) ? query_length - query_begin : segment_length_to_use;
            const uint16_t segment_length_to_use_t =
                    (h == text_length - 1) ? target_length - target_begin : segment_length_to_use;

            auto *aln = new alignment_t();
            const bool alignment_performed =
                    do_wfa_segment_alignment(
                            *wflign.query_name,
                            wflign.query,
                            query_sketches[v],
                            wflign.query_length,
                            query_begin,
                            *wflign.target_name,
                            wflign.target,
                            target_sketches[h],
                            wflign.target_length,
                            target_begin,
                            segment_length_to_use_q,
                            segment_length_to_use_t,
                            step_size,
                            extend_data,
                            *aln);
#ifdef WFA_PNG_TSV_TIMING
            if (wflign.emit_tsv) {
                // 0) Mis-match, alignment skipped
                // 1) Mis-match, alignment performed
                // 2) Match, alignment performed
                *(wflign.out_tsv) << v << "\t" << h << "\t"
                                  << (alignment_performed ? (aln->ok ? 2 : 1) : 0)
                                  << std::endl;
            }
#endif
#ifdef WFA_PNG_TSV_TIMING
            ++(extend_data->num_alignments);
#endif
            if (alignment_performed) {
#ifdef WFA_PNG_TSV_TIMING
                ++(extend_data->num_alignments_performed);
#endif
                if (aln->ok){
                    is_a_match = true;
                    alignments[k] = aln;
                } else {
                    alignments[k] = nullptr;
                }
            }
#ifdef WFA_PNG_TSV_TIMING
            else {
                if (emit_png) {
                    // Save only the mismatches, as they are not cached
                    high_order_dp_matrix_mismatch.insert(encode_pair(v, h));
                }
            }
#endif
            if (!is_a_match) {
                delete aln;
            }

            if (extend_data->num_sketches_allocated > extend_data->max_num_sketches_in_memory) {
                clean_up_sketches(query_sketches);
                clean_up_sketches(target_sketches);
                extend_data->num_sketches_allocated = 0;
            }

        }
    } else if (h < 0 || v < 0) { // It can be removed using an edit-distance
        // mode as high-level of WF-inception
        is_a_match = false;
    }
    return is_a_match;
}

int wflambda_trace_match(
    robin_hood::unordered_flat_map<uint64_t,alignment_t*>& alignments,
    wfa::WFAlignerGapAffine& wflambda_aligner,
    std::vector<alignment_t*>& trace,
    const int pattern_length,
    const int text_length) {
    // Retrieve CIGAR
    char* cigar_ops;
    int cigar_length;
    wflambda_aligner.getAlignment(&cigar_ops,&cigar_length);
    // Traverse CIGAR tracing matches
    int v = pattern_length-1;
    int h = text_length-1;
    int i, num_alignments = 0;
    for (i=cigar_length-1;i>=0;--i) {
        switch (cigar_ops[i]) {
            case 'I':      --h; break;
            case 'D': --v;      break;
            case 'X': --v; --h; break;
            case 'M': {
                // Add alignment to trace
                uint64_t k = encode_pair(v,h);
                alignment_t* aln = alignments[k];
                trace.push_back(aln);
                aln->keep = true;
                ++num_alignments;
                // Next
                --v;
                --h;
                break;
            }
        }
    }
    return num_alignments;
}
/*
* WFling align
*/
void WFlign::wflign_affine_wavefront(
    const std::string& query_name,
    char* const query,
    const uint64_t query_total_length,
    const uint64_t query_offset,
    const uint64_t query_length,
    const bool query_is_rev,
    const std::string& target_name,
    char* const target,
    const uint64_t target_total_length,
    const uint64_t target_offset,
    const uint64_t target_length) {
    // Set query
    this->query_name = &query_name;
    this->query = query;
    this->query_total_length = query_total_length;
    this->query_offset = query_offset;
    this->query_length = query_length;
    this->query_is_rev = query_is_rev;
    // Set target
    this->target_name = &target_name;
    this->target = target;
    this->target_total_length = target_total_length;
    this->target_offset = target_offset;
    this->target_length = target_length;

    if (query_offset + query_length > query_total_length ||
        target_offset + target_length > target_total_length) {
        return;
    }

    // Set penalties for wfa convex
    wflign_penalties_t wfa_convex_penalties;
    if (wfa_patching_mismatch_score > 0 && wfa_patching_gap_opening_score1 > 0 && wfa_patching_gap_extension_score1 > 0 && wfa_patching_gap_opening_score2 > 0 && wfa_patching_gap_extension_score2 > 0){
        wfa_convex_penalties.match = 0;
        wfa_convex_penalties.mismatch = wfa_patching_mismatch_score;
        wfa_convex_penalties.gap_opening1 = wfa_patching_gap_opening_score1;
        wfa_convex_penalties.gap_extension1 = wfa_patching_gap_extension_score1;
        wfa_convex_penalties.gap_opening2 = wfa_patching_gap_opening_score2;
        wfa_convex_penalties.gap_extension2 = wfa_patching_gap_extension_score2;
    } else {
        wfa_convex_penalties.match = 0;
        wfa_convex_penalties.mismatch = 6;
        wfa_convex_penalties.gap_opening1 = 6;
        wfa_convex_penalties.gap_extension1 = 2;
        wfa_convex_penalties.gap_opening2 = 26;
        wfa_convex_penalties.gap_extension2 = 1;
    }

    // Use biWFA for smaller sequences or very high identity matches

    if (!force_wflign && (query_length <= MAX_LEN_FOR_STANDARD_WFA || target_length <= MAX_LEN_FOR_STANDARD_WFA)) {
        do_biwfa_alignment(
            query_name, query, query_total_length, query_offset, query_length, query_is_rev,
            target_name, target, target_total_length, target_offset, target_length,
            *out, wfa_convex_penalties, emit_md_tag, paf_format_else_sam, no_seq_in_sam, disable_chain_patching,
            min_identity, wflign_max_len_minor, mashmap_estimated_identity,
            -1, 1, 1); // Not part of a chain when using direct biWFA
        return;
    }

    // Check if mashmap_estimated_identity == 1 to avoid division by zero, leading to a minhash_kmer_size of 8.
    // Such low value was leading to confusion in HORs alignments in the human centromeres (high runtime and memory usage, and wrong alignments)
    const int minhash_kmer_size = mashmap_estimated_identity == 1 ? 17 : std::max(8, std::min(17, (int)std::floor(1.0 / (1.0 - mashmap_estimated_identity))));

    // Set penalties for wfa affine
    wflign_penalties_t wfa_affine_penalties;
    if (wfa_mismatch_score > 0 && wfa_gap_opening_score > 0 && wfa_gap_extension_score > 0){
        wfa_affine_penalties.match = 0;
        wfa_affine_penalties.mismatch = wfa_mismatch_score;
        wfa_affine_penalties.gap_opening1 = wfa_gap_opening_score;
        wfa_affine_penalties.gap_extension1 = wfa_gap_extension_score;
    } else {
        wfa_affine_penalties.match = 0;
        wfa_affine_penalties.mismatch = 6;
        wfa_affine_penalties.gap_opening1 = 8;
        wfa_affine_penalties.gap_extension1 = 1;
        /*
        if (mashmap_estimated_identity >= 0.80) {
            // Polynomial fitting
            const int mismatch = (int)ceil(9190.56599005*pow(mashmap_estimated_identity, 3) -24087.9418638*pow(mashmap_estimated_identity, 2) + 21032.49248734*mashmap_estimated_identity -6111.50339641);
            const int gap_opening = (int)ceil(11826.71109956*pow(mashmap_estimated_identity, 3) -30851.33099739 * pow(mashmap_estimated_identity, 2) + 26827.95391065*mashmap_estimated_identity -7767.00185348);

            wfa_affine_penalties.match = 0;
            wfa_affine_penalties.mismatch = mismatch;
            wfa_affine_penalties.gap_opening = gap_opening;
            wfa_affine_penalties.gap_extension = 1;
        } else {
            wfa_affine_penalties.match = 0;
            wfa_affine_penalties.mismatch = 3;
            wfa_affine_penalties.gap_opening = 5;
            wfa_affine_penalties.gap_extension = 1;
        }
        */
    }


    // heuristic bound on the max mash dist, adaptive based on estimated
    // identity the goal here is to sparsify the set of alignments in the
    // wflambda layer we then patch up the gaps between them

    float inception_score_max_ratio = 1.0 + 0.5 / mashmap_estimated_identity;
    float max_mash_dist_to_evaluate = std::min(0.55, 0.05 / std::pow(mashmap_estimated_identity,13));
    float mash_sketch_rate = 1.0;
    int wf_max_dist_threshold = 256;

    if (mashmap_estimated_identity >= 0.99) {
        mash_sketch_rate = 0.1;
    } else if (mashmap_estimated_identity >= 0.98) {
        mash_sketch_rate = 0.15;
    } else if (mashmap_estimated_identity >= 0.97) {
        mash_sketch_rate = 0.2;
    } else if (mashmap_estimated_identity >= 0.95) {
        mash_sketch_rate = 0.25;
    } else if (mashmap_estimated_identity >= 0.9) {
        mash_sketch_rate = 0.5;
    } /*else if (mashmap_estimated_identity >= 0.85) {
    } else if (mashmap_estimated_identity >= 0.8) {
    } else if (mashmap_estimated_identity >= 0.75) {
    } else {
    }*/

    // override erosion if not given on input
    if (erode_k < 0) {
        // heuristic setting of erosion
        erode_k = std::min(127.0,std::round(1.0/(1.0-mashmap_estimated_identity)));
    }

    // override max mash dist if given on input
    if (wflign_max_mash_dist > 0) {
        max_mash_dist_to_evaluate = wflign_max_mash_dist;
    }

    // accumulate runs of matches in reverse order
    // then trim the cigars of successive mappings
    std::vector<alignment_t*> trace;
#ifdef WFA_PNG_TSV_TIMING
    const auto start_time = std::chrono::steady_clock::now();
#endif

    if (!force_wflign) {
        wfa::WFAlignerGapAffine2Pieces* wf_aligner =
                new wfa::WFAlignerGapAffine2Pieces(
                        0,
                        wfa_convex_penalties.mismatch,
                        wfa_convex_penalties.gap_opening1,
                        wfa_convex_penalties.gap_extension1,
                        wfa_convex_penalties.gap_opening2,
                        wfa_convex_penalties.gap_extension2,
                        wfa::WFAligner::Alignment,
                        wfa::WFAligner::MemoryUltralow);
        wf_aligner->setHeuristicNone();
        
        const int status = wf_aligner->alignEnd2End(target,(int)target_length,query,(int)query_length);

        auto *aln = new alignment_t();
        aln->j = 0;
        aln->i = 0;

        aln->ok = (status == 0); // WF_ALIGN_SUCCESSFUL

        // fill the alignment info if we aligned
        if (aln->ok) {
            aln->query_length = query_length;
            aln->target_length = target_length;
    #ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(wf_aligner->cigar, query, target,
                        segment_length_q, segment_length_t, aln.j, aln.i)) {
                std::cerr << "cigar failure at alignment " << aln.j << " "
                << aln.i << std::endl;
                unpack_display_cigar(wf_aligner->cigar, query,
                                     target, segment_length_q, segment_length_t,
                                     aln.j, aln.i);
                std::cerr << ">query" << std::endl
                << std::string(query + j, segment_length_q)
                << std::endl;
                std::cerr << ">target" << std::endl
                << std::string(target + i, segment_length_t)
                << std::endl;
                assert(false);
            }
    #endif

            wflign_edit_cigar_copy(*wf_aligner,&aln->edit_cigar);

    #ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(aln.edit_cigar, query, target, segment_length_q,
                                segment_length_t, aln.j, aln.i)) {
                std::cerr << "cigar failure after cigar copy in alignment "
                << aln.j << " " << aln.i << std::endl;
                assert(false);
            }
    #endif
        }

        trace.push_back(aln);

#ifdef WFA_PNG_TSV_TIMING
        const long elapsed_time_wflambda_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - start_time).count();
#endif


        // Free old aligner
        delete wf_aligner;

        // use biWFA for all patching
        wf_aligner =
                new wfa::WFAlignerGapAffine2Pieces(
                        0,
                        wfa_convex_penalties.mismatch,
                        wfa_convex_penalties.gap_opening1,
                        wfa_convex_penalties.gap_extension1,
                        wfa_convex_penalties.gap_opening2,
                        wfa_convex_penalties.gap_extension2,
                        wfa::WFAligner::Alignment,
                        wfa::WFAligner::MemoryUltralow);
        wf_aligner->setHeuristicNone();

        // write a merged alignment
        write_merged_alignment(
                *out,
                trace,
                *wf_aligner,
                wfa_convex_penalties,
                emit_md_tag,
                paf_format_else_sam,
                no_seq_in_sam,
                query,
                query_name,
                query_total_length,
                query_offset,
                query_length,
                query_is_rev,
                target,
                target_name,
                target_total_length,
                target_offset,
                target_length,
                min_identity,
#ifdef WFA_PNG_TSV_TIMING
                elapsed_time_wflambda_ms,
                1,
                1,
#endif
                mashmap_estimated_identity,
                wflign_max_len_major,
                wflign_max_len_minor,
                erode_k,
                chain_gap,
                max_patching_score,
                min_inversion_length,
                MIN_WF_LENGTH,
                wf_max_dist_threshold
#ifdef WFA_PNG_TSV_TIMING
                ,
                prefix_wavefront_plot_in_png,
                wfplot_max_size,
                emit_patching_tsv,
                out_patching_tsv
#endif
                );

        // Free biWFA aligner
        delete wf_aligner;
    } else {
#ifdef WFA_PNG_TSV_TIMING
        if (emit_tsv) {
            *out_tsv << "# query_name=" << query_name << std::endl;
            *out_tsv << "# query_start=" << query_offset << std::endl;
            *out_tsv << "# query_end=" << query_offset+query_length << std::endl;
            *out_tsv << "# target_name=" << target_name << std::endl;
            *out_tsv << "# target_start=" << target_offset << std::endl;
            *out_tsv << "# target_end=" << target_offset+target_length << std::endl;
            *out_tsv << "# info: 0) mismatch, mash-distance > threshold; 1) mismatch, WFA-score >= max_score; 2) match, WFA-score < max_score" << std::endl;
            *out_tsv << "v" << "\t" << "h" << "\t" << "info" << std::endl;
        }
#endif

        const uint16_t segment_length_to_use =
                (query_length < segment_length || target_length < segment_length)
                ? std::min(query_length, target_length)
                : segment_length;

        // set up our implicit matrix
        const uint8_t steps_per_segment = 2;
        const uint16_t step_size = segment_length_to_use / steps_per_segment;

        // If the query_length/target_length are not multiple of step_size, we count
        // a fragment less, and the last one will be longer than segment_length_to_use
        const int pattern_length = (int)query_length / step_size - (query_length % step_size != 0 ? 1 : 0);
        const int text_length = (int)target_length / step_size - (target_length % step_size != 0 ? 1 : 0);

        wflign_penalties_t wflambda_affine_penalties;
        if (wflign_mismatch_score > 0 && wflign_gap_opening_score > 0 && wflign_gap_extension_score > 0){
            wflambda_affine_penalties.match = 0;
            wflambda_affine_penalties.mismatch = wflign_mismatch_score;
            wflambda_affine_penalties.gap_opening1 = wflign_gap_opening_score;
            wflambda_affine_penalties.gap_extension1 = wflign_gap_extension_score;
        } else {
            wflambda_affine_penalties.match = 0;
            wflambda_affine_penalties.mismatch = 4;
            wflambda_affine_penalties.gap_opening1 = 6;
            wflambda_affine_penalties.gap_extension1 = 1;
        }

        //std::cerr << "wfa_affine_penalties.mismatch " << wfa_affine_penalties.mismatch << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_opening " << wfa_affine_penalties.gap_opening << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_extension " << wfa_affine_penalties.gap_extension << std::endl;
        //std::cerr << "wflambda_affine_penalties.mismatch " << wflambda_affine_penalties.mismatch << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_opening1 " << wflambda_affine_penalties.gap_opening1 << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_extension1 " << wflambda_affine_penalties.gap_extension1 << std::endl;
        //std::cerr << "max_mash_dist_to_evaluate " << max_mash_dist_to_evaluate << std::endl;

        // Configure the attributes of the wflambda-aligner
        wfa::WFAlignerGapAffine* wflambda_aligner =
                new wfa::WFAlignerGapAffine(
                        wflambda_affine_penalties.mismatch,
                        wflambda_affine_penalties.gap_opening1,
                        wflambda_affine_penalties.gap_extension1,
                        wfa::WFAligner::Alignment,
                        wfa::WFAligner::MemoryUltralow);
        wflambda_aligner->setHeuristicNone(); // It should help
        if (wflign_max_distance_threshold <= 0) {
            wflambda_aligner->setHeuristicWFmash(wflign_min_wavefront_length, (int) (2048.0 / (mashmap_estimated_identity*mashmap_estimated_identity)));
        } else {
            wflambda_aligner->setHeuristicWFmash(wflign_min_wavefront_length, wflign_max_distance_threshold);
        }

        // Save computed alignments in a pair-indexed map
        robin_hood::unordered_flat_map<uint64_t,alignment_t*> alignments;
        // Allocate vectors to store our sketches
        std::vector<std::vector<rkmh::hash_t>*> query_sketches(pattern_length,nullptr);
        std::vector<std::vector<rkmh::hash_t>*> target_sketches(text_length,nullptr);

        // Allocate subsidiary WFAligner
        wfa::WFAlignerGapAffine* wf_aligner =
                new wfa::WFAlignerGapAffine(
                        wfa_affine_penalties.mismatch,
                        wfa_affine_penalties.gap_opening1,
                        wfa_affine_penalties.gap_extension1,
                        wfa::WFAligner::Alignment,
                        wfa::WFAligner::MemoryHigh);
        wf_aligner->setHeuristicNone();

        // Save mismatches if wfplots are requested
        robin_hood::unordered_set<uint64_t> high_order_dp_matrix_mismatch;

        // Setup WFling extend data
        wflign_extend_data_t extend_data;
        extend_data.wflign = this;
        extend_data.pattern_length = pattern_length;
        extend_data.text_length = text_length;
        extend_data.step_size = step_size;
        extend_data.segment_length_to_use = segment_length_to_use;
        extend_data.minhash_kmer_size = minhash_kmer_size;
        extend_data.max_mash_dist_to_evaluate = max_mash_dist_to_evaluate;
        extend_data.mash_sketch_rate = mash_sketch_rate;
        extend_data.inception_score_max_ratio = inception_score_max_ratio;
        extend_data.alignments = &alignments;
        extend_data.query_sketches = &query_sketches;
        extend_data.target_sketches = &target_sketches;
        extend_data.wf_aligner = wf_aligner;
//        extend_data.wflambda_aligner = wflambda_aligner;
//        extend_data.last_breakpoint_v = 0;
//        extend_data.last_breakpoint_h = 0;
//        extend_data.wfa_affine_penalties = &wfa_affine_penalties;
#ifdef WFA_PNG_TSV_TIMING
        extend_data.num_alignments = 0;
        extend_data.num_alignments_performed = 0;
#endif
        extend_data.num_sketches_allocated = 0;
        // 128 MB of memory for sketches
        extend_data.max_num_sketches_in_memory = 128 * 1024 * 1024
            / (sizeof(std::vector<rkmh::hash_t>) + mash_sketch_rate * segment_length_to_use * sizeof(rkmh::hash_t));
#ifdef WFA_PNG_TSV_TIMING
        extend_data.emit_png = !prefix_wavefront_plot_in_png->empty() && wfplot_max_size > 0;
        extend_data.high_order_dp_matrix_mismatch = &high_order_dp_matrix_mismatch;
#endif

        // Align
        wflambda_aligner->alignEnd2End(
                wflambda_extend_match, (void*)&extend_data,
                pattern_length,text_length);

        // Extract the trace
        if (wflambda_aligner->getAlignmentStatus() == WF_STATUS_ALG_COMPLETED) {
#ifdef WFA_PNG_TSV_TIMING
            extend_data.num_alignments += wflambda_trace_match(
                    alignments,*wflambda_aligner,trace,pattern_length,text_length);
#else
            wflambda_trace_match(alignments,*wflambda_aligner,trace,pattern_length,text_length);
#endif
        }

        // Free
        delete wflambda_aligner;
        delete wf_aligner;

#ifdef WFA_PNG_TSV_TIMING
        if (extend_data.emit_png) {
            const int wfplot_vmin = 0, wfplot_vmax = pattern_length; //v_max;
            const int wfplot_hmin = 0, wfplot_hmax = text_length; //h_max

            int v_max = wfplot_vmax - wfplot_vmin;
            int h_max = wfplot_hmax - wfplot_hmin;
            const algorithms::color_t COLOR_MASH_MISMATCH = { 0xffefefef };
            const algorithms::color_t COLOR_WFA_MISMATCH = { 0xff0000ff };
            const algorithms::color_t COLOR_WFA_MATCH = { 0xff00ff00 };

            const double scale = std::min(1.0, (double)wfplot_max_size / (double)std::max(v_max, h_max));
            const uint64_t width = (uint64_t)(scale * (double)v_max);
            const uint64_t height = (uint64_t)(scale * (double)h_max);
            const double source_width = (double)width;
            const double source_height = (double)height;

            const double x_off = 0, y_off = 0;
            //const double line_width = 1.0;
            const double source_min_x = 0, source_min_y = 0;

            auto plot_point = (v_max <= wfplot_max_size && h_max <= wfplot_max_size)
                              ? [](const algorithms::xy_d_t &point, algorithms::atomic_image_buf_t& image, const algorithms::color_t &color) {
                        image.layer_pixel(point.x, point.y, color);
                    }
                              : [](const algorithms::xy_d_t &point, algorithms::atomic_image_buf_t& image, const algorithms::color_t &color) {
                        wflign::algorithms::wu_calc_wide_line(
                                point, point,
                                color,
                                image);
                    };

            // Plot with only fragments belonging to the best alignment (that is, only the anchors)
            if (!trace.empty()) {
                algorithms::atomic_image_buf_t image(width, height,
                                                     source_width, source_height,
                                                     source_min_x, source_min_y);

                for (const auto &p : alignments) {
                    if (p.second != nullptr && p.second->keep) {
                        int v, h;
                        decode_pair(p.first, &v, &h);

                        if (v >= wfplot_vmin & v <= wfplot_vmax && h >= wfplot_hmin && h <= wfplot_hmax) {
                            algorithms::xy_d_t xy0 = {
                                    (v * scale) - x_off,
                                    (h * scale) + y_off
                            };
                            xy0.into(source_min_x, source_min_y,
                                     source_width, source_height,
                                     0, 0,
                                     width, height);

                            plot_point(xy0, image, COLOR_WFA_MATCH);
                        }
                    }
                }

                auto bytes = image.to_bytes();
                const std::string filename = *prefix_wavefront_plot_in_png +
                                             query_name + "_" + std::to_string(query_offset) + "_" + std::to_string(query_offset+query_length) + " _ " + (query_is_rev ? "-" : "+") +
                                             "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".1.anchors.png";
                encodeOneStep(filename.c_str(), bytes, width, height);
            }

            // Full plot
            algorithms::atomic_image_buf_t image(width, height,
                                                 source_width, source_height,
                                                 source_min_x, source_min_y);

            for (const auto &p : alignments) {
                int v, h;
                decode_pair(p.first, &v, &h);

                if (v >= wfplot_vmin & v <= wfplot_vmax && h >= wfplot_hmin && h <= wfplot_hmax) {
                    algorithms::xy_d_t xy0 = {
                            (v * scale) - x_off,
                            (h * scale) + y_off
                    };
                    xy0.into(source_min_x, source_min_y,
                             source_width, source_height,
                             0, 0,
                             width, height);

                    plot_point(xy0, image, p.second != nullptr ? COLOR_WFA_MATCH : COLOR_WFA_MISMATCH);
                }
            }

            for (auto high_order_DP_cell: high_order_dp_matrix_mismatch) {
                int v, h;
                decode_pair(high_order_DP_cell, &v, &h);

                if (v >= wfplot_vmin & v <= wfplot_vmax && h >= wfplot_hmin && h <= wfplot_hmax){
                    algorithms::xy_d_t xy0 = {
                            (v * scale) - x_off,
                            (h * scale) + y_off
                    };
                    xy0.into(source_min_x, source_min_y,
                             source_width, source_height,
                             0, 0,
                             width, height);

                    plot_point(xy0, image, COLOR_MASH_MISMATCH);
                }
            }

            auto bytes = image.to_bytes();
            const std::string filename = *prefix_wavefront_plot_in_png +
                                         query_name + "_" + std::to_string(query_offset) + "_" + std::to_string(query_offset+query_length) + " _ " + (query_is_rev ? "-" : "+") +
                                         "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".0.wflign.png";
            encodeOneStep(filename.c_str(), bytes, width, height);
        }
#endif

        // Clean alignments not to be kept (do not belong to the optimal alignment)
        for (const auto &p : alignments) {
            if (p.second != nullptr && !p.second->keep) {
                delete p.second;
                //p.second = nullptr;
            }
        }
#ifdef WFA_PNG_TSV_TIMING
        const long elapsed_time_wflambda_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now()-start_time).count();
#endif
        //#define WFLIGN_DEBUG
    #ifdef WFLIGN_DEBUG
        // get alignment score
        const int score = wflambda::edit_cigar_score_gap_affine(
            &affine_wavefronts->edit_cigar, &wflambda_affine_penalties);

        std::cerr << "[wflign::wflign_affine_wavefront] alignment score " << score
        << " for query: " << query_name << " target: " << target_name
        << std::endl;
    #endif

        clean_up_sketches(query_sketches);
        clean_up_sketches(target_sketches);

        // todo: implement alignment identifier based on hash of the input, params,
        // and commit annotate each PAF record with it and the full alignment score

        // Trim alignments that overlap in the query
        if (!trace.empty()) {
    #ifdef VALIDATE_WFA_WFLIGN
            if (!trace.front()->validate(query, target)) {
                std::cerr << "first traceback is wrong" << std::endl;
                trace.front()->display();
                assert(false);
            }
    #endif

            auto x = trace.rbegin();
            auto end_trace = trace.rend();

            auto c = x + 1;
            while (x != end_trace && c != end_trace) {
                // establish our last and curr alignments to consider when trimming
                auto &last = **x;
                auto &curr = **c;

    #ifdef VALIDATE_WFA_WFLIGN
                if (curr.ok && !curr.validate(query, target)) {
                    std::cerr << "curr traceback is wrong before trimming @ "
                    << curr.j << " " << curr.i << std::endl;
                    curr.display();
                    assert(false);
                }
    #endif
                trace_pos_t last_pos(
                        last.j,last.i,&last.edit_cigar,
                        last.edit_cigar.begin_offset);
                trace_pos_t curr_pos(
                        curr.j,curr.i,&curr.edit_cigar,
                        curr.edit_cigar.begin_offset);

                // trace the last alignment until we overlap the next
                // to record our match
                trace_pos_t match_pos;

                // walk until they are matched at the query position
                while (!last_pos.at_end() && !curr_pos.at_end()) {
                    //std::cerr << "last_pos " << last_pos.j << " - " << last_pos.i << " - " << last_pos.offset << " - " << last_pos.curr() << std::endl;
                    //std::cerr << "curr_pos " << curr_pos.j << " - " << curr_pos.i << " - " << curr_pos.offset << " - " << curr_pos.curr() << std::endl;
                    if (last_pos.equal(curr_pos)) {
                        // they equal and we can splice them at the first match
                        match_pos = last_pos;
                        break;
                    }
                    if (last_pos.j == curr_pos.j) {
                        last_pos.incr();
                        curr_pos.incr();
                    } else if (last_pos.j < curr_pos.j) {
                        last_pos.incr();
                    } else {
                        curr_pos.incr();
                    }
                }

                // if we matched, we'll be able to splice the alignments together
                int trim_last, trim_curr;
                if (match_pos.assigned()) {
                    //std::cerr << "match_pos " << match_pos.j << " - " << match_pos.i << " - " << match_pos.offset << " - " << match_pos.curr() << std::endl;
                    // we'll use our match position to set up the trims
                    trim_last = (last.j + last.query_length) - match_pos.j;
                    trim_curr = match_pos.j - curr.j;
                } else {
                    // we want to remove any possible overlaps in query or target
                    // walk back last until we don't overlap in i or j
                    // recording the distance walked as an additional trim on last
                    bool flip = false;
                    while (last_pos.j > curr_pos.j || last_pos.i > curr_pos.i) {
                        if (flip) {
                            last_pos.decr();
                        } else {
                            curr_pos.incr();
                        }
                        flip ^= true;
                    }
                    trim_last = (last.j + last.query_length) - last_pos.j + 1;
                    trim_curr = curr_pos.j - curr.j + 1;
                    assert(last_pos.j <= curr_pos.j);
                    assert(last_pos.i <= curr_pos.i);
                }

                // assign our cigar trim
                //last.display();
                if (trim_last > 0) {
                    //std::cerr << "trim_last " << trim_last << std::endl;
                    last.trim_back(trim_last);
                    //last.display();
    #ifdef VALIDATE_WFA_WFLIGN
                    if (last.ok && !last.validate(query, target)) {
                        std::cerr << "traceback is wrong after last trimming @ "
                        << last.j << " " << last.i << std::endl;
                        last.display();
                        assert(false);
                    }
    #endif
                }
                //curr.display();
                if (trim_curr > 0) {
                    //std::cerr << "trim_curr " << trim_curr << std::endl;
                    curr.trim_front(trim_curr);
                    //curr.display();
    #ifdef VALIDATE_WFA_WFLIGN
                    if (curr.ok && !curr.validate(query, target)) {
                        std::cerr << "traceback is wrong after curr trimming @ "
                        << curr.j << " " << curr.i << std::endl;
                        curr.display();
                        assert(false);
                    }
    #endif
                }
                if (curr.ok) {
                    x = c;
                }
                ++c;
    #ifdef VALIDATE_WFA_WFLIGN
                auto distance_target = (curr.i - (last.i + last.target_length));
                auto distance_query = (curr.j - (last.j + last.query_length));
                if (last.ok && curr.ok &&
                (distance_query < 0 || distance_target < 0)) {
                    std::cerr << "distance_target_query " << distance_target << " "
                    << distance_query << std::endl;
                    std::cerr << "trimming failure at @ " << last.j << "," << last.i
                    << " -> " << curr.j << "," << curr.i << std::endl;
                    last.display();
                    curr.display();
                    exit(1);
                }
    #endif
            }

            if (merge_alignments) {
                // use biWFA for all patching
                wfa::WFAlignerGapAffine2Pieces* wf_aligner =
                        new wfa::WFAlignerGapAffine2Pieces(
                                0,
                                wfa_convex_penalties.mismatch,
                                wfa_convex_penalties.gap_opening1,
                                wfa_convex_penalties.gap_extension1,
                                wfa_convex_penalties.gap_opening2,
                                wfa_convex_penalties.gap_extension2,
                                wfa::WFAligner::Alignment,
                                wfa::WFAligner::MemoryUltralow);
                wf_aligner->setHeuristicNone();

                // write a merged alignment
                write_merged_alignment(
                        *out,
                        trace,
                        *wf_aligner,
                        wfa_convex_penalties,
                        emit_md_tag,
                        paf_format_else_sam,
                        no_seq_in_sam,
                        query,
                        query_name,
                        query_total_length,
                        query_offset,
                        query_length,
                        query_is_rev,
                        target,
                        target_name,
                        target_total_length,
                        target_offset,
                        target_length,
                        min_identity,
#ifdef WFA_PNG_TSV_TIMING
                        elapsed_time_wflambda_ms,
                        extend_data.num_alignments,
                        extend_data.num_alignments_performed,
#endif
                        mashmap_estimated_identity,
                        wflign_max_len_major,
                        wflign_max_len_minor,
                        erode_k,
                        chain_gap,
                        max_patching_score,
                        min_inversion_length,
                        MIN_WF_LENGTH,
                        wf_max_dist_threshold
#ifdef WFA_PNG_TSV_TIMING
                        ,
                        prefix_wavefront_plot_in_png,
                        wfplot_max_size,
                        emit_patching_tsv,
                        out_patching_tsv
#endif
                );

                delete wf_aligner;
            } else {
                // todo old implementation (and SAM format is not supported)
                for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
                    write_alignment_paf(
                            *out,
                            **x,
                            "",
                            query_name,
                            query_total_length,
                            query_offset,
                            query_length,
                            query_is_rev,
                            target_name,
                            target_total_length,
                            target_offset,
                            target_length,
                            min_identity,
                            mashmap_estimated_identity,
                            0, 0, 0);
                }
            }
        }
    }
}

} // namespace wavefront
} // namespace wflign
