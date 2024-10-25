#ifndef ALIGNMENT_PRINTER_HPP
#define ALIGNMENT_PRINTER_HPP

#include <string>
#include <sstream>
#include "wflign.hpp"

namespace wflign {
namespace wavefront {

void write_tag_and_md_string(
    std::ostream &out,
    const char *cigar_ops,
    const int cigar_start,
    const int cigar_end,
    const int target_start,
    const char *target,
    const int64_t target_offset,
    const int64_t target_pointer_shift);

void write_alignment_sam(
    std::ostream &out,
    const alignment_t& patch_aln,
    const std::string& query_name,
    const uint64_t& query_total_length,
    const uint64_t& query_offset,
    const uint64_t& query_length,
    const bool& query_is_rev,
    const std::string& target_name,
    const uint64_t& target_total_length,
    const uint64_t& target_offset,
    const uint64_t& target_length,
    const float& min_identity,
    const float& mashmap_estimated_identity,
    const bool& no_seq_in_sam,
    const bool& emit_md_tag,
    const char* query,
    const char* target,
    const int64_t& target_pointer_shift);

bool write_alignment_paf(
    std::ostream& out,
    const alignment_t& aln,
    const std::string& query_name,
    const uint64_t& query_total_length,
    const uint64_t& query_offset,
    const uint64_t& query_length,
    const bool& query_is_rev,
    const std::string& target_name,
    const uint64_t& target_total_length,
    const uint64_t& target_offset,
    const uint64_t& target_length,
    const float& min_identity,
    const float& mashmap_estimated_identity,
    const bool& with_endline = true,
    const bool& is_rev_patch = false);

void write_merged_alignment(
    std::ostream &out,
    const std::vector<alignment_t *> &trace,
    wfa::WFAlignerGapAffine2Pieces& wf_aligner,
    const wflign_penalties_t& convex_penalties,
    const bool& emit_md_tag,
    const bool& paf_format_else_sam,
    const bool& no_seq_in_sam,
    const char* query,
    const std::string& query_name,
    const uint64_t& query_total_length,
    const uint64_t& query_offset,
    const uint64_t& query_length,
    const bool& query_is_rev,
    const char* target,
    const std::string& target_name,
    const uint64_t& target_total_length,
    const uint64_t& target_offset,
    const uint64_t& target_length,
    const float& min_identity,
#ifdef WFA_PNG_TSV_TIMING
    const long& elapsed_time_wflambda_ms,
    const uint64_t& num_alignments,
    const uint64_t& num_alignments_performed,
#endif
    const float& mashmap_estimated_identity,
    const uint64_t& wflign_max_len_major,
    const uint64_t& wflign_max_len_minor,
    const int& erode_k,
    const int64_t& chain_gap,
    const int& max_patching_score,
    const uint64_t& min_inversion_length,
    const int& min_wf_length,
    const int& max_dist_threshold,
#ifdef WFA_PNG_TSV_TIMING
    const std::string* prefix_wavefront_plot_in_png,
    const uint64_t& wfplot_max_size,
    const bool& emit_patching_tsv,
    std::ostream* out_patching_tsv,
#endif
    const bool& with_endline = true);

} // namespace wavefront
} // namespace wflign

#endif
