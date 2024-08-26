#ifndef WFLIGN_PATCH_HPP_
#define WFLIGN_PATCH_HPP_

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>
#include <functional>
#include <fstream>
#include <sstream>
#include <lodepng/lodepng.h>

#include "dna.hpp"
#include "rkmh.hpp"
#include "wflign.hpp"
#include "wflign_alignment.hpp"

/*
 * Configuration
 */
//#define WFLIGN_DEBUG true // for debugging messages

/*
 * Namespaces
 */
namespace wflign {
    namespace wavefront {

        bool do_wfa_segment_alignment(
                const std::string& query_name,
                const char* query,
                std::vector<rkmh::hash_t>*& query_sketch,
                const uint64_t& query_length,
                const int64_t& j,
                const std::string& target_name,
                const char* target,
                std::vector<rkmh::hash_t>*& target_sketch,
                const uint64_t& target_length,
                const int64_t& i,
                const uint16_t& segment_length_q,
                const uint16_t& segment_length_t,
                const uint16_t& step_size,
                wflign_extend_data_t* extend_data,
                alignment_t& aln);
        void do_wfa_patch_alignment(
            const char* query,
            const uint64_t& j,
            const uint64_t& query_length,
            const char* target,
            const uint64_t& i,
            const uint64_t& target_length,
            wfa::WFAlignerGapAffine2Pieces& wf_aligner,
            const wflign_penalties_t& convex_penalties,
            alignment_t& aln,
            alignment_t& rev_aln,
            const int64_t& chain_gap,
            const int& max_patching_score,
            const uint64_t& min_inversion_length,
            bool ends_free);
        void trim_alignment(alignment_t& aln);
        std::vector<alignment_t> do_progressive_wfa_patch_alignment(
            const char* query,
            const uint64_t& query_start,
            const uint64_t& query_length,
            const char* target,
            const uint64_t& target_start,
            const uint64_t& target_length,
            wfa::WFAlignerGapAffine2Pieces& wf_aligner,
            const wflign_penalties_t& convex_penalties,
            const int64_t& chain_gap,
            const int& max_patching_score,
            const uint64_t& min_inversion_length,
            const int& erode_k);
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
        void write_alignment_paf(
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
                const uint64_t& target_length, // unused
                const float& min_identity,
                const float& mashmap_estimated_identity,
                const bool& with_endline = true,
                const bool& is_rev_patch = false);
        double float2phred(const double& prob);
        void sort_indels(std::vector<char>& v);

    } /* namespace wavefront */

    void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height);
} /* namespace wflign */

#endif /* WFLIGN_PATCH_HPP_ */
