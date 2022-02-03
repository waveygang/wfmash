#include "WFA/edit/edit_dp.h"
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

#include "robin-hood-hashing/robin_hood.h"

//#include "wfa_edit_callback.hpp"
#include "dna.hpp"
#include "rkmh.hpp"

//#define WFLIGN_DEBUG true // for debugging messages

namespace wflign {

namespace wavefront {

bool do_wfa_segment_alignment(
    const std::string &query_name, const char *query,
    std::vector<rkmh::hash_t> *&query_sketches, const uint64_t &query_length,
    const uint64_t &j, const std::string &target_name, const char *target,
    std::vector<rkmh::hash_t> *&target_sketches, const uint64_t &target_length,
    const uint64_t &i,
    const uint16_t &segment_length_q,
    const uint16_t &segment_length_t,
    const uint16_t &step_size, const uint64_t &minhash_kmer_size,
    const int &min_wavefront_length,
    const int &max_distance_threshold, const float &max_mash_dist, const float& mashmap_estimated_identity,
    wfa::wavefront_aligner_t *const wf_aligner,
    wfa::affine_penalties_t *const affine_penalties, alignment_t &aln);

void do_wfa_patch_alignment(const char *query, const uint64_t &j,
                            const uint64_t &query_length, const char *target,
                            const uint64_t &i, const uint64_t &target_length,
                            const int &segment_length,
                            const int &min_wavefront_length,
                            const int &max_distance_threshold,
                            wfa::wavefront_aligner_t *const wf_aligner,
                            wfa::affine_penalties_t *const affine_penalties,
                            alignment_t &aln);

void write_merged_alignment(
    std::ostream &out, const std::vector<alignment_t *> &trace,
    wfa::wavefront_aligner_t *const wf_aligner,
    wfa::affine_penalties_t *const affine_penalties,
    const bool &emit_md_tag,
    const bool &paf_format_else_sam, const bool &seq_in_sam,
    const char *query,
    const std::string &query_name, const uint64_t &query_total_length,
    const uint64_t &query_offset, const uint64_t &query_length,
    const bool &query_is_rev,
    const char *target,
    const std::string &target_name, const uint64_t &target_total_length,
    const uint64_t &target_offset, const uint64_t &target_length,
    const uint16_t &segment_length,
    const float &min_identity, const long &elapsed_time_wflambda_ms,
    const uint64_t &num_alignments, const uint64_t &num_alignments_performed,
    const float &mashmap_estimated_identity,
    const uint64_t &wflign_max_len_major, const uint64_t &wflign_max_len_minor,
    const uint16_t &erode_k,
    const int &min_wf_length, const int &max_dist_threshold,
    const bool &with_endline = true);

void write_alignment(std::ostream &out, const alignment_t &aln,
                     const std::string &query_name,
                     const uint64_t &query_total_length,
                     const uint64_t &query_offset, const uint64_t &query_length,
                     const bool &query_is_rev, const std::string &target_name,
                     const uint64_t &target_total_length,
                     const uint64_t &target_offset,
                     const uint64_t &target_length, const float &min_identity,
                     const float &mashmap_estimated_identity,
                     const bool &with_endline = true);

char *alignment_to_cigar(const std::vector<char> &edit_cigar,
                         const uint64_t &start_idx, const uint64_t &end_idx,
                         uint64_t &target_aligned_length,
                         uint64_t &query_aligned_length, uint64_t &matches,
                         uint64_t &mismatches, uint64_t &insertions,
                         uint64_t &inserted_bp, uint64_t &deletions,
                         uint64_t &deleted_bp);

char *wfa_alignment_to_cigar(const wfa::cigar_t *const edit_cigar,
                             uint64_t &target_aligned_length,
                             uint64_t &query_aligned_length, uint64_t &matches,
                             uint64_t &mismatches, uint64_t &insertions,
                             uint64_t &inserted_bp, uint64_t &deletions,
                             uint64_t &deleted_bp);

double float2phred(const double &prob);

void sort_indels(std::vector<char> &v);

} // namespace wavefront

} // namespace wflign
