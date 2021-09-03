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

//#include "WFA/gap_affine/affine_wavefront.hpp"
//#include "WFA/gap_affine/affine_wavefront_align.h"
#include "WFA/gap_affine/affine_matrix.h"
#include "WFA/gap_affine/swg.h"
#include "WFA/wavefront/wavefront_align.h"
#include "WFA/wavefront/wavefront_reduction.h"
#include "edlib.h"
#include "robin-hood-hashing/robin_hood.h"

//#include "wfa_edit_callback.hpp"
#include "dna.hpp"
#include "rkmh.hpp"
#include "wflambda/utils/commons.h"
#include "wflambda/gap_affine/affine_matrix.h"
#include "wflambda/gap_affine/swg.h"
#include "wflambda/wavefront/wavefront_align.h"
#include "wflambda/wavefront/wavefront_reduction.h"

//#define WFLIGN_DEBUG true // for debugging messages

namespace wflign {

namespace wavefront {

/*bool hack_cigar(wfa::cigar_t &cigar, const char *query, const char *target,
                const uint64_t &query_aln_len, const uint64_t &target_aln_len,
                uint64_t j, uint64_t i);*/

bool validate_cigar(const wfa::cigar_t &cigar, const char *query,
                    const char *target, const uint64_t &query_aln_len,
                    const uint64_t &target_aln_len, uint64_t j, uint64_t i);

bool validate_trace(const std::vector<char> &tracev, const char *query,
                    const char *target, const uint64_t &query_aln_len,
                    const uint64_t &target_aln_len, uint64_t j, uint64_t i);

bool unpack_display_cigar(const wfa::cigar_t &cigar, const char *query,
                          const char *target, const uint64_t &query_aln_len,
                          const uint64_t &target_aln_len, uint64_t j,
                          uint64_t i);

struct alignment_t {
    int j = 0;
    int i = 0;
    uint16_t query_length = 0;
    uint16_t target_length = 0;
    bool ok = false;
    bool keep = false;
    //int score = std::numeric_limits<int>::max();
    // float mash_dist = 1;
    wfa::cigar_t edit_cigar{};
    void display(void) {
        std::cerr << j << " " << i << " " << query_length << " "
                  << target_length << " " << ok << std::endl;
        for (int x = edit_cigar.begin_offset; x < edit_cigar.end_offset; ++x) {
            std::cerr << edit_cigar.operations[x];
        }
        std::cerr << std::endl;
    }
    bool validate(const char *query, const char *target) {
        return validate_cigar(edit_cigar, query, target, query_length,
                              target_length, j, i);
    }
    void trim_front(int query_trim) {
        // this kills the alignment
        if (query_trim >= query_length) {
            ok = false;
            return;
        }
        // increment j and i appropriately
        int trim_to_j = j + query_trim;
        int x = edit_cigar.begin_offset;
        while (x < edit_cigar.end_offset && j < trim_to_j) {
            switch (edit_cigar.operations[x++]) {
            case 'M':
            case 'X':
                --query_length;
                --target_length;
                ++j;
                ++i;
                break;
            case 'I':
                --query_length;
                ++j;
                break;
            case 'D':
                --target_length;
                ++i;
                break;
            default:
                break;
            }
            if (target_length <= 0 || query_length <= 0) {
                ok = false;
                return;
            }
        }
        while (x < edit_cigar.end_offset && edit_cigar.operations[x] == 'D') {
            ++x;
            --target_length;
            ++i;
        }
        if (x == edit_cigar.end_offset)
            ok = false;
        edit_cigar.begin_offset = x;
    }
    void trim_back(int query_trim) {
        if (query_trim >= query_length) {
            ok = false;
            return;
        }
        int x = edit_cigar.end_offset;
        int q = 0;
        while (x > edit_cigar.begin_offset && q < query_trim) {
            switch (edit_cigar.operations[--x]) {
            case 'M':
            case 'X':
                --query_length;
                --target_length;
                ++q;
                break;
            case 'I':
                --query_length;
                ++q;
                break;
            case 'D':
                --target_length;
                break;
            default:
                break;
            }
            if (target_length <= 0 || query_length <= 0) {
                ok = false;
                return;
            }
        }
        while (x >= edit_cigar.begin_offset &&
               edit_cigar.operations[x - 1] == 'D') {
            --x;
            --target_length;
        }
        if (x == edit_cigar.begin_offset)
            ok = false;
        edit_cigar.end_offset = x;
    }
    ~alignment_t() { free(edit_cigar.operations); }
};

// link a position in a traceback matrix to its edit
struct trace_pos_t {
    int j = 0;
    int i = 0;
    wfa::cigar_t *edit_cigar = nullptr;
    int offset = 0;
    bool incr() {
        if (offset < edit_cigar->end_offset) {
            switch (curr()) {
            case 'M':
            case 'X':
                ++j;
                ++i;
                break;
            case 'I':
                ++j;
                break;
            case 'D':
                ++i;
                break;
            default:
                break;
            }
            ++offset;
            return true;
        } else {
            return false;
        }
    }
    bool decr() {
        if (offset > edit_cigar->begin_offset) {
            --offset;
            switch (curr()) {
            case 'M':
            case 'X':
                --j;
                --i;
                break;
            case 'I':
                --j;
                break;
            case 'D':
                --i;
                break;
            default:
                break;
            }
            return true;
        } else {
            return false;
        }
    }
    bool at_end() const { return offset == edit_cigar->end_offset; }
    char curr() const {
        assert(!at_end());
        return edit_cigar->operations[offset];
    }
    bool equal(const trace_pos_t &other) const {
        return j == other.j && i == other.i && curr() == 'M' &&
               curr() == other.curr();
    }
    bool assigned() const { return edit_cigar != nullptr; }
};

void wflign_edit_cigar_copy(wfa::cigar_t *const edit_cigar_dst,
                            wfa::cigar_t *const edit_cigar_src);

void copy_wfa_alignment_into_trace(const wfa::cigar_t *const edit_cigar,
                                   std::vector<char> &trace);

/*
void edlib_to_wflign_edit_cigar_copy(
    wfa::cigar_t* const edit_cigar_dst,
    char* const edlib_cigar_src,
    const uint64_t& edit_distance,
    const uint64_t& edlib_cigar_len);
*/

inline uint64_t encode_pair(int v, int h) {
    return ((uint64_t)v << 32) | (uint64_t)h;
}

wfa::wavefront_aligner_t* get_wavefront_aligner(
    const wfa::affine_penalties_t& wfa_affine_penalties,
    const uint64_t& target_length,
    const uint64_t& query_length,
    const bool& low_memory);

void wflign_affine_wavefront(
    std::ostream &out,
    const bool &emit_tsv, std::ostream &out_tsv,
    const bool &merge_alignments, const bool &emit_md_tag,
    const bool &paf_format_else_sam, const std::string &query_name,
    const char *query, const uint64_t &query_total_length,
    const uint64_t &query_offset, const uint64_t &query_length,
    const bool &query_is_rev, const std::string &target_name,
    const char *target, const uint64_t &target_total_length,
    const uint64_t &target_offset, const uint64_t &target_length,
    const uint16_t &segment_length, const float &min_identity, const int& minhash_kmer_size,
    const int &wfa_mismatch_score,
    const int &wfa_gap_opening_score,
    const int &wfa_gap_extension_score,
    const int &wflambda_min_wavefront_length, // with these set at 0 we do exact
                                              // WFA for wflambda
    const int &wflambda_max_distance_threshold, const float &mashmap_estimated_identity,
    const int &wflign_mismatch_score,
    const int &wflign_gap_opening_score,
    const int &wflign_gap_extension_score,
    const float &wflign_max_mash_dist,
    const uint64_t &wflign_max_len_major, const uint64_t &wflign_max_len_minor,
    const uint16_t &erode_k);
// const int& wfa_min_wavefront_length, // with these set at 0 we do exact WFA
// for WFA itself const int& wfa_max_distance_threshold);

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

//EdlibAlignResult do_edlib_patch_alignment(const char *query, const uint64_t &j,
//                                          const uint64_t &query_length,
//                                          const char *target, const uint64_t &i,
//                                          const uint64_t &target_length,
//                                          const EdlibAlignMode &align_mode);

void write_merged_alignment(
    std::ostream &out, const std::vector<alignment_t *> &trace,
    wfa::wavefront_aligner_t *const wf_aligner,
    wfa::affine_penalties_t *const affine_penalties, const bool &emit_md_tag,
    const bool &paf_format_else_sam, const char *query,
    const std::string &query_name, const uint64_t &query_total_length,
    const uint64_t &query_offset, const uint64_t &query_length,
    const bool &query_is_rev, const char *target,
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

/*char* edlib_alignment_to_cigar(
    const unsigned char* const alignment,
    const int alignment_length,
    uint64_t& target_aligned_length,
    uint64_t& query_aligned_length,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& inserted_bp,
    uint64_t& deletions,
    uint64_t& deleted_bp);*/

double float2phred(const double &prob);

void sort_indels(std::vector<char> &v);

} // namespace wavefront

} // namespace wflign
