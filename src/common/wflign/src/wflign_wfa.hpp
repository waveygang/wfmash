#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cctype>
#include <charconv>
#include "WFA/edit/edit_cigar.hpp"
//#include "WFA/gap_affine/affine_wavefront.hpp"
#include "WFA/gap_affine/affine_wavefront_align.hpp"
#include "patchmap.hpp"
//#include "wfa_edit_callback.hpp"
#include "wflambda/gap_affine/affine_wavefront_align.hpp"
#include "wflambda/gap_affine/affine_wavefront_backtrace.hpp"
#include "dna.hpp"
#include "rkmh.hpp"

//#define WFLIGN_DEBUG true // for debugging messages

namespace wflign {

namespace wavefront {

struct alignment_t {
    int j = 0;
    int i = 0;
    int length = 0;
    bool ok = false;
    int score = std::numeric_limits<int>::max();
    double mash_dist = 1;
    wfa::edit_cigar_t edit_cigar;
    void trim_front(int query_trim) {
        // increment j and i appropriately
        int trim_to_j = j + query_trim;
        int x = edit_cigar.begin_offset;
        while (x < edit_cigar.end_offset
               && j < trim_to_j) {
            switch (edit_cigar.operations[x++]) {
            case 'M': case 'X': --length; ++j; ++i; break;
            case 'I': --length; ++j; break;
            case 'D': ++i; break;
            default: break;
            }
        }
        if (x == edit_cigar.end_offset) ok = false;
        edit_cigar.begin_offset = x;
    }
    void trim_back(int query_trim) {
        int x = edit_cigar.end_offset;
        int q = 0;
        while (x > edit_cigar.begin_offset
               && q < query_trim) {
            switch (edit_cigar.operations[--x]) {
            case 'M': case 'X': --length; ++q; break;
            case 'I': --length; ++q; break;
            case 'D': break;
            default: break;
            }
        }
        if (x == edit_cigar.begin_offset) ok = false;
        edit_cigar.end_offset = x;
    }
    ~alignment_t(void) {
        free(edit_cigar.operations);
    }
};

// link a position in a traceback matrix to its edit
struct trace_pos_t {
    int j = 0;
    int i = 0;
    wfa::edit_cigar_t* edit_cigar = nullptr;
    int offset = 0;
    bool incr(void) {
        if (offset < edit_cigar->end_offset) {
            switch (curr()) {
            case 'M': case 'X': ++j; ++i; break;
            case 'I': ++j; break;
            case 'D': ++i; break;
            default: break;
            }
            ++offset;
            return true;
        } else {
            return false;
        }
    }
    bool decr(void) {
        if (offset > edit_cigar->begin_offset) {
            --offset;
            switch (curr()) {
            case 'M': case 'X': --j; --i; break;
            case 'I': --j; break;
            case 'D': --i; break;
            default: break;
            }
            return true;
        } else {
            return false;
        }
    }
    bool at_end(void) const {
        return offset == edit_cigar->end_offset;
    }
    char curr(void) const {
        assert(!at_end());
        return edit_cigar->operations[offset];
    }
    bool equal(const trace_pos_t& other) const {
        return j == other.j
            && i == other.i
            && curr() == 'M' && curr() == other.curr();
    }
    bool assigned(void) const {
        return edit_cigar != nullptr;
    }
};

void wflign_edit_cigar_copy(
    wfa::edit_cigar_t* const edit_cigar_dst,
    wfa::edit_cigar_t* const edit_cigar_src);

inline uint64_t encode_pair(int v, int h) {
    return ((uint64_t)v << 32) | (uint64_t)h;
}

void wflign_affine_wavefront(
    std::ostream& out,
    const bool& merge_alignments,
    const std::string& query_name,
    const char* query,
    const uint64_t& query_total_length,
    const uint64_t& query_offset,
    const uint64_t& query_length,
    const bool& query_is_rev,
    const std::string& target_name,
    const char* target,
    const uint64_t& target_total_length,
    const uint64_t& target_offset,
    const uint64_t& target_length,
    const uint64_t& segment_length,
    const float& min_identity,
    const int& wflambda_min_wavefront_length, // with these set at 0 we do exact WFA for wflambda
    const int& wflambda_max_distance_threshold);
    //const int& wfa_min_wavefront_length, // with these set at 0 we do exact WFA for WFA itself
    //const int& wfa_max_distance_threshold);

bool do_alignment(
    const std::string& query_name,
    const char* query,
    std::vector<rkmh::hash_t>*& query_sketches,
    const uint64_t& query_length,
    const uint64_t& j,
    const std::string& target_name,
    const char* target,
    std::vector<rkmh::hash_t>*& target_sketches,
    const uint64_t& target_length,
    const uint64_t& i,
    const uint64_t& segment_length,
    const uint64_t& step_size,
    const uint64_t& minhash_kmer_size,
    const int min_wavefront_length,
    const int max_distance_threshold,
    wfa::mm_allocator_t* const mm_allocator,
    wfa::affine_penalties_t* const affine_penalties,
    const float& min_identity,
    alignment_t& aln);

void merge_alignments(
    alignment_t& base,
    const alignment_t& ext);

void write_merged_alignment(
    std::ostream& out,
    const std::vector<alignment_t*> trace,
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
    const bool& with_endline = true);

void write_alignment(
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
    const bool& with_endline = true);

char* alignmentToCigar(
    const wfa::edit_cigar_t* const edit_cigar,
    uint64_t& refAlignedLength,
    uint64_t& qAlignedLength,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& deletions,
    uint64_t& softclips);

double float2phred(const double& prob);

}

}
