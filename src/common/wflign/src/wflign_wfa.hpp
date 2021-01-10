#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
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
    bool ok = false;
    int score = std::numeric_limits<int>::max();
    double mash_dist = 1;
    double identity = 0;
    int skip_query_start = 0;
    int keep_query_length = 0;
    int skip_target_start = 0;
    int keep_target_length = 0;
    wfa::edit_cigar_t edit_cigar;
    ~alignment_t(void) {
        free(edit_cigar.operations);
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
    const int skip_query_start,
    const int keep_query_length,
    int& skipped_target_start,
    int& kept_target_length,
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
