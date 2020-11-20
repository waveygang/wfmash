#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include "edlib.h"
#include "patchmap.hpp"
#include "wfa_edit_callback.hpp"
#include "wflambda/gap_affine/affine_wavefront_align.hpp"
#include "wflambda/gap_affine/affine_wavefront_backtrace.hpp"

//#define WFLIGN_DEBUG true // for debugging messages

namespace wflign {

struct alignment_t {
    int j = 0;
    int i = 0;
    double identity = 0;
    int skip_query_start = 0;
    int keep_query_length = 0;
    int skip_target_start = 0;
    int keep_target_length = 0;
    const std::string* query_name;
    uint64_t query_size;
    const std::string* target_name;
    uint64_t target_size;
    EdlibAlignResult result;
    ~alignment_t(void) {
        edlibFreeAlignResult(result);
    }
};

inline uint64_t encode_pair(int v, int h) {
    return ((uint64_t)v << 32) | (uint64_t)h;
}

void wflign_full(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length);

void wflign_wavefront(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length);

void wflign_affine_wavefront(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length,
    const float& min_identity,
    const int& min_wavefront_length = 0, // with these set at 0 we do exact WFA
    const int& max_distance_threshold = 0);

bool do_alignment(
    const std::string& query_name,
    const std::string& query,
    const uint64_t& j,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& i,
    const uint64_t& segment_length,
    const uint64_t& step_size,
    alignment_t& alignment);

std::ostream& operator<<(
    std::ostream& out,
    const alignment_t& aln);

void write_alignment(
    std::ostream& out,
    const alignment_t& aln,
    const float& min_identity,
    const bool& with_endline = true);

char* alignmentToCigar(
    const unsigned char* const alignment,
    const int alignmentLength,
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
