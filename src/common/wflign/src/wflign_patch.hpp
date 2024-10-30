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
#include "alignment_printer.hpp"

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
            const uint64_t& min_inversion_length);

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

        double float2phred(const double& prob);
        void sort_indels(std::vector<char>& v);

    } /* namespace wavefront */

    void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height);
} /* namespace wflign */

#endif /* WFLIGN_PATCH_HPP_ */
