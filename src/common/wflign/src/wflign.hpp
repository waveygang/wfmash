#ifndef WFLIGN_HPP_
#define WFLIGN_HPP_

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

#include "wflign_alignment.hpp"

#include "robin-hood-hashing/robin_hood.h"
#include "dna.hpp"
#include "rkmh.hpp"

/*
 * Configuration
 */
//#define WFLIGN_DEBUG true // for debugging messages
//#define VALIDATE_WFA_WFLIGN

#include "atomic_image.hpp"
#include "lodepng/lodepng.h"


/*
 * Namespaces
 */
namespace wflign {
    namespace wavefront {

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
            const float min_identity,
            const uint64_t wflign_max_len_minor,
            const float mashmap_estimated_identity,
            const struct align::MappingBoundaryRow& currentRecord);

        class WFlign {
        public:
            // WFlambda parameters
            uint16_t segment_length;
            float min_identity;
            int _minhash_kmer_size;

            int wfa_mismatch_score;
            int wfa_gap_opening_score;
            int wfa_gap_extension_score;

            int wfa_patching_mismatch_score;
            int wfa_patching_gap_opening_score1;
            int wfa_patching_gap_extension_score1;
            int wfa_patching_gap_opening_score2;
            int wfa_patching_gap_extension_score2;
            
            float mashmap_estimated_identity;
            // WFlign parameters
            int wflign_mismatch_score;
            int wflign_gap_opening_score;
            int wflign_gap_extension_score;
            float wflign_max_mash_dist;
            int wflign_min_wavefront_length;
            int wflign_max_distance_threshold;
            uint64_t wflign_max_len_major;
            uint64_t wflign_max_len_minor;
            int erode_k;
            int64_t chain_gap;
            int min_inversion_length;
            int max_patching_score;
            // Query
            const std::string* query_name;
            char* query;
            uint64_t query_total_length;
            uint64_t query_offset;
            uint64_t query_length;
            bool query_is_rev;
            // Target
            const std::string* target_name;
            char* target;
            uint64_t target_total_length;
            uint64_t target_offset;
            uint64_t target_length;
            // Output
            std::ostream* out;
#ifdef WFA_PNG_TSV_TIMING
            bool emit_tsv;
            std::ostream* out_tsv;
            const std::string* prefix_wavefront_plot_in_png;
            uint64_t wfplot_max_size;
            bool emit_patching_tsv;
            std::ostream* out_patching_tsv;
#endif
            bool force_wflign;
            bool merge_alignments;
            bool emit_md_tag;
            bool paf_format_else_sam;
            bool no_seq_in_sam;
            bool force_biwfa_alignment;
            // Setup
            WFlign(
                    const uint16_t segment_length,
                    const float min_identity,
                    const bool force_wflign,
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
                    const int max_patching_score);
            // Set output configuration
            void set_output(
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
                    const bool no_seq_in_sam);
            // WFling affine
            void wflign_affine_wavefront(
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
                    const uint64_t target_length);
        };

    } /* namespace wavefront */

} /* namespace wflign */

/*
* DTO ()
*/
typedef struct {
    // WFlign
    wflign::wavefront::WFlign* wflign;
    // Parameters
    int pattern_length;
    int text_length;
    uint16_t step_size;
    uint16_t segment_length_to_use;
    int minhash_kmer_size;
    float max_mash_dist_to_evaluate;
    float mash_sketch_rate;
    float inception_score_max_ratio;
    // Alignments and sketches
    robin_hood::unordered_flat_map<uint64_t,alignment_t*>* alignments;
    std::vector<std::vector<rkmh::hash_t>*>* query_sketches;
    std::vector<std::vector<rkmh::hash_t>*>* target_sketches;
    // Subsidiary WFAligner
    wfa::WFAlignerGapAffine* wf_aligner;
//    // Bidirectional
//    wfa::WFAlignerGapAffine* wflambda_aligner;
//    int last_breakpoint_v;
//    int last_breakpoint_h;
//    wflign_penalties_t* wfa_affine_penalties;
    // Stats
#ifdef WFA_PNG_TSV_TIMING
    uint64_t num_alignments;
    uint64_t num_alignments_performed;
#endif
    // For performance improvements
    uint64_t max_num_sketches_in_memory;
    uint64_t num_sketches_allocated;
#ifdef WFA_PNG_TSV_TIMING
    // wfplot
    bool emit_png;
    robin_hood::unordered_set<uint64_t>* high_order_dp_matrix_mismatch;
#endif
} wflign_extend_data_t;

#endif /* WFLIGN_HPP_ */
