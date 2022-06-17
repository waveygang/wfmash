#include <stddef.h>
#include <stdlib.h>
#include <cassert>
#include <chrono>
#include <iterator>
#include <string>

#include "wflign.hpp"
#include "wflign_patch.hpp"

// Namespaces
namespace wflign {
namespace wavefront {

/*
* Configuration
*/
#define MAX_LEN_FOR_PURE_WFA    20000 // only for low-divergence, otherwise disabled
#define MIN_WF_LENGTH           256
//#define MAX_DIST_THRESHOLD      1024

/*
* DTO ()
*/
typedef struct {
    // WFlign
    WFlign* wflign;
    // Parameters
    int pattern_length;
    int text_length;
    uint16_t step_size;
    uint16_t segment_length_to_use;
    int minhash_kmer_size;
    int wfa_min_wavefront_length;
    int wfa_max_distance_threshold;
    float max_mash_dist_to_evaluate;
    float mash_sketch_rate;
    float inception_score_max_ratio;
    // Sequences data
    int v_max;
    int h_max;
    // Alignments and sketches
    robin_hood::unordered_flat_map<uint64_t,alignment_t*>* alignments;
    std::vector<std::vector<rkmh::hash_t>*>* query_sketches;
    std::vector<std::vector<rkmh::hash_t>*>* target_sketches;
    // Subsidiary WFAligner
    wfa::WFAlignerGapAffine* wf_aligner;
    wflign_penalties_t* wfa_affine_penalties;
    // Stats
    uint64_t num_alignments;
    uint64_t num_alignments_performed;
    // wfplot
    bool emit_png;
    robin_hood::unordered_set<uint64_t>* high_order_dp_matrix_mismatch;
} wflign_extend_data_t;

/*
* Utils
*/
inline uint64_t encode_pair(
    int v,
    int h) {
return ((uint64_t)v << 32) | (uint64_t)h;
}

inline void decode_pair(uint64_t pair, int *v, int *h) {
*v = pair >> 32;
*h = pair & 0x00000000FFFFFFFF;
}

/*
* Setup
*/
WFlign::WFlign(
    const uint16_t segment_length,
    const float min_identity,
    const int _minhash_kmer_size,
    const int wfa_mismatch_score,
    const int wfa_gap_opening_score,
    const int wfa_gap_extension_score,
    const int wflambda_min_wavefront_length,
    const int wflambda_max_distance_threshold,
    const float mashmap_estimated_identity,
    const int wflign_mismatch_score,
    const int wflign_gap_opening_score,
    const int wflign_gap_extension_score,
    const float wflign_max_mash_dist,
    const uint64_t wflign_max_len_major,
    const uint64_t wflign_max_len_minor,
    const int erode_k) {
    // Parameters
    this->segment_length = segment_length;
    this->min_identity = min_identity;
    this->_minhash_kmer_size = _minhash_kmer_size;
    this->wfa_mismatch_score = wfa_mismatch_score;
    this->wfa_gap_opening_score = wfa_gap_opening_score;
    this->wfa_gap_extension_score = wfa_gap_extension_score;
    this->wflambda_min_wavefront_length = wflambda_min_wavefront_length;
    this->wflambda_max_distance_threshold = wflambda_max_distance_threshold;
    this->mashmap_estimated_identity = mashmap_estimated_identity;
    this->wflign_mismatch_score = wflign_mismatch_score;
    this->wflign_gap_opening_score = wflign_gap_opening_score;
    this->wflign_gap_extension_score = wflign_gap_extension_score;
    this->wflign_max_mash_dist = wflign_max_mash_dist;
    this->wflign_max_len_major = wflign_max_len_major;
    this->wflign_max_len_minor = wflign_max_len_minor;
    this->erode_k = erode_k;
    // Query
    this->query_name = NULL;
    this->query = NULL;
    this->query_total_length = 0;
    this->query_offset = 0;
    this->query_length = 0;
    this->query_is_rev = false;
    // Target
    this->target_name = NULL;
    this->target = NULL;
    this->target_total_length = 0;
    this->target_offset = 0;
    this->target_length = 0;
    // Output
    this->out = NULL;
    this->emit_tsv = false;
    this->out_tsv = NULL;
    this->prefix_wavefront_plot_in_png = NULL;
    this->wfplot_max_size = 0;
    this->merge_alignments = false;
    this->emit_md_tag = false;
    this->paf_format_else_sam = false;
    this->no_seq_in_sam = false;
}
/*
* Output configuration
*/
void WFlign::set_output(
    std::ostream* const out,
    const bool emit_tsv,
    std::ostream* const out_tsv,
    const std::string &wfplot_filepath,
    const uint64_t wfplot_max_size,
    const bool merge_alignments,
    const bool emit_md_tag,
    const bool paf_format_else_sam,
    const bool no_seq_in_sam) {
    this->out = out;
    this->emit_tsv = emit_tsv;
    this->out_tsv = out_tsv;
    this->prefix_wavefront_plot_in_png = &wfplot_filepath;
    this->wfplot_max_size = wfplot_max_size;
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
    const int query_length = wflign.query_length;
    const int target_length = wflign.target_length;
    const int step_size = extend_data->step_size;
    const int segment_length_to_use = extend_data->segment_length_to_use;
    const int pattern_length = extend_data->pattern_length;
    const int text_length = extend_data->text_length;
    robin_hood::unordered_flat_map<uint64_t,alignment_t*>& alignments = *(extend_data->alignments);
    std::vector<std::vector<rkmh::hash_t>*>& query_sketches = *(extend_data->query_sketches);
    std::vector<std::vector<rkmh::hash_t>*>& target_sketches = *(extend_data->target_sketches);
    // wfplots
    const bool emit_png = extend_data->emit_png;
    robin_hood::unordered_set<uint64_t>& high_order_dp_matrix_mismatch = *(extend_data->high_order_dp_matrix_mismatch);
    // Check match
    bool is_a_match = false;
    if (v >= 0 && h >= 0 && v < pattern_length && h < text_length) {
        const uint64_t k = encode_pair(v, h);
        const auto f = alignments.find(k); // high-level of WF-inception
        if (f != alignments.end()) {
            is_a_match = (alignments[k] != nullptr);
        } else {
            const int query_begin = v * step_size;
            const int target_begin = h * step_size;

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
                            extend_data->minhash_kmer_size,
                            extend_data->wfa_min_wavefront_length,
                            extend_data->wfa_max_distance_threshold,
                            extend_data->max_mash_dist_to_evaluate,
                            extend_data->mash_sketch_rate,
                            extend_data->inception_score_max_ratio,
                            *extend_data->wf_aligner,
                            *extend_data->wfa_affine_penalties,
                            *aln);
            if (wflign.emit_tsv) {
                // 0) Mis-match, alignment skipped
                // 1) Mis-match, alignment performed
                // 2) Match, alignment performed
                *(wflign.out_tsv) << v << "\t" << h << "\t"
                                  << (alignment_performed ? (aln->ok ? 2 : 1) : 0)
                                  << std::endl;
            }
            ++(extend_data->num_alignments);
            if (alignment_performed) {
                ++(extend_data->num_alignments_performed);
                if (aln->ok){
                    is_a_match = true;
                    alignments[k] = aln;
                } else {
                    alignments[k] = nullptr;
                }
            } else {
                if (emit_png) {
                    // Save only the mismatches, as they are not cached
                    high_order_dp_matrix_mismatch.insert(encode_pair(v, h));
                }
            }
            if (!is_a_match) {
                delete aln;
            }

            // cleanup old sketches
            if (v > extend_data->v_max) {
                extend_data->v_max = v;
                if (v >= wflign.wflambda_max_distance_threshold) {
                    auto &s = query_sketches[v - wflign.wflambda_max_distance_threshold];
                    // The C++ language guarantees that delete p will do
                    // nothing if p is equal to NULL
                    delete s;
                    s = nullptr;
                }
            }
            if (h > extend_data->h_max) {
                extend_data->h_max = h;
                if (h >= wflign.wflambda_max_distance_threshold) {
                    auto &s = target_sketches[h - wflign.wflambda_max_distance_threshold];
                    // The C++ language guarantees that delete p will do
                    // nothing if p is equal to NULL
                    delete s;
                    s = nullptr;
                }
            }
        }
    } else if (h < 0 || v < 0) { // It can be removed using an edit-distance
        // mode as high-level of WF-inception
        is_a_match = true;
    }
    return is_a_match;
}
//auto trace_match = [&](const int &v, const int &h) {
//    if (v >= 0 && h >= 0 && v < pattern_length && h < text_length) {
//        const uint64_t k = encode_pair(v, h);
//        const auto f = alignments.find(k);
//        if (f != alignments.end() && alignments[k] != nullptr) {
//            auto *aln = alignments[k];
//            trace.push_back(aln);
//            aln->keep = true;
//            ++num_alignments;
//            return true;
//        }
//        return false;
//    } else {
//        return false;
//    }
//};
int wflambda_trace_match(
    robin_hood::unordered_flat_map<uint64_t,alignment_t*>& alignments,
    wfa::WFAlignerGapAffine& wflambda_aligner,
    std::vector<alignment_t*>& trace,
    const int pattern_length,
    const int text_length) {
// Retrieve CIGAR
char* cigar_ops;
int cigar_length;
wflambda_aligner.getAlignmentCigar(&cigar_ops,&cigar_length);
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

    //auto minhash_kmer_size = _minhash_kmer_size;
    auto minhash_kmer_size = std::max(8, std::min(19, (int)std::round(1.0 / (1.0 - mashmap_estimated_identity))));

    // Set penalties
    wflign_penalties_t wfa_affine_penalties;
    if (wfa_mismatch_score > 0 && wfa_gap_opening_score > 0 && wfa_gap_extension_score > 0){
        wfa_affine_penalties.match = 0;
        wfa_affine_penalties.mismatch = wfa_mismatch_score;
        wfa_affine_penalties.gap_opening = wfa_gap_opening_score;
        wfa_affine_penalties.gap_extension = wfa_gap_extension_score;
    } else {
        wfa_affine_penalties.match = 0;
        wfa_affine_penalties.mismatch = 4;
        wfa_affine_penalties.gap_opening = 6;
        wfa_affine_penalties.gap_extension = 1;
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

    int _erode_k = 0;
    float inception_score_max_ratio = 1.0 + 0.5 / std::pow(mashmap_estimated_identity,3);
    float max_mash_dist_to_evaluate = std::min(0.95, 0.05 / std::pow(mashmap_estimated_identity,9));
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
    } else if (mashmap_estimated_identity >= 0.85) {
    } else if (mashmap_estimated_identity >= 0.8) {
    } else if (mashmap_estimated_identity >= 0.75) {
    } else {
    }

    // heuristic setting of erosion
    _erode_k = std::max(127.0,std::round(1.0/(1.0-mashmap_estimated_identity)));

    // override erosion if not given on input
    if (erode_k < 0) {
        erode_k = _erode_k;
    }

    // override max mash dist if given on input
    if (wflign_max_mash_dist > 0) {
        max_mash_dist_to_evaluate = wflign_max_mash_dist;
    }

    // accumulate runs of matches in reverse order
    // then trim the cigars of successive mappings
    std::vector<alignment_t*> trace;
    const auto start_time = std::chrono::steady_clock::now();

    if (
            (query_length <= segment_length * 8 || target_length <= segment_length * 8) ||
            (mashmap_estimated_identity >= 0.99 // about the limit of what our reduction thresholds allow
             && query_length <= MAX_LEN_FOR_PURE_WFA && target_length <= MAX_LEN_FOR_PURE_WFA)
            ) {
        uint64_t num_alignments = 0;
        uint64_t num_alignments_performed = 0;
        wfa::WFAlignerGapAffine* wf_aligner =
                new wfa::WFAlignerGapAffine(
                        wfa_affine_penalties.mismatch,
                        wfa_affine_penalties.gap_opening,
                        wfa_affine_penalties.gap_extension,
                        wfa::WFAligner::Alignment,
                        wfa::WFAligner::MemoryUltralow);
        //wf_aligner->setHeuristicWFadaptive(MIN_WF_LENGTH,wf_max_dist_threshold);
        const int status = wf_aligner->alignEnd2End(target,target_length,query,query_length);

        auto *aln = new alignment_t();
        aln->j = 0;
        aln->i = 0;

        // aln.mash_dist = mash_dist;
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

        ++num_alignments;
        ++num_alignments_performed;

        trace.push_back(aln);

        const long elapsed_time_wflambda_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - start_time).count();

        // Free old aligner
        delete wf_aligner;

        // use biWFA for all patching
        wf_aligner = new wfa::WFAlignerGapAffine(
            wfa_affine_penalties.mismatch,
            wfa_affine_penalties.gap_opening,
            wfa_affine_penalties.gap_extension,
            wfa::WFAligner::Alignment,
            wfa::WFAligner::MemoryUltralow);
        wf_aligner->setHeuristicNone();

        // write a merged alignment
        write_merged_alignment(
                *out,
                trace,
                *wf_aligner,
                wfa_affine_penalties,
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
                MAX_LEN_FOR_PURE_WFA,
                min_identity,
                elapsed_time_wflambda_ms,
                num_alignments,
                num_alignments_performed,
                mashmap_estimated_identity,
                wflign_max_len_major,
                wflign_max_len_minor,
                erode_k,
                MIN_WF_LENGTH,
                wf_max_dist_threshold,
                prefix_wavefront_plot_in_png,
                wfplot_max_size);

        // Free biWFA aligner
        delete wf_aligner;

    } else {
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

        const uint16_t segment_length_to_use =
                (query_length < segment_length || target_length < segment_length)
                ? std::min(query_length, target_length)
                : segment_length;

        // set up our implicit matrix
        const uint8_t steps_per_segment = 2;
        const uint16_t step_size = segment_length_to_use / steps_per_segment;

        // Pattern & Text
        // If the query_length/target_length are not multiple of step_size, we count
        // a fragment less, and the last one will be longer than segment_length_to_use
        const int pattern_length = (int)query_length / step_size - (query_length % step_size != 0 ? 1 : 0);
        const int text_length = (int)target_length / step_size - (target_length % step_size != 0 ? 1 : 0);

        // uncomment to use reduced WFA locally
        // currently not supported due to issues with traceback when applying
        // WF-reduction on small problems
        const int wfa_min_wavefront_length = 0; // segment_length_to_use / 16;
        const int wfa_max_distance_threshold = 0; // segment_length_to_use / 8;

        wflign_penalties_t wflambda_affine_penalties;
        if (wflign_mismatch_score > 0 && wflign_gap_opening_score > 0 && wflign_gap_extension_score > 0){
            wflambda_affine_penalties.match = 0;
            wflambda_affine_penalties.mismatch = wflign_mismatch_score;
            wflambda_affine_penalties.gap_opening = wflign_gap_opening_score;
            wflambda_affine_penalties.gap_extension = wflign_gap_extension_score;
        } else {
            wflambda_affine_penalties.match = 0;
            wflambda_affine_penalties.mismatch = 4;
            wflambda_affine_penalties.gap_opening = 6;
            wflambda_affine_penalties.gap_extension = 1;
        }

        //std::cerr << "wfa_affine_penalties.mismatch " << wfa_affine_penalties.mismatch << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_opening " << wfa_affine_penalties.gap_opening << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_extension " << wfa_affine_penalties.gap_extension << std::endl;
        //std::cerr << "wflambda_affine_penalties.mismatch " << wflambda_affine_penalties.mismatch << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_opening " << wflambda_affine_penalties.gap_opening << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_extension " << wflambda_affine_penalties.gap_extension << std::endl;
        //std::cerr << "max_mash_dist_to_evaluate " << max_mash_dist_to_evaluate << std::endl;

        // Configure the attributes of the wflambda-aligner
        wfa::WFAlignerGapAffine* wflambda_aligner =
                new wfa::WFAlignerGapAffine(
                        wflambda_affine_penalties.mismatch,
                        wflambda_affine_penalties.gap_opening,
                        wflambda_affine_penalties.gap_extension,
                        wfa::WFAligner::Alignment,
                        wfa::WFAligner::MemoryMed);

        uint64_t _wflambda_max_distance_threshold =
                std::min((uint64_t)std::max(query_length,target_length)/10,
                         (uint64_t)wflambda_max_distance_threshold) / step_size;

        //std::cerr << "wflambda_max_distance_threshold = "
        //          << wflambda_max_distance_threshold * step_size << std::endl;

        if (wflambda_min_wavefront_length || _wflambda_max_distance_threshold) {
            wflambda_aligner->setHeuristicWFadaptive(wflambda_min_wavefront_length,_wflambda_max_distance_threshold);
        } else {
            wflambda_aligner->setHeuristicNone();
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
                        wfa_affine_penalties.gap_opening,
                        wfa_affine_penalties.gap_extension,
                        wfa::WFAligner::Alignment,
                        wfa::WFAligner::MemoryHigh);
        wf_aligner->setHeuristicNone();

        // Save mismatches if wfplots are requsted
        robin_hood::unordered_set<uint64_t> high_order_dp_matrix_mismatch;

        // Setup WFling extend data
        wflign_extend_data_t extend_data;
        extend_data.wflign = this;
        extend_data.pattern_length = pattern_length;
        extend_data.text_length = text_length;
        extend_data.step_size = step_size;
        extend_data.segment_length_to_use = segment_length_to_use;
        extend_data.minhash_kmer_size = minhash_kmer_size;
        extend_data.wfa_min_wavefront_length = wfa_min_wavefront_length;
        extend_data.wfa_max_distance_threshold = wfa_max_distance_threshold;
        extend_data.max_mash_dist_to_evaluate = max_mash_dist_to_evaluate;
        extend_data.mash_sketch_rate = mash_sketch_rate;
        extend_data.inception_score_max_ratio = inception_score_max_ratio;
        extend_data.v_max = 0;
        extend_data.h_max = 0;
        extend_data.alignments = &alignments;
        extend_data.query_sketches = &query_sketches;
        extend_data.target_sketches = &target_sketches;
        extend_data.wf_aligner = wf_aligner;
        extend_data.wfa_affine_penalties = &wfa_affine_penalties;
        extend_data.num_alignments = 0;
        extend_data.num_alignments_performed = 0;
        extend_data.emit_png = !prefix_wavefront_plot_in_png->empty() && wfplot_max_size > 0;
        extend_data.high_order_dp_matrix_mismatch = &high_order_dp_matrix_mismatch;
        wflambda_aligner->setMatchFunct(wflambda_extend_match,(void*)&extend_data);

        // HERE WAS auto extend_match = [&](const int &v, const int &h);
        // HERE WAS auto trace_match = [&](const int &v, const int &h)

        // Align
        wflambda_aligner->alignEnd2EndLambda(pattern_length,text_length);

        // Extract the trace
        extend_data.num_alignments += wflambda_trace_match(
                alignments,*wflambda_aligner,trace,pattern_length,text_length);

        // Free
        delete wflambda_aligner;

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
            {
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
                                             "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".anchors.png";
                encodeOneStep(filename.c_str(), bytes, width, height);
            }

            // Full plot
            {
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
                                             "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".png";
                encodeOneStep(filename.c_str(), bytes, width, height);
            }
        }

        // FIXME: Can we put this inside wflambda_trace_match?
        for (const auto &p : alignments) {
            if (p.second != nullptr && !p.second->keep) {
                delete p.second;
                //p.second = nullptr;
            }
        }

        const long elapsed_time_wflambda_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now()-start_time).count();

        //#define WFLIGN_DEBUG
    #ifdef WFLIGN_DEBUG
        // get alignment score
    const int score = wflambda::edit_cigar_score_gap_affine(
        &affine_wavefronts->edit_cigar, &wflambda_affine_penalties);

    std::cerr << "[wflign::wflign_affine_wavefront] alignment score " << score
    << " for query: " << query_name << " target: " << target_name
    << std::endl;
    #endif

        // clean up sketches
        // The C++ language guarantees that delete p will do nothing if p is equal
        // to NULL
        for (auto &s : query_sketches) {
            delete s;
            s = nullptr;
        }
        for (auto &s : target_sketches) {
            delete s;
            s = nullptr;
        }

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
                int trim_last = 0, trim_curr = 0;
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

                // Free old aligner
                delete wf_aligner;

                // use biWFA for all patching
                wf_aligner = new wfa::WFAlignerGapAffine(
                    wfa_affine_penalties.mismatch,
                    wfa_affine_penalties.gap_opening,
                    wfa_affine_penalties.gap_extension,
                    wfa::WFAligner::Alignment,
                    wfa::WFAligner::MemoryUltralow);
                wf_aligner->setHeuristicNone();

                // write a merged alignment
                write_merged_alignment(
                        *out,
                        trace,
                        *wf_aligner,
                        wfa_affine_penalties,
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
                        segment_length_to_use,
                        min_identity,
                        elapsed_time_wflambda_ms,
                        extend_data.num_alignments,
                        extend_data.num_alignments_performed,
                        mashmap_estimated_identity,
                        wflign_max_len_major,
                        wflign_max_len_minor,
                        erode_k,
                        MIN_WF_LENGTH,
                        wf_max_dist_threshold,
                        prefix_wavefront_plot_in_png,
                        wfplot_max_size);
            } else {
                // todo old implementation (and SAM format is not supported)
                for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
                    write_alignment(
                            *out,
                            **x,
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
                            mashmap_estimated_identity);
                }
            }
        }

        // Free
        delete wf_aligner;
    }
}

} // namespace wavefront
} // namespace wflign
