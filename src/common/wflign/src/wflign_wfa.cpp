#include "wflign_wfa.hpp"
#include <chrono>

// not doing this results in a linker error
#include "WFA/wavefront/wavefront_penalties.c" // TODO fixme

#include "atomic_image.hpp"
#include "lodepng/lodepng.h"

void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height) {
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

namespace wflign {

namespace wavefront {

#define MAX_LEN_FOR_PURE_WFA    20000 // only for low-divergence, otherwise disabled
#define MIN_WF_LENGTH           256
#define MAX_DIST_THRESHOLD      256

wfa::wavefront_aligner_t* get_wavefront_aligner(
    const wfa::affine_penalties_t& wfa_affine_penalties,
    //ToDo: to remove if wavefront_aligner_new will not re-take the seqs' lens in input
    const uint64_t& target_length,
    const uint64_t& query_length,
    const bool& low_memory) {
    // Configure the attributes of the wf-aligner
    wfa::wavefront_aligner_attr_t attributes =
        wfa::wavefront_aligner_attr_default;
    attributes.distance_metric = wfa::gap_affine;
    attributes.affine_penalties = wfa_affine_penalties;
    // attributes.distance_metric = gap_affine2p;
    // attributes.affine2p_penalties = affine2p_penalties;
    attributes.reduction.reduction_strategy =
        wfa::wavefront_reduction_none; // wavefront_reduction_dynamic
    // attributes.reduction.min_wavefront_length = 10;
    // attributes.reduction.max_distance_threshold = 50;
    attributes.alignment_scope =
        wfa::compute_alignment; // alignment_scope_score
    attributes.low_memory = low_memory;
    //wfa::wavefront_aligner_t *const wf_aligner =
    //return wfa::wavefront_aligner_new(target_length, query_length, &attributes);
    return wfa::wavefront_aligner_new(&attributes);
}

void wflign_affine_wavefront(
    std::ostream &out,
    const bool &emit_tsv, std::ostream &out_tsv,
    const std::string &prefix_wavefront_plot_in_png, const uint64_t &wfplot_max_size,
    const bool &merge_alignments,
    const bool &emit_md_tag,
    const bool &paf_format_else_sam, const bool &no_seq_in_sam,
    const std::string &query_name,
    const char *query, const uint64_t &query_total_length,
    const uint64_t &query_offset, const uint64_t &query_length,
    const bool &query_is_rev,
    const std::string &target_name,
    const char *target, const uint64_t &target_total_length,
    const uint64_t &target_offset, const uint64_t &target_length,
    const uint16_t &segment_length,
    const float &min_identity, const int& _minhash_kmer_size,
    const int &wfa_mismatch_score,
    const int &wfa_gap_opening_score,
    const int &wfa_gap_extension_score,
    const int &wflambda_min_wavefront_length, // with these set at 0 we do exact
    // WFA for wflambda
    const int &_wflambda_max_distance_threshold, const float &mashmap_estimated_identity,
    const int &wflign_mismatch_score,
    const int &wflign_gap_opening_score,
    const int &wflign_gap_extension_score,
    const float &wflign_max_mash_dist,
    const uint64_t &wflign_max_len_major, const uint64_t &wflign_max_len_minor,
    const int& _erode_k) {

    if (query_offset + query_length > query_total_length ||
        target_offset + target_length > target_total_length) {
        return;
    }

    auto minhash_kmer_size = _minhash_kmer_size;

    // Set penalties
    wfa::affine_penalties_t wfa_affine_penalties;
    if (wfa_mismatch_score > 0 && wfa_gap_opening_score > 0 && wfa_gap_extension_score > 0){
        wfa_affine_penalties = {
                .match = 0,
                .mismatch = wfa_mismatch_score,
                .gap_opening = wfa_gap_opening_score,
                .gap_extension = wfa_gap_extension_score
        };
        minhash_kmer_size = 17;
    } else {
        if (mashmap_estimated_identity >= 0.99) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 19,
                .gap_opening = 31,
                .gap_extension = 1,
            };
            minhash_kmer_size = 19;
        } else if (mashmap_estimated_identity >= 0.98) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 15,
                .gap_opening = 25,
                .gap_extension = 1,
            };
            minhash_kmer_size = 17;
        } else if (mashmap_estimated_identity >= 0.97) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 13,
                .gap_opening = 21,
                .gap_extension = 1,
            };
            minhash_kmer_size = 17;
        } else if (mashmap_estimated_identity >= 0.95) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 11,
                .gap_opening = 17,
                .gap_extension = 1,
            };
            minhash_kmer_size = 15;
        } else if (mashmap_estimated_identity >= 0.90) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 7,
                .gap_opening = 11,
                .gap_extension = 1,
            };
            minhash_kmer_size = 13;
        } else if (mashmap_estimated_identity >= 0.85) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 6,
                .gap_opening = 9,
                .gap_extension = 1,
            };
            minhash_kmer_size = 11;
        } else if (mashmap_estimated_identity >= 0.80) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 5,
                .gap_opening = 8,
                .gap_extension = 1,
            };
            minhash_kmer_size = 11;
        } else if (mashmap_estimated_identity >= 0.75) {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 4,
                .gap_opening = 6,
                .gap_extension = 1,
            };
            minhash_kmer_size = 11;
        } else {
            wfa_affine_penalties = {
                .match = 0,
                .mismatch = 3,
                .gap_opening = 5,
                .gap_extension = 1,
            };
            minhash_kmer_size = 9;
        }
    }

    // heuristic bound on the max mash dist, adaptive based on estimated
    // identity the goal here is to sparsify the set of alignments in the
    // wflambda layer we then patch up the gaps between them

    int erode_k = 0;
    float inception_score_max_ratio = 3;
    float max_mash_dist_to_evaluate = 1;
    float mash_sketch_rate = 1;

    if (mashmap_estimated_identity >= 0.99) {
        max_mash_dist_to_evaluate = 0.05;
        mash_sketch_rate = 0.125;
        inception_score_max_ratio = 2;
        erode_k = 13;
    } else if (mashmap_estimated_identity >= 0.98) {
        max_mash_dist_to_evaluate = 0.05;
        mash_sketch_rate = 0.125;
        inception_score_max_ratio = 2;
        erode_k = 13;
    } else if (mashmap_estimated_identity >= 0.97) {
        max_mash_dist_to_evaluate = 0.075;
        mash_sketch_rate = 0.125;
        inception_score_max_ratio = 3;
        erode_k = 11;
    } else if (mashmap_estimated_identity >= 0.95) {
        max_mash_dist_to_evaluate = 0.15;
        mash_sketch_rate = 0.25;
        inception_score_max_ratio = 3;
        erode_k = 9;
    } else if (mashmap_estimated_identity >= 0.9) {
        max_mash_dist_to_evaluate = 0.3;
        mash_sketch_rate = 0.3;
        inception_score_max_ratio = 4;
        erode_k = 7;
    } else if (mashmap_estimated_identity >= 0.85) {
        max_mash_dist_to_evaluate = 0.4;
        mash_sketch_rate = 0.35;
        inception_score_max_ratio = 5;
        erode_k = 0;
    } else if (mashmap_estimated_identity >= 0.8) {
        max_mash_dist_to_evaluate = 0.5;
        mash_sketch_rate = 0.4;
        inception_score_max_ratio = 6;
        erode_k = 0;
    } else if (mashmap_estimated_identity >= 0.75) {
        max_mash_dist_to_evaluate = 0.6;
        mash_sketch_rate = 0.45;
        inception_score_max_ratio = 7;
        erode_k = 0;
    } else {
        max_mash_dist_to_evaluate = 0.7;
        mash_sketch_rate = 0.5;
        inception_score_max_ratio = 8;
        erode_k = 0;
    }

    // override erosion if given on input
    if (_erode_k >= 0) {
        erode_k = _erode_k;
    }

    // override max mash dist if given on input
    if (wflign_max_mash_dist > 0) {
        max_mash_dist_to_evaluate = wflign_max_mash_dist;
    }

    // accumulate runs of matches in reverse order
    // then trim the cigars of successive mappings
    std::vector<alignment_t *> trace;

    uint64_t num_alignments = 0;
    uint64_t num_alignments_performed = 0;
    const auto start_time = std::chrono::steady_clock::now();

    // if we expect the alignment to be low divergence, and the mapping is less than 50kb
    // it's faster to just align directly with WFA
    if (mashmap_estimated_identity >= 0.99 // about the limit of what our reduction thresholds allow
        && query_length <= MAX_LEN_FOR_PURE_WFA && target_length <= MAX_LEN_FOR_PURE_WFA) {
        wfa::wavefront_aligner_t* const wf_aligner = get_wavefront_aligner(wfa_affine_penalties,
                                                                           target_length,
                                                                           query_length,
                                                                           true);
        wfa::wavefront_reduction_set_adaptive(&wf_aligner->reduction,
                                              MIN_WF_LENGTH,
                                              MAX_DIST_THRESHOLD);

        auto *aln = new alignment_t();
        wfa::wavefront_aligner_resize(wf_aligner, target_length, query_length);

        const int status = wfa::wavefront_align(wf_aligner,
                                                target, target_length,
                                                query, query_length);

        aln->j = 0;
        aln->i = 0;

        // aln.mash_dist = mash_dist;
        aln->ok = status == WF_ALIGN_SUCCESSFUL;

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

            wflign_edit_cigar_copy(&aln->edit_cigar, &wf_aligner->cigar);

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
                        std::chrono::steady_clock::now() - start_time)
                        .count();

        // write a merged alignment
        write_merged_alignment(
                out, trace, wf_aligner, &wfa_affine_penalties,
                emit_md_tag,
                paf_format_else_sam, no_seq_in_sam,
                query,
                query_name, query_total_length,
                query_offset, query_length, query_is_rev,
                target,
                target_name, target_total_length, target_offset, target_length,
                segment_length,
                MAX_LEN_FOR_PURE_WFA, min_identity,
                elapsed_time_wflambda_ms, num_alignments,
                num_alignments_performed, mashmap_estimated_identity,
                wflign_max_len_major, wflign_max_len_minor,
                erode_k,
                inception_score_max_ratio,
                MIN_WF_LENGTH, MAX_DIST_THRESHOLD,
                prefix_wavefront_plot_in_png, wfplot_max_size);

        // Free
        wfa::wavefront_aligner_delete(wf_aligner);

    } else { // regular wflign
        if (emit_tsv) {
            out_tsv << "# query_name=" << query_name << std::endl;
            out_tsv << "# query_start=" << query_offset << std::endl;
            out_tsv << "# query_end=" << query_offset+query_length << std::endl;
            out_tsv << "# target_name=" << target_name << std::endl;
            out_tsv << "# target_start=" << target_offset << std::endl;
            out_tsv << "# target_end=" << target_offset+target_length << std::endl;
            out_tsv << "# info: 0) mismatch, mash-distance > threshold; 1) mismatch, WFA-score >= max_score; 2) match, WFA-score < max_score" << std::endl;
            out_tsv << "v" << "\t" << "h" << "\t" << "info" << std::endl;
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

        wflambda::affine_penalties_t wflambda_affine_penalties;
        if (wflign_mismatch_score > 0 && wflign_gap_opening_score > 0 && wflign_gap_extension_score > 0){
            wflambda_affine_penalties = {
                    .match = 0,
                    .mismatch = wflign_mismatch_score,
                    .gap_opening = wflign_gap_opening_score,
                    .gap_extension = wflign_gap_extension_score
            };
        } else {
            wflambda_affine_penalties = {
                .match = wfa_affine_penalties.match,
                .mismatch = wfa_affine_penalties.mismatch,
                .gap_opening = wfa_affine_penalties.gap_opening,
                .gap_extension = wfa_affine_penalties.gap_extension
            };
        }

        //std::cerr << "wfa_affine_penalties.mismatch " << wfa_affine_penalties.mismatch << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_opening " << wfa_affine_penalties.gap_opening << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_extension " << wfa_affine_penalties.gap_extension << std::endl;
        //std::cerr << "wflambda_affine_penalties.mismatch " << wflambda_affine_penalties.mismatch << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_opening " << wflambda_affine_penalties.gap_opening << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_extension " << wflambda_affine_penalties.gap_extension << std::endl;
        //std::cerr << "max_mash_dist_to_evaluate " << max_mash_dist_to_evaluate << std::endl;

        // Configure the attributes of the wflambda-aligner
        wflambda::wavefront_aligner_attr_t attributes =
                wflambda::wavefront_aligner_attr_default;
        attributes.distance_metric = wflambda::gap_affine;
        attributes.affine_penalties = wflambda_affine_penalties;
        // attributes.distance_metric = gap_affine2p;
        // attributes.affine2p_penalties = affine2p_penalties;

        uint64_t wflambda_max_distance_threshold =
            std::min((uint64_t)_wflambda_max_distance_threshold,
                     std::max((uint64_t)query_length, (uint64_t)target_length)/10) / step_size;

        if (wflambda_min_wavefront_length || wflambda_max_distance_threshold) {
            attributes.reduction.reduction_strategy =
                    wflambda::wavefront_reduction_dynamic; // wavefront_reduction_dynamic
                    attributes.reduction.min_wavefront_length = wflambda_min_wavefront_length;
                    attributes.reduction.max_distance_threshold = wflambda_max_distance_threshold;
        } else {
            attributes.reduction.reduction_strategy =
                    wflambda::wavefront_reduction_none; // wavefront_reduction_dynamic
        }
        attributes.alignment_scope = wflambda::alignment_scope_alignment; // alignment_scope_score
        attributes.low_memory = true;
        wflambda::wavefront_aligner_t *const wflambda_aligner = wflambda::wavefront_aligner_new(
                pattern_length, text_length, &attributes);


        // save computed alignments in a pair-indexed map
        robin_hood::unordered_flat_map<uint64_t, alignment_t *> alignments;

        // allocate vectors to store our sketches
        std::vector<std::vector<rkmh::hash_t> *> query_sketches(pattern_length,
                                                                nullptr);
        std::vector<std::vector<rkmh::hash_t> *> target_sketches(text_length,
                                                                 nullptr);

        wfa::wavefront_aligner_t* const wf_aligner = get_wavefront_aligner(wfa_affine_penalties,
                                                                           segment_length_to_use,
                                                                           segment_length_to_use,
                                                                           false);

        bool emit_png = !prefix_wavefront_plot_in_png.empty() && wfplot_max_size > 0;
        robin_hood::unordered_set<uint64_t> high_order_dp_matrix_mismatch;

        int v_max = 0;
        int h_max = 0;

        auto extend_match = [&](const int &v, const int &h) {
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
                    const auto segment_length_to_use_q = (uint16_t) (v == pattern_length - 1 ? query_length - query_begin : segment_length_to_use);
                    const auto segment_length_to_use_t = (uint16_t) (h == text_length - 1 ? target_length - target_begin : segment_length_to_use);

                    auto *aln = new alignment_t();
                    const bool alignment_performed = do_wfa_segment_alignment(
                            query_name, query, query_sketches[v], query_length,
                            query_begin, target_name, target, target_sketches[h],
                            target_length, target_begin,
                            segment_length_to_use_q,
                            segment_length_to_use_t,
                            step_size, minhash_kmer_size, wfa_min_wavefront_length,
                            wfa_max_distance_threshold,
                            max_mash_dist_to_evaluate, mash_sketch_rate,
                            inception_score_max_ratio,
                            wf_aligner, &wfa_affine_penalties, *aln);
                    if (emit_tsv) {
                        // 0) Mis-match, alignment skipped
                        // 1) Mis-match, alignment performed
                        // 2) Match, alignment performed
                        out_tsv << v << "\t" << h << "\t" << (alignment_performed ? (aln->ok ? 2 : 1) : 0) << std::endl;
                    }
                    //std::cerr << v << "\t" << h << "\t" << (alignment_performed ? (aln->ok ? 2 : 1) : 0) << std::endl;

                    ++num_alignments;
                    if (alignment_performed) {
                        ++num_alignments_performed;
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
                    if (v > v_max) {
                        v_max = v;
                        if (v >= wflambda_max_distance_threshold) {
                            auto &s =
                                    query_sketches[v - wflambda_max_distance_threshold];
                            // The C++ language guarantees that delete p will do
                            // nothing if p is equal to NULL
                            delete s;
                            s = nullptr;
                        }
                    }
                    if (h > h_max) {
                        h_max = h;
                        if (h >= wflambda_max_distance_threshold) {
                            auto &s =
                                    target_sketches[h -
                                    wflambda_max_distance_threshold];
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
        };

        auto trace_match = [&](const int &v, const int &h) {
            if (v >= 0 && h >= 0 && v < pattern_length && h < text_length) {
                const uint64_t k = encode_pair(v, h);
                const auto f = alignments.find(k);
                if (f != alignments.end() && alignments[k] != nullptr) {
                    auto *aln = alignments[k];
                    trace.push_back(aln);
                    aln->keep = true;
                    ++num_alignments;
                    return true;
                }
                return false;
            } else {
                return false;
            }
        };

        // Align
        wflambda::wavefront_aligner_clear__resize(wflambda_aligner, pattern_length,
                                                  text_length);
        wflambda::wavefront_align(wflambda_aligner, extend_match,
                                  trace_match, pattern_length, text_length);
        wflambda::wavefront_aligner_delete(wflambda_aligner);

        if (emit_png) {
            const int wfplot_vmin = 0, wfplot_vmax = pattern_length; //v_max;
            const int wfplot_hmin = 0, wfplot_hmax = text_length; //h_max

            int v_max = wfplot_vmax - wfplot_vmin;
            int h_max = wfplot_hmax - wfplot_hmin;

            const algorithms::color_t COLOR_MASH_MISMATCH = { 0xffefefef };
            const algorithms::color_t COLOR_WFA_MISMATCH = { 0xffff0000 };
            const algorithms::color_t COLOR_WFA_MATCH = { 0xff00ff00 };

            const double scale = std::min(1.0, (double)wfplot_max_size / (double)std::max(v_max, h_max));

            const uint64_t width = (int)(scale * (double)v_max);
            const uint64_t height = (int)(scale * (double)h_max);
            const double source_width = (double)width;
            const double source_height = (double)height;

            const double x_off = 0, y_off = 0;
            const double line_width = 1.0;
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

                        if (v >= wfplot_vmin && v <= wfplot_vmax && h >= wfplot_hmin && h <= wfplot_hmax) {
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
                const std::string filename = prefix_wavefront_plot_in_png +
                        query_name + "_" + std::to_string(query_offset) + "_" + std::to_string(query_offset+query_length) + " _ " + (query_is_rev ? "-" : "+") +
                        "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".1.anchors.png";
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

                    if (v >= wfplot_vmin && v <= wfplot_vmax && h >= wfplot_hmin && h <= wfplot_hmax) {
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

                    if (v >= wfplot_vmin && v <= wfplot_vmax && h >= wfplot_hmin && h <= wfplot_hmax){
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
                const std::string filename = prefix_wavefront_plot_in_png +
                        query_name + "_" + std::to_string(query_offset) + "_" + std::to_string(query_offset+query_length) + " _ " + (query_is_rev ? "-" : "+") +
                        "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".0.full.png";
                encodeOneStep(filename.c_str(), bytes, width, height);
            }
        }

        for (const auto &p : alignments) {
            if (p.second != nullptr && !p.second->keep) {
                delete p.second;
                //p.second = nullptr;
            }
        }

        const long elapsed_time_wflambda_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - start_time)
                        .count();

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
            //#define VALIDATE_WFA_WFLIGN
            #ifdef VALIDATE_WFA_WFLIGN
            if (!trace.front()->validate(query, target)) {
                std::cerr << "first traceback is wrong" << std::endl;
                trace.front()->display();
                assert(false);
            }
            #endif

            if(attributes.low_memory) {
                std::reverse(trace.begin(), trace.end());
            }

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
            trace_pos_t last_pos = {last.j, last.i, &last.edit_cigar,
                                    last.edit_cigar.begin_offset};
                trace_pos_t curr_pos = {curr.j, curr.i, &curr.edit_cigar,
                                        curr.edit_cigar.begin_offset};

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
                // write a merged alignment
                write_merged_alignment(
                        out, trace, wf_aligner, &wfa_affine_penalties,
                        emit_md_tag,
                        paf_format_else_sam, no_seq_in_sam,
                        query,
                        query_name,
                        query_total_length,
                        query_offset, query_length, query_is_rev, target, target_name,
                        target_total_length, target_offset, target_length,
                        segment_length_to_use,
                        MAX_LEN_FOR_PURE_WFA,
                        min_identity,
                        elapsed_time_wflambda_ms, num_alignments,
                        num_alignments_performed, mashmap_estimated_identity,
                        wflign_max_len_major, wflign_max_len_minor,
                        erode_k,
                        inception_score_max_ratio,
                        MIN_WF_LENGTH, MAX_DIST_THRESHOLD,
                        prefix_wavefront_plot_in_png, wfplot_max_size);
            } else {
                // todo old implementation (and SAM format is not supported)
                for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
                    write_alignment(out, **x, query_name, query_total_length,
                                    query_offset, query_length, query_is_rev,
                                    target_name, target_total_length, target_offset,
                                    target_length, min_identity, mashmap_estimated_identity);
                }
            }
        }

        // Free
        wfa::wavefront_aligner_delete(wf_aligner);
    }
}

// accumulate alignment objects
// run the traceback determine which are part of the main chain
// order them and write them out
// needed--- 0-cost deduplication of alignment regions (how????)
//     --- trim the alignment back to the first 1/2 of the query
bool do_wfa_segment_alignment(
    const std::string &query_name, const char *query,
    std::vector<rkmh::hash_t> *&query_sketch, const uint64_t &query_length,
    const uint64_t &j, const std::string &target_name, const char *target,
    std::vector<rkmh::hash_t> *&target_sketch, const uint64_t &target_length,
    const uint64_t &i,
    const uint16_t &segment_length_q,
    const uint16_t &segment_length_t,
    const uint16_t &step_size,
    const uint64_t &minhash_kmer_size,
    const int &min_wavefront_length,
    const int &max_distance_threshold,
    const float &max_mash_dist,
    const float &mash_sketch_rate,
    const float &inception_score_max_ratio,
    wfa::wavefront_aligner_t *const wf_aligner,
    wfa::affine_penalties_t *const affine_penalties,
    alignment_t &aln) {

    // first make the sketches if we haven't yet
    if (query_sketch == nullptr) {
        query_sketch = new std::vector<rkmh::hash_t>();
        *query_sketch = rkmh::hash_sequence(
            query + j, segment_length_q, minhash_kmer_size, segment_length_q * mash_sketch_rate);
    }
    if (target_sketch == nullptr) {
        target_sketch = new std::vector<rkmh::hash_t>();
        *target_sketch = rkmh::hash_sequence(
            target + i, segment_length_t, minhash_kmer_size, segment_length_t * mash_sketch_rate);
    }

    // first check if our mash dist is inbounds
    const float mash_dist =
        rkmh::compare(*query_sketch, *target_sketch, minhash_kmer_size);
    //std::cerr << "mash_dist is " << mash_dist << std::endl;

    // this threshold is set low enough that we tend to randomly sample wflambda
    // matrix cells for alignment the threshold is adaptive, based on the mash
    // distance of the mapping we are aligning we should obtain enough
    // alignments that we can still patch between them
    if (mash_dist > max_mash_dist) {
        // if it isn't, return false
        return false;
    } else {
        // if it is, we'll align

        const int max_score = std::max(segment_length_q, segment_length_t) * inception_score_max_ratio;

        wfa::wavefront_aligner_resize(wf_aligner, segment_length_t,
                                             segment_length_q);

        wfa::wavefront_aligner_set_max_alignment_score(wf_aligner, max_score);
        const int status =
            wfa::wavefront_align(wf_aligner, target + i, segment_length_t,
                                         query + j, segment_length_q);

        aln.j = j;
        aln.i = i;

        // aln.mash_dist = mash_dist;
        aln.ok = status == WF_ALIGN_SUCCESSFUL && wf_aligner->cigar.score < max_score;

        // fill the alignment info if we aligned
        if (aln.ok) {
            aln.query_length = segment_length_q;
            aln.target_length = segment_length_t;
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

            wflign_edit_cigar_copy(&aln.edit_cigar, &wf_aligner->cigar);

#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(aln.edit_cigar, query, target, segment_length_q,
                                segment_length_t, aln.j, aln.i)) {
                std::cerr << "cigar failure after cigar copy in alignment "
                          << aln.j << " " << aln.i << std::endl;
                assert(false);
            }
#endif
        }

        return true;
    }
}

/*
// No more necessary
bool hack_cigar(wfa::cigar_t &cigar, const char *query, const char *target,
                const uint64_t &query_aln_len, const uint64_t &target_aln_len,
                uint64_t j, uint64_t i) {
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    bool ok = true;
    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        if (j >= j_max && i >= i_max) {
            cigar.end_offset = c;
            ok = false;
            break;
        }
        // if new sequence of same moves started
        switch (cigar.operations[c]) {
        case 'M':
            // check that we match
            if (j < j_max && i < i_max && query[j] != target[i]) {
                // std::cerr << "mismatch @ " << j << " " << i << " " <<
                // query[j] << " " << target[i] << std::endl;
                cigar.operations[c] = 'X';
                ok = false;
            }
            if (j >= j_max) {
                // std::cerr << "query out of bounds @ " << j << " " << i << " "
                // << query[j] << " " << target[i] << std::endl;
                cigar.operations[c] = 'D';
                ok = false;
            }
            if (i >= i_max) {
                // std::cerr << "target out of bounds @ " << j << " " << i << "
                // " << query[j] << " " << target[i] << std::endl;
                cigar.operations[c] = 'I';
                ok = false;
            }
            ++j;
            ++i;
            break;
        case 'X':
            if (j < j_max && i < i_max && query[j] == target[i]) {
                cigar.operations[c] = 'M';
                ok = false;
            }
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
    }
    return ok;
}*/

void do_wfa_patch_alignment(const char *query, const uint64_t &j,
                            const uint64_t &query_length, const char *target,
                            const uint64_t &i, const uint64_t &target_length,
                            const int &segment_length,
                            const int &min_wavefront_length,
                            const int &max_distance_threshold,
                            const float& inception_score_max_ratio,
                            wfa::wavefront_aligner_t *const _wf_aligner,
                            wfa::affine_penalties_t *const affine_penalties,
                            alignment_t &aln) {
    const long max_seg_len = 3 * segment_length;
    const bool big_wave = (query_length > max_seg_len || target_length > max_seg_len);
    wfa::wavefront_aligner_t* const wf_aligner
        = big_wave ?
            get_wavefront_aligner(*affine_penalties,
                                  target_length,
                                  query_length,
                                  true)
        : _wf_aligner;

    /*
     std::cerr << "do_wfa_patch q:" << j << " qlen:" << query_length
               << " t:" << i << " tlen:" << target_length << std::endl;

     {
         //std::hash
         std::string query_seq(query+j, query_length);
         //query + j, query_length
         std::string target_seq(query+j, query_length);
         //target + i, target_length
         std::stringstream namess;
         auto target_hash = std::hash<std::string>{}(target_seq);
         auto query_hash = std::hash<std::string>{}(query_seq);
         namess << "wfpatch_" << target_length << "x" << query_length << "_" << target_hash << "_" << query_hash << ".fa";
         std::ofstream out(namess.str());
         out << ">" << target_hash << std::endl << target_seq << std::endl
             << ">" << query_hash << std::endl << query_seq << std::endl;
     }
    */

    // Reduction strategy
    if (query_length < max_distance_threshold &&
        target_length < max_distance_threshold) {
        // wavefront_reduction_none
        wfa::wavefront_reduction_set_none(&wf_aligner->reduction);
    } else {
        wfa::wavefront_reduction_set_adaptive(&wf_aligner->reduction,
                                             min_wavefront_length,
                                             max_distance_threshold);
    }

    const int max_score = affine_penalties->gap_opening +
        (std::min(target_length, query_length)
         * affine_penalties->gap_extension * inception_score_max_ratio * 4);

    wfa::wavefront_aligner_resize(wf_aligner, target_length,
                                         query_length);

    wfa::wavefront_aligner_set_max_alignment_score(wf_aligner, max_score);
    const int status =
        wfa::wavefront_align(wf_aligner, target + i, target_length,
                                     query + j, query_length);

    aln.ok = status == WF_ALIGN_SUCCESSFUL && wf_aligner->cigar.score < max_score;
    if (aln.ok) {
        // No more necessary: correct X/M errors in the cigar
        //hack_cigar(wf_aligner->cigar, query, target, query_length, target_length, j, i);

#ifdef VALIDATE_WFA_WFLIGN
        if (!validate_cigar(wf_aligner->cigar, query, target, query_length,
                            target_length, j, i)) {
            std::cerr << "cigar failure at alignment " << aln.j << " " << aln.i
                      << std::endl;
            unpack_display_cigar(wf_aligner->cigar, query, target, query_length,
                                 target_length, aln.j, aln.i);
            std::cerr << ">query" << std::endl
                      << std::string(query + j, query_length) << std::endl;
            std::cerr << ">target" << std::endl
                      << std::string(target + i, target_length) << std::endl;
            assert(false);
        }
#endif

        wflign_edit_cigar_copy(&aln.edit_cigar, &wf_aligner->cigar);
    }

    if (big_wave) {
        wfa::wavefront_aligner_delete(wf_aligner);
    }
}

bool validate_cigar(const wfa::cigar_t &cigar, const char *query,
                    const char *target, const uint64_t &query_aln_len,
                    const uint64_t &target_aln_len, uint64_t j, uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    bool ok = true;
    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (cigar.operations[c]) {
        case 'M':
            // check that we match
            if (query[j] != target[i]) {
                std::cerr << "mismatch @ " << j << " " << i << " " << query[j]
                          << " " << target[i] << std::endl;
                ok = false;
            }
            if (j >= j_max) {
                std::cerr << "query out of bounds @ " << j << " " << i << " "
                          << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            if (i >= i_max) {
                std::cerr << "target out of bounds @ " << j << " " << i << " "
                          << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            ++j;
            ++i;
            break;
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
    }
    return ok;
}

bool validate_trace(const std::vector<char> &tracev, const char *query,
                    const char *target, const uint64_t &query_aln_len,
                    const uint64_t &target_aln_len, uint64_t j, uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = 0;
    const int end_idx = tracev.size();
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    bool ok = true;
     //std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
     //std::cerr << "j_max " << j_max << " - i_max " << i_max << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (tracev[c]) {
            case 'M':
                if (j < j_max && i < i_max) {
                    // check that we match
                    if (query[j] != target[i]) {
                        std::cerr << "mismatch @ " << tracev[c] << " " << j << " " << i
                                  << " " << query[j] << " " << target[i] << std::endl;
                        ok = false;
                    }
                } else  {
                    if (j >= j_max) {
                        std::cerr << "M - query out of bounds @ " << j << " " << i << std::endl;//" "
                        //<< query[j] << " " << target[i] << std::endl;
                        ok = false;
                    }

                    if (i >= i_max) {
                        std::cerr << "M - target out of bounds @ " << j << " " << i << std::endl;//" "
                        //<< query[j] << " " << target[i] << std::endl;
                        ok = false;
                    }
                }

                ++j;
                ++i;
                break;
            case 'X':
                if (j < j_max && i < i_max) {
                    // check that we don't match
                    if (query[j] == target[i]) {
                        std::cerr << "match @ " << tracev[c] << " " << j << " " << i
                                  << " " << query[j] << " " << target[i] << std::endl;
                        ok = false;
                    }
                } else {
                    if (j >= j_max) {
                        std::cerr << "X - query out of bounds @ " << j << " " << i << std::endl;//" "
                        //<< query[j] << " " << target[i] << std::endl;
                        ok = false;
                    }

                    if (i >= i_max) {
                        std::cerr << "X - target out of bounds @ " << j << " " << i << std::endl;//" "
                        //<< query[j] << " " << target[i] << std::endl;
                        ok = false;
                    }
                }

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
    }
    return ok;
}

bool unpack_display_cigar(const wfa::cigar_t &cigar, const char *query,
                          const char *target, const uint64_t &query_aln_len,
                          const uint64_t &target_aln_len, uint64_t j,
                          uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (cigar.operations[c]) {
        case 'M':
            // check that we match
            std::cerr << "M"
                      << " " << j << " " << i << " "
                      << "q:" << query[j] << " "
                      << "t:" << target[i] << " "
                      << (query[j] == target[i] ? " " : "👾")
                      << (j >= j_max ? "⚡" : " ") << (i >= i_max ? "🔥" : " ")
                      << std::endl;
            ++j;
            ++i;
            break;
        case 'X':
            std::cerr << "X"
                      << " " << j << " " << i << " "
                      << "q:" << query[j] << " "
                      << "t:" << target[i] << " "
                      << (query[j] != target[i] ? " " : "👾")
                      << (j >= j_max ? "⚡" : " ") << (i >= i_max ? "🔥" : " ")
                      << std::endl;
            ++j;
            ++i;
            break;
        case 'I':
            std::cerr << "I"
                      << " " << j << " " << i << " "
                      << "q:" << query[j] << " "
                      << "t:"
                      << "|"
                      << "  " << (j >= j_max ? "⚡" : " ")
                      << (i >= i_max ? "🔥" : " ") << std::endl;
            ++j;
            break;
        case 'D':
            std::cerr << "D"
                      << " " << j << " " << i << " "
                      << "q:"
                      << "|"
                      << " "
                      << "t:" << target[i] << "  " << (j >= j_max ? "⚡" : " ")
                      << (i >= i_max ? "🔥" : " ") << std::endl;
            ++i;
            break;
        default:
            break;
        }
    }
    return true;
}

void write_merged_alignment(
    std::ostream &out, const std::vector<alignment_t *> &trace,
    wfa::wavefront_aligner_t *const wf_aligner,
    wfa::affine_penalties_t *const affine_penalties,
    const bool &emit_md_tag,
    const bool &paf_format_else_sam, const bool &no_seq_in_sam,
    const char *query,
    const std::string &query_name, const uint64_t &query_total_length,
    const uint64_t &query_offset, const uint64_t &query_length,
    const bool &query_is_rev,
    const char *target,
    const std::string &target_name, const uint64_t &target_total_length,
    const uint64_t &target_offset, const uint64_t &target_length,
    const uint16_t &segment_length,
    const uint64_t &max_pure_wfa,
    const float &min_identity, const long &elapsed_time_wflambda_ms,
    const uint64_t &num_alignments, const uint64_t &num_alignments_performed,
    const float &mashmap_estimated_identity,
    const uint64_t &wflign_max_len_major, const uint64_t &wflign_max_len_minor,
    const int &erode_k,
    const float &inception_score_max_ratio,
    const int &min_wf_length, const int &max_dist_threshold,
    const std::string &prefix_wavefront_plot_in_png, const uint64_t &wfplot_max_size,
    const bool &with_endline) {

    int64_t target_pointer_shift = 0;

    uint64_t target_length_mut = target_length;

    // patching parameters
    // we will nibble patching back to this length
    const uint64_t min_wfa_patch_length = 0; //128;

    // we need to get the start position in the query and target
    // then run through the whole alignment building up the cigar
    // finally emitting it
    // our final cigar
    //
    // std::string cigarstr;
    uint64_t matches = 0;
    uint64_t mismatches = 0;
    uint64_t insertions = 0;
    uint64_t inserted_bp = 0;
    uint64_t deletions = 0;
    uint64_t deleted_bp = 0;
    uint64_t query_start = 0;
    uint64_t target_start = 0;
    uint64_t total_query_aligned_length = 0;
    uint64_t total_target_aligned_length = 0;
    uint64_t query_end = 0;
    uint64_t target_end = 0;
    uint64_t total_score = 0;

    // double mash_dist_sum = 0;
    uint64_t ok_alns = 0;

    auto start_time = std::chrono::steady_clock::now();

#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] processing traceback"
              << std::endl;
#endif
    // write trace into single cigar vector
    std::vector<char> tracev;
    {
        // patch: walk the cigar, patching directly when we have simultaneous
        // gaps in query and ref and adding our results to the final trace as we
        // go

#define MAX_NUM_INDELS_TO_LOOK_AT 3
        auto distance_close_big_enough_indels =
            [](const uint32_t indel_len, auto iterator,
               const std::vector<char> &trace) {
                const uint32_t min_indel_len_to_find = indel_len / 3;
                const uint16_t max_dist_to_look_at =
                    std::min(indel_len * 64, (uint32_t)4096);

                // std::cerr << "min_indel_len_to_find " <<
                // min_indel_len_to_find << std::endl; std::cerr <<
                // "max_dist_to_look_at " << max_dist_to_look_at << std::endl;

                auto q = iterator;

                uint8_t num_indels_to_find = MAX_NUM_INDELS_TO_LOOK_AT;
                uint32_t curr_size_close_indel = 0;
                int32_t dist_close_indels = 0;

                //std::cerr << "indel_len " << indel_len << std::endl;
                //std::cerr << "min_indel_len_to_find " << min_indel_len_to_find << std::endl;
                //std::cerr << "max_dist_to_look_at " << max_dist_to_look_at << std::endl;

                while (q != trace.end() &&
                       dist_close_indels < max_dist_to_look_at) {
                    curr_size_close_indel = 0;
                    while (q != trace.end() && (*q == 'I' || *q == 'D')) {
                        ++curr_size_close_indel;

                        ++dist_close_indels;
                        ++q;
                    }
                    // std::cerr << "\t\tcurr_size_close_indel " <<
                    // curr_size_close_indel << std::endl;
                    if (curr_size_close_indel >= min_indel_len_to_find) {
                        // std::cerr << "\t\tnum_indels_to_find " <<
                        // (uint16_t)num_indels_to_find << std::endl;
                        if (--num_indels_to_find == 0) {
                            break;
                        }
                    }

                    while (q != trace.end() &&
                           (dist_close_indels < max_dist_to_look_at) &&
                           *q != 'I' && *q != 'D') {
                        ++dist_close_indels;
                        ++q;
                    }
                }

                //std::cerr << "dist_close_indels " << dist_close_indels << std::endl;
                //std::cerr << "num_indels_found " << MAX_NUM_INDELS_TO_LOOK_AT - num_indels_to_find << std::endl;

                return num_indels_to_find < MAX_NUM_INDELS_TO_LOOK_AT
                           ? dist_close_indels
                           : -1;
            };

        auto patching = [&query, &query_name, &query_length, &query_start,
                         &query_offset, &target, &target_name,
                         &target_length_mut, &target_start, &target_offset,
                         &target_total_length, &target_end,
                         &target_pointer_shift,
                         &segment_length,
                         &wflign_max_len_major,
                         &wflign_max_len_minor,
                         &distance_close_big_enough_indels, &min_wf_length,
                         &inception_score_max_ratio,
                         &max_dist_threshold, &wf_aligner,
                         &affine_penalties](std::vector<char> &unpatched,
                                            std::vector<char> &patched) {
            auto q = unpatched.begin();

            uint64_t query_pos = query_start;
            uint64_t target_pos = target_start;

            uint64_t query_delta = 0;
            uint64_t target_delta = 0;

            bool got_alignment = false;

            // Head patching
            {
                // how long a gap?
                while (q != unpatched.end() && *q == 'I') {
                    ++query_delta;
                    ++q;
                }
                while (q != unpatched.end() && *q == 'D') {
                    ++target_delta;
                    ++q;
                }

                if (query_delta > 0 && query_delta < wflign_max_len_minor) {
                    // Semi-global mode for patching the heads

                    // TODO: when we will have semi-global WFA
                    // nibble forward if we're below the correct length
                    // ...
                    // TODO: when we will have semi-global WFA

//                    std::cerr << "HEAD patching in"
//                              << query_name << " "
//                              << query_offset << "@ " << query_pos <<
//                              " - " << query_delta
//                              << " --- "
//                              << target_name << " " << target_offset
//                              << " @ " <<
//                              target_pos << " - "
//                              << target_delta
//                              << std::endl;

                    // min_wfa_patch_length-bps of margins to manage insertions in the query
                    const uint64_t delta_to_ask = query_delta + min_wfa_patch_length <= target_delta ? 0 : query_delta + min_wfa_patch_length - target_delta;

                    uint64_t target_delta_to_shift = 0;
                    uint64_t target_pos_x, target_start_x;
                    int64_t target_pointer_shift_x;

                    // Note that target_pos >= target_start
                    if (target_pos >= delta_to_ask) {
                        // std::cerr << "A\n";
                        // Easy, we don't have to manage 'negative' indexes for
                        // the target array
                        target_delta_to_shift = delta_to_ask;

                        target_pos_x = target_pos - target_delta_to_shift;
                        target_start_x = target_pos_x < target_start
                                             ? target_pos_x
                                             : target_start;

                        target_pointer_shift_x = target_pointer_shift; // + 0;
                    } else {
                        // std::cerr << "B\n";
                        target_pos_x = 0;
                        target_start_x = 0;

                        const uint64_t positions_to_get =
                            delta_to_ask - target_pos;

                        // Manage negative indexes
                        if (target_offset - target_pointer_shift >= positions_to_get) {
                            // std::cerr << "B.1\n";
                            target_pointer_shift_x = target_pointer_shift +
                                                     (int64_t)positions_to_get;

                            target_delta_to_shift =
                                delta_to_ask; // we can get all the positions we
                                              // need
                        } else {
                            // std::cerr << "B.2\n";
                            // We can't get all the positions we need

                            target_pointer_shift_x = (int64_t)
                                target_offset; //<=> target_pointer_shift +
                                               //(int64_t)target_offset -
                                               //target_pointer_shift;

                            target_delta_to_shift = target_offset + target_pos -
                                                    target_pointer_shift;
                        }
                    }

                    //                            std::cerr << "target_start "
                    //                            << target_start << std::endl;
                    //                            std::cerr << "target_start_x "
                    //                            << target_start_x <<
                    //                            std::endl; std::cerr <<
                    //                            "target_pos " << target_pos <<
                    //                            std::endl; std::cerr <<
                    //                            "target_pos_x " <<
                    //                            target_pos_x << std::endl;
                    //                            std::cerr << "target_offset "
                    //                            << target_offset << std::endl;
                    //                            std::cerr <<
                    //                            "target_pointer_shift_x " <<
                    //                            target_pointer_shift_x <<
                    //                            std::endl; std::cerr <<
                    //                            "target_delta_to_shift " <<
                    //                            target_delta_to_shift <<
                    //                            std::endl;

                    const uint64_t target_delta_x =
                        target_delta + target_delta_to_shift;

                    if (target_delta_x > 0) {
//                        std::cerr << "B HEAD patching in "
//                                  << query_name << " "
//                                  << query_offset << "@ " << query_pos <<
//                                  " - " << query_delta
//                                << " --- "
//                                << target_name << " " << target_offset
//                                << " @ " <<
//                                target_pos << " - "
//                                << target_delta_x
//                                << std::endl;

                        std::string query_rev(query + query_pos, query_delta);
                        std::reverse(query_rev.begin(), query_rev.end());

                        std::string target_rev(target - target_pointer_shift_x +
                                                   target_pos_x,
                                               target_delta_x);
                        std::reverse(target_rev.begin(), target_rev.end());

//                        std::cerr << "query: ";
//                        for (int i = 0; i < query_delta; ++i) {
//                            std::cerr << query_rev[i];
//                        }
//                        std::cerr << "\ntarget: ";;
//                        for (int i = 0; i < target_delta_x; ++i) {
//                            std::cerr << target_rev[i];
//                        }
//                        std::cerr << std::endl;

                        wfa::wavefront_aligner_t* const wf_aligner_heads
                                = get_wavefront_aligner(*affine_penalties,
                                                        target_rev.size()+1,
                                                        query_rev.size()+1, true);
                        wavefront_aligner_set_alignment_free_ends(
                                wf_aligner_heads,
                                0,
                                0,
                                0,
                                query_rev.size());
                        wfa::wavefront_reduction_set_adaptive(&wf_aligner_heads->reduction,
                                                              min_wf_length,
                                                              max_dist_threshold);
                        const int status =
                                wfa::wavefront_align(wf_aligner_heads, target_rev.c_str(), target_rev.size(),
                                                     query_rev.c_str(), query_rev.size());

                        if (status == WF_ALIGN_SUCCESSFUL) {
                            //hack_cigar(wf_aligner_heads->cigar, query_rev.c_str(), target_rev.c_str(), query_rev.size(), target_rev.size(), 0, 0);

#ifdef VALIDATE_WFA_WFLIGN
                            if (!validate_cigar(wf_aligner_heads->cigar, query_rev.c_str(), target_rev.c_str(), query_rev.size(),
                                                target_rev.size(), 0, 0)) {
                                std::cerr << "cigar failure at head alignment " << std::endl;
                                unpack_display_cigar(wf_aligner_heads->cigar, query_rev.c_str(), target_rev.c_str(), query_rev.size(),
                                                     target_rev.size(), 0, 0);
                                std::cerr << ">query" << std::endl
                                          << std::string(query_rev.c_str(), query_rev.size()) << std::endl;
                                std::cerr << ">target" << std::endl
                                          << std::string(target_rev.c_str(), target_rev.size()) << std::endl;
                                assert(false);
                            }
#endif

                            //std::cerr << "Head patching\n";
                            got_alignment = true;

                            target_pos = target_pos_x;
                            target_delta = target_delta_x;
                            target_pointer_shift = target_pointer_shift_x;

                            target_start = target_start_x;
                            target_length_mut += target_delta_to_shift;

                            for(int xxx = wf_aligner_heads->cigar.end_offset - 1; xxx >= wf_aligner_heads->cigar.begin_offset; --xxx) {
                                //std::cerr << wf_aligner_heads->cigar.operations[xxx];
                                patched.push_back(wf_aligner_heads->cigar.operations[xxx]);
                            }
                            //std::cerr << "\n";
                        }
                        wfa::wavefront_aligner_delete(wf_aligner_heads);
                    }
                }

                // add in stuff if we didn't align
                if (!got_alignment) {
                    for (uint64_t i = 0; i < query_delta; ++i) {
                        patched.push_back('I');
                    }
                    for (uint64_t i = 0; i < target_delta; ++i) {
                        patched.push_back('D');
                    }
                }

                query_pos += query_delta;
                target_pos += target_delta;

                query_delta = 0;
                target_delta = 0;
            }

            // Patching in the middle
            // get to the first match ... we'll not yet try to patch the
            // alignment tips
            while (q != unpatched.end()) {
                while (q != unpatched.end() && (*q == 'M' || *q == 'X')) {
                    /*
                    std::cerr << "q: " << query[query_pos] << " "
                    << "t: " << target[target_pos - target_pointer_shift] <<
                    std::endl;
                    */
                    if (query_pos >= query_length ||
                        target_pos >= target_length_mut) {
                        std::cerr << "[wflign::wflign_affine_wavefront] "
                                     "corrupted traceback (out of bounds) for "
                                  << query_name << " " << query_offset << " "
                                  << target_name << " " << target_offset
                                  << std::endl;
                        exit(1);
                    }

                    if (*q == 'M') {
                        if (query[query_pos] !=
                            target[target_pos - target_pointer_shift]) {
                            std::cerr << "[wflign::wflign_affine_wavefront] "
                                         "corrupted traceback (M, but there is "
                                         "a mismatch) for "
                                      << query_name << " " << query_offset
                                      << " " << target_name << " "
                                      << target_offset << std::endl;
                            exit(1);
                        }
                    } else {
                        if (query[query_pos] ==
                            target[target_pos - target_pointer_shift]) {
                            std::cerr << "[wflign::wflign_affine_wavefront] "
                                         "corrupted traceback (X, but there is "
                                         "a match) for "
                                      << query_name << " " << query_offset
                                      << " " << target_name << " "
                                      << target_offset << std::endl;
                            exit(1);
                        }
                    }

                    patched.push_back(*q);
                    ++query_pos;
                    ++target_pos;
                    ++q;
                }

                // how long a gap?
                while (q != unpatched.end() && *q == 'I') {
                    ++query_delta;
                    ++q;
                }
                while (q != unpatched.end() && *q == 'D') {
                    ++target_delta;
                    ++q;
                }

                // how long was our last gap?
                // if it's long enough, patch it
                int32_t size_region_to_repatch = 0;

                do {
                    got_alignment = false;

                    if (size_region_to_repatch > 0 ||
                        (query_delta > 0 && target_delta > 0) ||
                        (query_delta > 2 || target_delta > 2) &&
                        (query_delta < wflign_max_len_major &&
                         target_delta < wflign_max_len_major) &&
                        (query_delta < wflign_max_len_minor ||
                         target_delta < wflign_max_len_minor)) {

                        { //if (false) {

                            int32_t distance_close_indels = -1; /* (query_delta > 3 || target_delta > 3) ?
                                distance_close_big_enough_indels(std::max(query_delta, target_delta), q, unpatched) :
                                -1;*/
                            // std::cerr << "distance_close_indels " <<
                            // distance_close_indels << std::endl;
                            // Trigger the patching if there is a dropout
                            // (consecutive Is and Ds) or if there is a close and
                            // big enough indel forward
                            if (size_region_to_repatch > 0 ||
                                (query_delta > 0 && target_delta > 0) ||
                                (distance_close_indels > 0)) {
#ifdef WFLIGN_DEBUG
                                // std::cerr << "query_delta " << query_delta <<
                                // "\n"; std::cerr << "target_delta " << target_delta
                                // << "\n"; std::cerr << "distance_close_indel " <<
                                // distance_close_indel << "\n";

                                std::cerr << "[wflign::wflign_affine_wavefront] "
                                    "patching in "
                                          << query_name << " " << query_offset
                                          << " @ " << query_pos << " - "
                                          << query_delta << " " << target_name
                                          << " " << target_offset << " @ "
                                          << target_pos << " - " << target_delta
                                          << std::endl;
#endif
                                /*
                                  std::cerr << "A patching in "
                                  << query_name << " " << query_offset << " @ " << query_pos << " - " << query_delta << " "
                                  << target_name << " " << target_offset << " @ " << target_pos << " - " << target_delta
                                  << std::endl;
                                */

                                // if we are continuing a patch, we can't nibble
                                // backward too much to avoid the risk of going in
                                // endless loop
                                if (size_region_to_repatch > 0) {
                                    // nibble backward
                                    while (!patched.empty() &&
                                           size_region_to_repatch > 0) {
                                        const auto &c = patched.back();
                                        switch (c) {
                                        case 'M':
                                        case 'X':
                                            --query_pos;
                                            --target_pos;
                                            ++query_delta;
                                            ++target_delta;
                                            break;
                                        case 'I':
                                            ++query_delta;
                                            --query_pos;
                                            break;
                                        case 'D':
                                            ++target_delta;
                                            --target_pos;
                                            break;
                                        default:
                                            break;
                                        }
                                        patched.pop_back();
                                        --size_region_to_repatch;
                                    }

                                    //distance_close_indels = distance_close_big_enough_indels(std::max(query_delta, target_delta), q, unpatched);
                                } else {
                                    // nibble backward if we're below the correct
                                    // length
                                    while (
                                        !patched.empty() &&
                                        (query_delta < (min_wfa_patch_length / 2) ||
                                         target_delta <
                                         (min_wfa_patch_length / 2))) {
                                        const auto &c = patched.back();
                                        switch (c) {
                                        case 'M':
                                        case 'X':
                                            --query_pos;
                                            --target_pos;
                                            ++query_delta;
                                            ++target_delta;
                                            break;
                                        case 'I':
                                            ++query_delta;
                                            --query_pos;
                                            break;
                                        case 'D':
                                            ++target_delta;
                                            --target_pos;
                                            break;
                                        default:
                                            break;
                                        }
                                        patched.pop_back();
                                    }
                                }

                                // nibble forward if we're below the correct length
                                while (q != unpatched.end() &&
                                       (query_delta < min_wfa_patch_length ||
                                        target_delta < min_wfa_patch_length)) {
                                    const auto &c = *q++;
                                    switch (c) {
                                    case 'M':
                                    case 'X':
                                        ++query_delta;
                                        ++target_delta;
                                        break;
                                    case 'I':
                                        ++query_delta;
                                        break;
                                    case 'D':
                                        ++target_delta;
                                        break;
                                    default:
                                        break;
                                    }

                                    --distance_close_indels;
                                }

                                /*
                                  std::cerr
                                  << "B patching in "
                                  << query_name << " " << query_offset << " @ " << query_pos << " - " << query_delta << " "
                                  << target_name << " " << target_offset << " @ " << target_pos << " - " << target_delta
                                  << std::endl;
                                */

                                // Nibble until the close, big enough indel is
                                // reached Important when the patching can't be
                                // computed correctly without including the next
                                // indel
                                while (q != unpatched.end() &&
                                       distance_close_indels > 0) {
                                    const auto &c = *q++;
                                    switch (c) {
                                    case 'M':
                                    case 'X':
                                        ++query_delta;
                                        ++target_delta;
                                        break;
                                    case 'I':
                                        ++query_delta;
                                        break;
                                    case 'D':
                                        ++target_delta;
                                        break;
                                    default:
                                        break;
                                    }

                                    --distance_close_indels;
                                }


                                /*
                                  std::cerr << "C patching in "
                                  << query_name << " " << query_offset << " @ " << query_pos << " - " << query_delta << " "
                                  << target_name << " " << target_offset << " @ " << target_pos << " - " << target_delta
                                  << std::endl;
                                */

                                // check forward if there are other Is/Ds to merge
                                // in the current patch
                                while (q != unpatched.end() &&
                                       (*q == 'I' || *q == 'D') &&
                                       ((query_delta < wflign_max_len_major &&
                                         target_delta < wflign_max_len_major) &&
                                        (query_delta < wflign_max_len_minor ||
                                         target_delta < wflign_max_len_minor))) {
                                    const auto &c = *q++;
                                    if (c == 'I') {
                                        ++query_delta;
                                    } else {
                                        ++target_delta;
                                    }
                                }

                                /*
                                  std::cerr << "D patching in "
                                  << query_name << " " << query_offset << " @ " << query_pos << " - " << query_delta << " "
                                  << target_name << " " << target_offset << " @ " << target_pos << " - " << target_delta
                                  << std::endl;
                                */

                                // check backward if there are other Is/Ds to merge
                                // in the current patch it will eventually nibble
                                // the Is/Ds left from the last patch
                                while (!patched.empty() &&
                                       (patched.back() == 'I' ||
                                        patched.back() == 'D') &&
                                       ((query_delta < wflign_max_len_major &&
                                         target_delta < wflign_max_len_major) &&
                                        (query_delta < wflign_max_len_minor ||
                                         target_delta < wflign_max_len_minor))) {
                                    const auto &c = patched.back();
                                    if (c == 'I') {
                                        ++query_delta;
                                        --query_pos;
                                    } else {
                                        ++target_delta;
                                        --target_pos;
                                    }
                                    patched.pop_back();
                                }

                            }
                            /*
                              std::cerr << "E patching in "
                              << query_name << " " << query_offset << " @ " << query_pos << " - " << query_delta << " "
                              << target_name << " " << target_offset << " @ " << target_pos << " - " << target_delta
                              << std::endl;
                            */

                            size_region_to_repatch = 0;
                            // we need to be sure that our nibble made the
                            // problem long enough For affine WFA to be correct
                            // (to avoid trace-back errors), it must be at least
                            // 10 nt
                            { //if (query_delta >= 10 && target_delta >= 10) {
                                alignment_t patch_aln;
                                // WFA is only global
                                do_wfa_patch_alignment(
                                    query, query_pos, query_delta,
                                    target - target_pointer_shift, target_pos,
                                    target_delta, segment_length,
                                    min_wf_length, max_dist_threshold,
                                    inception_score_max_ratio,
                                    wf_aligner, affine_penalties, patch_aln);
                                if (patch_aln.ok) {
                                    // std::cerr << "got an ok patch aln" <<
                                    // std::endl;
                                    got_alignment = true;
                                    const int start_idx =
                                        patch_aln.edit_cigar.begin_offset;
                                    const int end_idx =
                                        patch_aln.edit_cigar.end_offset;
                                    for (int i = start_idx; i < end_idx; i++) {
                                        // std::cerr <<
                                        // patch_aln.edit_cigar.operations[i];
                                        patched.push_back(
                                            patch_aln.edit_cigar.operations[i]);
                                    }
                                    // std::cerr << "\n";

                                    // Check if there are too many indels in the
                                    // patch
                                    uint32_t size_indel = 0;
                                    for (int i = end_idx - 1; i >= start_idx;
                                         --i) {
                                        // std::cerr <<
                                        // patch_aln.edit_cigar.operations[i];
                                        if (patch_aln.edit_cigar
                                            .operations[i] == 'I' ||
                                            patch_aln.edit_cigar
                                            .operations[i] == 'D') {
                                            ++size_indel;
                                            ++size_region_to_repatch;
                                        } else {
                                            // Not too big, to avoid repatching
                                            // structural variants boundaries
                                            if (size_indel > 7 &&
                                                size_indel <= 4096 &&
                                                size_indel <
                                                (end_idx - start_idx)) {
                                                break;
                                            }

                                            ++size_region_to_repatch;
                                            size_indel = 0;
                                        }
                                    }
                                    // std::cerr << std::endl;

                                    // Not too big, to avoid repatching
                                    // structural variants boundaries
                                    //std::cerr << "size_region_to_repatch " << size_region_to_repatch << std::endl;
                                    //std::cerr << "end_idx - start_idx " << end_idx - start_idx << std::endl;
                                    if (size_indel > 7 && size_indel <= 4096 &&
                                        size_region_to_repatch <
                                        (end_idx - start_idx)) {
                                        //std::cerr << "REPATCH " << std::endl;
                                    } else {
                                        size_region_to_repatch = 0;
                                    }
                                }
                            }
                        } // if false --- to disable patching
                    }

                    // add in stuff if we didn't align
                    if (!got_alignment) {
                        for (uint64_t i = 0; i < query_delta; ++i) {
                            patched.push_back('I');
                        }
                        for (uint64_t i = 0; i < target_delta; ++i) {
                            patched.push_back('D');
                        }
                    }

                    // std::cerr << "query_delta " << query_delta << std::endl;
                    // std::cerr << "target_delta " << target_delta <<
                    // std::endl;
                    query_pos += query_delta;
                    target_pos += target_delta;

                    query_delta = 0;
                    target_delta = 0;
                } while (size_region_to_repatch > 0);
            }

            // Tail patching
            {
                // TODO: when we will have semi-global WFA
                // nibble backward if we're below the correct length
                /*bool nibble_fwd = true;
                while (!patched.empty() && query_delta < min_wfa_patch_length) {
                    const auto& c = patched.back();
                    switch (c) {
                        case 'M': case 'X':
                            --query_pos; --target_pos;
                            ++query_delta; ++target_delta; break;
                        case 'I': ++query_delta; --query_pos; break;
                        case 'D': ++target_delta; --target_pos; break;
                        default: break;
                    }
                    patched.pop_back();
                }
                */
                // TODO: when we will have semi-global WFA

                // Important: the last patch (in the middle of the traceback)
                // can generate a tail check backward if there are other Is/Ds
                // to merge in the current patch
                while (!patched.empty() &&
                       (patched.back() == 'I' || patched.back() == 'D') &&
                       ((query_delta < wflign_max_len_major &&
                         target_delta < wflign_max_len_major) &&
                        (query_delta < wflign_max_len_minor ||
                         target_delta < wflign_max_len_minor))) {
                    const auto &c = patched.back();
                    if (c == 'I') {
                        ++query_delta;
                        --query_pos;
                    } else {
                        ++target_delta;
                        --target_pos;
                    }
                    patched.pop_back();
                }

                got_alignment = false;

                if (query_delta > 0 && query_delta < wflign_max_len_minor) {
//                                            std::cerr << "A TAIL patching in "
//                                                      << query_name << " " <<
//                                                      query_offset << " @ " <<
//                                                      query_pos << " - " <<
//                                                      query_delta << " "
//                                                      << target_name << " " <<
//                                                      target_offset << " @ "
//                                                      << target_pos - target_pointer_shift << " - "
//                                                      << target_delta
//                                                      << std::endl;

                    // min_wfa_patch_length-bps of margins to manage insertions in the query
                    const uint64_t delta_to_ask = query_delta + min_wfa_patch_length <= target_delta ? 0 : query_delta + min_wfa_patch_length - target_delta;

                    // there is a piece of query
                    auto target_delta_x =
                        target_delta +
                        ((target_offset - target_pointer_shift) + target_pos +
                                     target_delta + delta_to_ask <
                                 target_total_length
                             ? delta_to_ask
                             : target_total_length -
                                   ((target_offset - target_pointer_shift) +
                                    target_pos + target_delta));

                    if (target_delta_x > 0) {
//                        std::cerr << "B TAIL patching in "
//                                  << query_name << " " <<
//                                  query_offset << " @ " <<
//                                  query_pos << " - " <<
//                                  query_delta << " "
//                                  << target_name << " " <<
//                                  target_offset << " @ "
//                                  << target_pos - target_pointer_shift << " - "
//                                  << target_delta_x
//                                  << std::endl;

                        wfa::wavefront_aligner_t* const wf_aligner_tails
                                = get_wavefront_aligner(*affine_penalties,
                                                        target_delta_x+1,
                                                        query_delta+1, true);
                        wavefront_aligner_set_alignment_free_ends(
                                wf_aligner_tails,
                                0,
                                0,
                                0,
                                query_delta);
                        wfa::wavefront_reduction_set_adaptive(&wf_aligner_tails->reduction,
                                                              min_wf_length,
                                                              max_dist_threshold);
                        const int status =
                                wfa::wavefront_align(wf_aligner_tails, target - target_pointer_shift + target_pos, target_delta_x,
                                                     query + query_pos, query_delta);

                        if (status == WF_ALIGN_SUCCESSFUL) {
                            //hack_cigar(wf_aligner_tails->cigar, query + query_pos, target - target_pointer_shift + target_pos, query_delta, target_delta_x, 0, 0);

#ifdef VALIDATE_WFA_WFLIGN
                            if (!validate_cigar(wf_aligner_tails->cigar, query + query_pos, target - target_pointer_shift + target_pos, query_delta,
                                                target_delta_x, 0, 0)) {
                                std::cerr << "cigar failure at head alignment " << std::endl;
                                unpack_display_cigar(wf_aligner_tails->cigar, query + query_pos, target - target_pointer_shift + target_pos, query_delta,
                                                     target_delta_x, 0, 0);
                                std::cerr << ">query" << std::endl
                                          << std::string(query + query_pos, query_delta) << std::endl;
                                std::cerr << ">target" << std::endl
                                          << std::string(target - target_pointer_shift + target_pos, target_delta_x) << std::endl;
                                assert(false);
                            }
#endif

                            //std::cerr << "Tail patching\n";
                            got_alignment = true;

                            {
                                // if (target_pos + target_delta_x >
                                // target_length_mut) {
                                //    target_end += (target_pos + target_delta_x
                                //    - target_length_mut); target_length_mut =
                                //    target_pos + target_delta_x;
                                //}
                                const uint32_t inc =
                                    target_delta_x - target_delta;
                                target_end += inc;
                                target_length_mut += inc;
                            }

                            target_delta = target_delta_x;

                            for(int xxx = wf_aligner_tails->cigar.begin_offset; xxx < wf_aligner_tails->cigar.end_offset; ++xxx) {
                                //std::cerr << wf_aligner_tails->cigar.operations[xxx];
                                patched.push_back(wf_aligner_tails->cigar.operations[xxx]);
                            }
                            //std::cerr << "\n";
                        }
                        wfa::wavefront_aligner_delete(wf_aligner_tails);
                    }
                }

                if (!got_alignment) {
                    // add in our tail gap / softclip
                    for (uint64_t i = 0; i < query_delta; ++i) {
                        patched.push_back('I');
                    }
                    for (uint64_t i = 0; i < target_delta; ++i) {
                        patched.push_back('D');
                    }
                }
                // query_pos += query_delta; // not used
                // target_pos += target_delta;
            }

#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] got unsorted "
                         "patched traceback: ";
            for (auto c : patched) {
                std::cerr << c;
            }
            std::cerr << std::endl;
#endif
        };

        std::vector<char> pre_tracev;
        {
            std::vector<char> erodev;
            {
                std::vector<char> rawv;

                // copy
#ifdef WFLIGN_DEBUG
                std::cerr
                    << "[wflign::wflign_affine_wavefront] copying traceback"
                    << std::endl;
#endif
                for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
                    auto &aln = **x;
                    if (aln.ok) {
                        if (ok_alns == 0) {
                            query_start = aln.j;
                            target_start = aln.i;
                        }
                        ++ok_alns;
                        if (query_end && aln.j > query_end) {
                            const int len = aln.j - query_end;
                            for (uint64_t i = 0; i < len; ++i) {
                                rawv.push_back('I');
                            }
                        }
                        if (target_end && aln.i > target_end) {
                            const int len = aln.i - target_end;
                            for (uint64_t i = 0; i < len; ++i) {
                                rawv.push_back('D');
                            }
                        }
                        uint64_t target_aligned_length = 0;
                        uint64_t query_aligned_length = 0;
                        const int start_idx = aln.edit_cigar.begin_offset;
                        const int end_idx = aln.edit_cigar.end_offset;
                        for (int i = start_idx; i < end_idx; i++) {
                            const auto &c = aln.edit_cigar.operations[i];
                            switch (c) {
                            case 'M':
                            case 'X':
                                ++query_aligned_length;
                                ++target_aligned_length;
                                break;
                            case 'I':
                                ++query_aligned_length;
                                break;
                            case 'D':
                                ++target_aligned_length;
                                break;
                            default:
                                break;
                            }
                            rawv.push_back(c);
                        }
                        query_end = aln.j + query_aligned_length;
                        target_end = aln.i + target_aligned_length;
                    }
                    // clean up
                    delete *x;
                }

#ifdef VALIDATE_WFA_WFLIGN
                if (!validate_trace(rawv, query, target - target_pointer_shift,
                                    query_length, target_length_mut,
                                    query_start, target_start)) {
                    std::cerr
                        << "cigar failure in rawv (at end) "
                        << "\t" << query_name << "\t" << query_total_length
                        << "\t"
                        << query_offset + (query_is_rev
                                               ? query_length - query_end
                                               : query_start)
                        << "\t"
                        << query_offset + (query_is_rev
                                               ? query_length - query_start
                                               : query_end)
                        << "\t" << (query_is_rev ? "-" : "+") << "\t"
                        << target_name << "\t" << target_total_length << "\t"
                        << target_offset - target_pointer_shift + target_start
                        << "\t" << target_offset + target_end << std::endl;
                    exit(1);
                }
#endif

#ifdef WFLIGN_DEBUG
                std::cerr << "[wflign::wflign_affine_wavefront] eroding "
                             "traceback at k="
                          << erode_k << std::endl;
#endif

                // erode by removing matches < k
                for (uint64_t i = 0; i < rawv.size();) {
                    if (rawv[i] == 'M' || rawv[i] == 'X') {
                        uint64_t j = i;
                        while (++j < rawv.size() &&
                               (rawv[j] == 'M' || rawv[j] == 'X')) {
                        }
                        if (j - i < erode_k) {
                            while (i < j) {
                                erodev.push_back('D');
                                erodev.push_back('I');
                                ++i;
                            }
                        } else {
                            while (i < j) {
                                erodev.push_back(rawv[i++]);
                            }
                        }
                    } else {
                        erodev.push_back(rawv[i++]);
                    }
                }
            }

#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_trace(erodev, query, target - target_pointer_shift,
                                query_length, target_length_mut, query_start,
                                target_start)) {
                std::cerr << "cigar failure in erodev "
                          << "\t" << query_name << "\t" << query_total_length
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_end
                                                 : query_start)
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_start
                                                 : query_end)
                          << "\t" << (query_is_rev ? "-" : "+") << "\t"
                          << target_name << "\t" << target_total_length << "\t"
                          << target_offset - target_pointer_shift + target_start
                          << "\t" << target_offset + target_end << std::endl;
                exit(1);
            }
#endif

#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] normalizing eroded "
                         "traceback"
                      << std::endl;
#endif

            // normalize: sort so that I<D and otherwise leave it as-is
            sort_indels(erodev);

#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] got normalized "
                         "eroded traceback: ";
            for (auto c : erodev) {
                std::cerr << c;
            }
            std::cerr << std::endl;
#endif

            // std::cerr << "FIRST PATCH ROUND
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            if (erode_k > 0) {
                patching(erodev, pre_tracev);
            } else {
                pre_tracev = erodev;
            }

#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_trace(tracev, query,
                                target - target_pointer_shift, query_length,
                                target_length_mut, query_start, target_start)) {
                std::cerr << "cigar failure in pre_tracev "
                          << "\t" << query_name << "\t" << query_total_length
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_end
                                                 : query_start)
                          << "\t"
                          << query_offset + (query_is_rev
                                                 ? query_length - query_start
                                                 : query_end)
                          << "\t" << (query_is_rev ? "-" : "+") << "\t"
                          << target_name << "\t" << target_total_length << "\t"
                          << target_offset - target_pointer_shift + target_start
                          << "\t" << target_offset + target_end << std::endl;
                exit(1);
            }
#endif

            // normalize: sort so that I<D and otherwise leave it as-is
            sort_indels(tracev);
        }

        // std::cerr << "SECOND PATCH ROUND
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        if (erode_k > 0) {
            patching(pre_tracev, tracev);
        } else {
            tracev = pre_tracev;
        }
    }

    // normalize the indels
    //sort_indels(tracev);

#ifdef WFLIGN_DEBUG
    std::cerr
        << "[wflign::wflign_affine_wavefront] got full patched traceback: ";
    for (auto c : tracev) {
        std::cerr << c;
    }
    std::cerr << std::endl;
#endif

#ifdef VALIDATE_WFA_WFLIGN
    //            std::cerr << "query_length: " << query_length << std::endl;
    //            std::cerr << "target_length_mut: " << target_length_mut <<
    //            std::endl; std::cerr << "query_start: " << query_start <<
    //            std::endl; std::cerr << "target_start: " << target_start <<
    //            std::endl;

    if (!validate_trace(tracev, query, target - target_pointer_shift,
                        query_length, target_length_mut, query_start,
                        target_start)) {
        std::cerr
            << "cigar failure at alignment (before head/tail del trimming) "
            << "\t" << query_name << "\t" << query_total_length << "\t"
            << query_offset +
                   (query_is_rev ? query_length - query_end : query_start)
            << "\t"
            << query_offset +
                   (query_is_rev ? query_length - query_start : query_end)
            << "\t" << (query_is_rev ? "-" : "+") << "\t" << target_name << "\t"
            << target_total_length << "\t"
            << target_offset - target_pointer_shift + target_start << "\t"
            << target_offset + target_end << std::endl;
        exit(1);
    }
#endif

    // trim deletions at start and end of tracev
    uint64_t begin_offset = 0;
    uint64_t end_offset = 0;
    {
        uint64_t trim_del_first = 0;
        uint64_t trim_del_last = 0;

        // 1.) sort initial ins/del to put del < ins
        auto first_non_indel = tracev.begin();
        while (first_non_indel != tracev.end() &&
               (*first_non_indel == 'D' || *first_non_indel == 'I')) {
            ++first_non_indel;
        }
        std::sort(tracev.begin(), first_non_indel,
                  [](char a, char b) { return a < b; });
        // 2.) find first non-D in tracev --> tracev_begin
        //   a.) add to target_start this count
        auto first_non_del = tracev.begin();
        while (first_non_del != tracev.end() && *first_non_del == 'D') {
            ++first_non_del;
        }
        trim_del_first = std::distance(tracev.begin(), first_non_del);
        target_start += trim_del_first;
        // target_length_mut -= trim_del_first;

        // 3.) count D's at end of tracev --> tracev_end
        //   b.) subtract from target_end this count
        auto last_non_del = tracev.rbegin();
        while (last_non_del != tracev.rend() && *last_non_del == 'D') {
            ++last_non_del;
        }
        trim_del_last = std::distance(tracev.rbegin(), last_non_del);
        target_end -= trim_del_last;
        // target_length_mut -= trim_del_last;

        begin_offset = trim_del_first;
        end_offset = tracev.size() - trim_del_last;
    }

    /*
#ifdef VALIDATE_WFA_WFLIGN
    if (!validate_trace(tracev, query, target - target_pointer_shift,
query_length, target_length_mut, query_start, target_start)) { std::cerr <<
"cigar failure at alignment (after head/tail del trimming) "
                  << "\t" << query_name
                  << "\t" << query_total_length
                  << "\t" << query_offset + (query_is_rev ? query_length -
query_end : query_start)
                  << "\t" << query_offset + (query_is_rev ? query_length -
query_start : query_end)
                  << "\t" << (query_is_rev ? "-" : "+")
                  << "\t" << target_name
                  << "\t" << target_total_length
                  << "\t" << target_offset - target_pointer_shift + target_start
                  << "\t" << target_offset + target_end << std::endl;
        exit(1);
    }
#endif
    */

    bool emit_png = !prefix_wavefront_plot_in_png.empty() && wfplot_max_size > 0;
    if (emit_png) {

        const int step_size = (segment_length / 2);

        //const int pattern_length = (int)query_length;
        //const int text_length = (int)target_length;
        const int pattern_length = (int)query_length / step_size - (query_length % step_size != 0 ? 1 : 0);
        const int text_length = (int)target_length / step_size - (target_length % step_size != 0 ? 1 : 0);

        const int wfplot_vmin = 0, wfplot_vmax = pattern_length;
        const int wfplot_hmin = 0, wfplot_hmax = text_length;

        int v_max = wfplot_vmax - wfplot_vmin;
        int h_max = wfplot_hmax - wfplot_hmin;

        const algorithms::color_t COLOR_MASH_MISMATCH = { 0xffefefef };
        const algorithms::color_t COLOR_WFA_MISMATCH = { 0xffff0000 };
        const algorithms::color_t COLOR_WFA_MATCH = { 0xff00ff00 };

        const double scale = std::min(1.0, (double)wfplot_max_size / (double)std::max(v_max, h_max));

        const uint64_t width = (int)(scale * (double)v_max);
        const uint64_t height = (int)(scale * (double)h_max);
        const double source_width = (double)width;
        const double source_height = (double)height;

        const double x_off = 0, y_off = 0;
        const double line_width = 1.0;
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

        // Plot the traceback
        {
            algorithms::atomic_image_buf_t image(width, height,
                                                 source_width, source_height,
                                                 source_min_x, source_min_y);

            
            uint64_t v = query_start; // position in the pattern
            uint64_t h = target_start; // position in the text
            int64_t last_v = -1;
            int64_t last_h = -1;
            for (const auto& c : tracev) {
                switch (c) {
                case 'M':
                case 'X':
                    ++v;
                    ++h;
                    {
                        uint64_t _v = (v / step_size);
                        uint64_t _h = (h / step_size);
                        if ((_v != last_v && _h != last_h)
                            && _v >= wfplot_vmin && _v <= wfplot_vmax
                            && _h >= wfplot_hmin && _h <= wfplot_hmax) {
                            algorithms::xy_d_t xy0 = {
                                (_v * scale) - x_off,
                                (_h * scale) + y_off
                            };
                            xy0.into(source_min_x, source_min_y,
                                     source_width, source_height,
                                     0, 0,
                                     width, height);
                            plot_point(xy0, image, COLOR_WFA_MATCH);
                            last_v = _v;
                            last_h = _h;
                        }
                    }
                    break;
                case 'I':
                    ++v;
                    break;
                case 'D':
                    ++h;
                    break;
                default:
                    break;
                }
                //std::cerr << "plot cell " << v << "," << h << std::endl;
            }

            auto bytes = image.to_bytes();
            const std::string filename = prefix_wavefront_plot_in_png +
                query_name + "_" + std::to_string(query_offset) + "_" + std::to_string(query_offset+query_length) + " _ " + (query_is_rev ? "-" : "+") +
                "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".2.trace.png";
            encodeOneStep(filename.c_str(), bytes, width, height);
        }
    }

    // convert trace to cigar, get correct start and end coordinates
    char *cigarv = alignment_to_cigar(
        tracev, begin_offset, end_offset,
        total_target_aligned_length, total_query_aligned_length, matches,
        mismatches, insertions, inserted_bp, deletions, deleted_bp);

    const double gap_compressed_identity =
        (double)matches /
        (double)(matches + mismatches + insertions + deletions);

    if (gap_compressed_identity >= min_identity) {
        const uint64_t edit_distance = mismatches + inserted_bp + deleted_bp;

        // identity over the full block
        const double block_identity =
                (double)matches / (double)(matches + edit_distance);

        auto write_tag_and_md_string = [&](std::ostream &out, const char *c,
                const int target_start) {
            out << "MD:Z:";

            char last_op = '\0';
            int last_len = 0;
            int t_off = target_start, l_MD = 0;
            int l = 0;
            int x = 0;
            while (c[x] != '\0') {
                while (isdigit(c[x]))
                    ++x;
                char op = c[x];
                int len = 0;
                std::from_chars(c + l, c + x, len);
                l = ++x;
                if (last_len) {
                    if (last_op == op) {
                        len += last_len;
                    } else {
                        // std::cerr << t_off << "   " << last_len << last_op <<
                        // std::endl;

                        if (last_op == '=') {
                            l_MD += last_len;
                            t_off += last_len;
                        } else if (last_op == 'X') {
                            for (uint64_t ii = 0; ii < last_len; ++ii) {
                                out << l_MD
                                << target[t_off + ii - target_pointer_shift];
                                l_MD = 0;
                            }

                            t_off += last_len;
                        } else if (last_op == 'D') {
                            out << l_MD << "^";
                            for (uint64_t ii = 0; ii < last_len; ++ii) {
                                out << target[t_off + ii - target_pointer_shift];
                            }

                            l_MD = 0;
                            t_off += last_len;
                        }
                    }
                }
                last_op = op;
                last_len = len;
            }

            if (last_len) {
                // std::cerr << t_off << "   " << last_len << last_op << std::endl;

                if (last_op == '=') {
                    out << last_len + l_MD;
                } else if (last_op == 'X') {
                    for (uint64_t ii = 0; ii < last_len; ++ii) {
                        out << l_MD << target[t_off + ii - target_pointer_shift];
                        l_MD = 0;
                    }
                    out << "0";
                } else if (last_op == 'I') {
                    out << l_MD;
                } else if (last_op == 'D') {
                    out << l_MD << "^";
                    for (uint64_t ii = 0; ii < last_len; ++ii) {
                        out << target[t_off + ii - target_pointer_shift];
                    }
                    out << "0";
                }
            }
        };

        const long elapsed_time_patching_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_time)
                .count();

        const std::string timings_and_num_alignements =
            "wt:i:" + std::to_string(elapsed_time_wflambda_ms) +
            "\tpt:i:" + std::to_string(elapsed_time_patching_ms) +
            "\taa:i:" + std::to_string(num_alignments) +
            "\tap:i:" + std::to_string(num_alignments_performed);

        if (paf_format_else_sam) {
            out << query_name << "\t" << query_total_length << "\t"
                << query_offset +
                       (query_is_rev ? query_length - query_end : query_start)
                << "\t"
                << query_offset +
                       (query_is_rev ? query_length - query_start : query_end)
                << "\t" << (query_is_rev ? "-" : "+") << "\t" << target_name
                << "\t" << target_total_length << "\t"
                << target_offset - target_pointer_shift + target_start << "\t"
                << target_offset + target_end << "\t" << matches << "\t"
                << std::max(total_target_aligned_length,
                            total_query_aligned_length)
                << "\t"
                << std::round(float2phred(1.0 - block_identity))
                //<< "\t" << "as:i:" << total_score
                << "\t"
                << "gi:f:" << gap_compressed_identity * 100 << "\t"
                << "bi:f:"
                << block_identity * 100
                //<< "\t" << "md:f:" << mash_dist_sum / trace.size()
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "ii:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "dd:i:" << deleted_bp
                << "\t"
                << "md:f:" << mashmap_estimated_identity * 100;

            if (emit_md_tag) {
                out << "\t";

                write_tag_and_md_string(out, cigarv, target_start);
            }

            out << "\t" << timings_and_num_alignements << "\t"
                << "cg:Z:" << cigarv << "\n";
        } else {
            out << query_name                          // Query template NAME
                << "\t" << (query_is_rev ? "16" : "0") // bitwise FLAG
                << "\t" << target_name // Reference sequence NAME
                << "\t"
                << target_offset - target_pointer_shift + target_start +
                       1 // 1-based leftmost mapping POSition
                << "\t"
                << std::round(
                       float2phred(1.0 - block_identity)) // MAPping Quality
                << "\t";

            // CIGAR
            const uint64_t query_start_pos =
                    query_offset +
                    (query_is_rev ? query_length - query_end : query_start);
            const uint64_t query_end_pos =
                    query_offset +
                    (query_is_rev ? query_length - query_start : query_end);

            if (query_is_rev) {
                if (query_length > query_end_pos) {
                    out << (query_length - query_end_pos) << "H";
                }
            } else {
                if (query_start_pos > 0) {
                    out << query_start_pos << "H";
                }
            }
            out << cigarv;
            if (query_is_rev) {
                if (query_start_pos > 0) {
                    out << query_start_pos << "H";
                }
            } else {
                if (query_length > query_end_pos) {
                    out << (query_length - query_end_pos) << "H";
                }
            }

            out << "\t"
                << "*" // Reference name of the mate/next read
                << "\t"
                << "0" // Position of the mate/next read
                << "\t"
                << "0" // observed Template LENgth
                << "\t";

            // segment SEQuence
            if (no_seq_in_sam) {
                out << "*";
            } else {
                for (uint64_t p = query_start; p < query_end; ++p) {
                    out << query[p];
                }
            }

            out << "\t"
                << "*" // ASCII of Phred-scaled base QUALity+33
                << "\t"
                << "NM:i:"
                << edit_distance
                //<< "\t" << "AS:i:" << total_score
                << "\t"
                << "gi:f:" << gap_compressed_identity << "\t"
                << "bi:f:"
                << block_identity
                //<< "\t" << "md:f:" << mash_dist_sum / trace.size()
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "ii:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "dd:i:" << deleted_bp
                << "";

            if (emit_md_tag) {
                out << "\t";

                write_tag_and_md_string(out, cigarv, target_start);
            }

            out << "\t" << timings_and_num_alignements << "\n";
        }
    }

    // always clean up
    free(cigarv);
}

void write_alignment(
    std::ostream &out, const alignment_t &aln, const std::string &query_name,
    const uint64_t &query_total_length,
    const uint64_t &query_offset, // query offset on the forward strand
    const uint64_t &
        query_length, // used to compute the coordinates for reversed alignments
    const bool &query_is_rev, const std::string &target_name,
    const uint64_t &target_total_length, const uint64_t &target_offset,
    const uint64_t &target_length, // unused
    const float &min_identity, const float &mashmap_estimated_identity,
    const bool &with_endline) {

    if (aln.ok) {
        uint64_t matches = 0;
        uint64_t mismatches = 0;
        uint64_t insertions = 0;
        uint64_t inserted_bp = 0;
        uint64_t deletions = 0;
        uint64_t deleted_bp = 0;
        uint64_t refAlignedLength = 0;
        uint64_t qAlignedLength = 0;

        char *cigar = wfa_alignment_to_cigar(
            &aln.edit_cigar, refAlignedLength, qAlignedLength, matches,
            mismatches, insertions, inserted_bp, deletions, deleted_bp);

        size_t alignmentRefPos = aln.i;
        double gap_compressed_identity =
            (double)matches /
            (double)(matches + mismatches + insertions + deletions);
        double block_identity =
            (double)matches /
            (double)(matches + mismatches + inserted_bp + deleted_bp);
        // convert our coordinates to be relative to forward strand (matching
        // PAF standard)

        if (gap_compressed_identity >= min_identity) {
            uint64_t q_start;
            if (query_is_rev) {
                q_start =
                    query_offset + (query_length - (aln.j + qAlignedLength));
            } else {
                q_start = query_offset + aln.j;
            }
            out << query_name << "\t" << query_total_length << "\t" << q_start
                << "\t" << q_start + qAlignedLength << "\t"
                << (query_is_rev ? "-" : "+") << "\t" << target_name << "\t"
                << target_total_length << "\t"
                << target_offset + alignmentRefPos << "\t"
                << target_offset + alignmentRefPos + refAlignedLength << "\t"
                << matches << "\t" << std::max(refAlignedLength, qAlignedLength)
                << "\t" << std::round(float2phred(1.0 - block_identity)) << "\t"
                //<< "as:i:" << aln.score << "\t"
                << "gi:f:" << gap_compressed_identity << "\t"
                << "bi:f:"
                << block_identity
                //<< "\t" << "md:f:" << aln.mash_dist
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "bi:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "bd:i:" << deleted_bp
                << "\t"
                << "cg:Z:" << cigar << "\t"
                << "md:f:" << mashmap_estimated_identity;
            if (with_endline) {
                out << std::endl;
            }
        }
        free(cigar);
    }
}

char *alignment_to_cigar(const std::vector<char> &edit_cigar,
                         const uint64_t &start_idx, const uint64_t &end_idx,
                         uint64_t &target_aligned_length,
                         uint64_t &query_aligned_length, uint64_t &matches,
                         uint64_t &mismatches, uint64_t &insertions,
                         uint64_t &inserted_bp, uint64_t &deletions,
                         uint64_t &deleted_bp) {

    // the edit cigar contains a character string of ops
    // here we compress them into the standard cigar representation

    auto *cigar = new std::vector<char>();
    char lastMove = 0; // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;

    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int i = start_idx; i <= end_idx; i++) {
        // if new sequence of same moves started
        if (i == end_idx || (edit_cigar[i] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
            case 'M':
                matches += numOfSameMoves;
                query_aligned_length += numOfSameMoves;
                target_aligned_length += numOfSameMoves;
                break;
            case 'X':
                mismatches += numOfSameMoves;
                query_aligned_length += numOfSameMoves;
                target_aligned_length += numOfSameMoves;
                break;
            case 'I':
                ++insertions;
                inserted_bp += numOfSameMoves;
                query_aligned_length += numOfSameMoves;
                break;
            case 'D':
                ++deletions;
                deleted_bp += numOfSameMoves;
                target_aligned_length += numOfSameMoves;
                break;
            default:
                break;
            }

            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            std::reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            // reassign 'M' to '=' for convenience
            lastMove = lastMove == 'M' ? '=' : lastMove;
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < end_idx) {
                numOfSameMoves = 0;
            }
        }
        if (i < end_idx) {
            lastMove = edit_cigar[i];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0); // Null character termination.

    char *cigar_ = (char *)malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}

char *wfa_alignment_to_cigar(const wfa::cigar_t *const edit_cigar,
                             uint64_t &target_aligned_length,
                             uint64_t &query_aligned_length, uint64_t &matches,
                             uint64_t &mismatches, uint64_t &insertions,
                             uint64_t &inserted_bp, uint64_t &deletions,
                             uint64_t &deleted_bp) {

    // the edit cigar contains a character string of ops
    // here we compress them into the standard cigar representation

    auto *cigar = new std::vector<char>();
    char lastMove = 0; // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    const int start_idx = edit_cigar->begin_offset;
    const int end_idx = edit_cigar->end_offset;

    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int i = start_idx; i <= end_idx; i++) {
        // if new sequence of same moves started
        if (i == end_idx ||
            (edit_cigar->operations[i] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
            case 'M':
                matches += numOfSameMoves;
                query_aligned_length += numOfSameMoves;
                target_aligned_length += numOfSameMoves;
                break;
            case 'X':
                mismatches += numOfSameMoves;
                query_aligned_length += numOfSameMoves;
                target_aligned_length += numOfSameMoves;
                break;
            case 'I':
                ++insertions;
                inserted_bp += numOfSameMoves;
                query_aligned_length += numOfSameMoves;
                break;
            case 'D':
                ++deletions;
                deleted_bp += numOfSameMoves;
                target_aligned_length += numOfSameMoves;
                break;
            default:
                break;
            }

            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            std::reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            // reassign 'M' to '=' for convenience
            lastMove = lastMove == 'M' ? '=' : lastMove;
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < end_idx) {
                numOfSameMoves = 0;
            }
        }
        if (i < end_idx) {
            lastMove = edit_cigar->operations[i];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0); // Null character termination.

    char *cigar_ = (char *)malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}

double float2phred(const double &prob) {
    if (prob == 1)
        return 255; // guards against "-0"
    double p = -10 * (double)log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
}

void wflign_edit_cigar_copy(wfa::cigar_t *const edit_cigar_dst,
                            wfa::cigar_t *const edit_cigar_src) {
    edit_cigar_dst->max_operations = edit_cigar_src->max_operations;
    edit_cigar_dst->begin_offset = 0;
    edit_cigar_dst->end_offset =
        edit_cigar_src->end_offset - edit_cigar_src->begin_offset;
    edit_cigar_dst->score = edit_cigar_src->score;
    // alloc our ops
    edit_cigar_dst->operations = (char *)malloc(edit_cigar_dst->end_offset);
    memcpy(edit_cigar_dst->operations,
           edit_cigar_src->operations + edit_cigar_src->begin_offset,
           edit_cigar_dst->end_offset);
}

void sort_indels(std::vector<char> &v) {
    auto f = v.begin();
    while (f != v.end()) {
        auto j = f;
        while (j != v.end() && (*j == 'D' || *j == 'I')) {
            ++j;
        }
        if (j != f) {
            std::sort(f, j, [](char a, char b) { return b < a; });
            f = j;
        } else {
            ++f;
        }
    }
}

} // namespace wavefront

} // namespace wflign
