#include <cstddef>
#include <chrono>
#include <cstdlib>
#include <string>
#include <atomic_image.hpp>
#include "rkmh.hpp"
#include "wflign_patch.hpp"

namespace wflign {

    void encodeOneStep(const char *filename, std::vector<unsigned char> &image, unsigned width, unsigned height) {
        //Encode the image
        unsigned error = lodepng::encode(filename, image, width, height);

        //if there's an error, display it
        if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    }

    namespace wavefront {

// accumulate alignment objects
// run the traceback determine which are part of the main chain
// order them and write them out
// needed--- 0-cost deduplication of alignment regions (how????)
//     --- trim the alignment back to the first 1/2 of the query
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
        alignment_t& aln) {

    // if our i or j index plus segment length in the query or target is too long we'll make a memory access error and weird stuff will happen
    if (i + segment_length_t > target_length || j + segment_length_q > query_length) {
        // display function parameters
        std::cerr << "query_name: " << query_name << " query_length: " << query_length << " target_name: " << target_name << " target_length: " << target_length << std::endl;
        std::cerr << "i: " << i << " j: " << j << " segment_length_t: " << segment_length_t << " segment_length_q: " << segment_length_q << std::endl;
    }
    
    // first make the sketches if we haven't yet
    if (query_sketch == nullptr) {
        query_sketch = new std::vector<rkmh::hash_t>();
        *query_sketch = rkmh::hash_sequence(
                query + j, segment_length_q, extend_data->minhash_kmer_size, (uint64_t)((float)segment_length_q * extend_data->mash_sketch_rate));
        ++extend_data->num_sketches_allocated;
    }
    if (target_sketch == nullptr) {
        target_sketch = new std::vector<rkmh::hash_t>();        
        *target_sketch = rkmh::hash_sequence(
                target + i, segment_length_t, extend_data->minhash_kmer_size, (uint64_t)((float)segment_length_t * extend_data->mash_sketch_rate));
        ++extend_data->num_sketches_allocated;
    }

    // first check if our mash dist is inbounds
    const float mash_dist =
            rkmh::compare(*query_sketch, *target_sketch, extend_data->minhash_kmer_size);
    //std::cerr << "mash_dist is " << mash_dist << std::endl;

    // this threshold is set low enough that we tend to randomly sample wflambda
    // matrix cells for alignment the threshold is adaptive, based on the mash
    // distance of the mapping we are aligning we should obtain enough
    // alignments that we can still patch between them
    if (mash_dist > extend_data->max_mash_dist_to_evaluate) {
        // if it isn't, return false
        return false;
    } else {
        // if it is, we'll align
        const int max_score = (int)((float)std::max(segment_length_q, segment_length_t) * extend_data->inception_score_max_ratio);

        extend_data->wf_aligner->setMaxAlignmentSteps(max_score);
        const int status = extend_data->wf_aligner->alignEnd2End(
                target + i,segment_length_t,
                query + j,segment_length_q);

        aln.j = j;
        aln.i = i;

        aln.ok = (status == WF_STATUS_ALG_COMPLETED);

        // fill the alignment info if we aligned
        if (aln.ok) {
            aln.query_length = segment_length_q;
            aln.target_length = segment_length_t;

            /*
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
             */

            wflign_edit_cigar_copy(*extend_data->wf_aligner,&aln.edit_cigar);

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
        const int& max_patching_score) {

    const int max_score =
        max_patching_score ? max_patching_score :
            convex_penalties.gap_opening2 +
            (convex_penalties.gap_extension1 * std::min(
                    (int)chain_gap,
                    (int)std::max(target_length, query_length)
                )) + 64;

    /*
    std::cerr << "doing wfa patch alignment with parameters:"
              << " query_length " << query_length
              << " target_length " << target_length
              << " chain_gap " << chain_gap
              << " max_patching_score " << max_patching_score
              << " max_score " << max_score
              << std::endl
              << "query " << std::string(query + j, query_length)
              << " target " << std::string(target + i, target_length)
              << std::endl;
    */

    wf_aligner.setMaxAlignmentSteps(max_score);

    int fwd_score = std::numeric_limits<int>::max();
    int rev_score = std::numeric_limits<int>::max();
    
    const int status = wf_aligner.alignEnd2End(target + i, target_length, query + j, query_length);
    aln.ok = (status == WF_STATUS_ALG_COMPLETED);

    //std::cerr << "score is " << wf_aligner.getAlignmentScore() << std::endl;

    if (aln.ok) {
#ifdef VALIDATE_WFA_WFLIGN
        if (!validate_cigar(wf_aligner.cigar, query + j, target + i, query_length,
                            target_length, 0, 0)) {
            std::cerr << "cigar failure at alignment " << aln.j << " " << aln.i << std::endl;
            unpack_display_cigar(wf_aligner.cigar, query + j, target + i, query_length,
                                 target_length, 0, 0);
            std::cerr << ">query" << std::endl
                      << std::string(query + j, query_length) << std::endl;
            std::cerr << ">target" << std::endl
                      << std::string(target + i, target_length) << std::endl;
            assert(false);
        }
#endif

        wflign_edit_cigar_copy(wf_aligner, &aln.edit_cigar);
        fwd_score = calculate_alignment_score(aln.edit_cigar, convex_penalties);
        //std::cerr << "forward score is " << fwd_score << std::endl;
    }
    // if the score is -int32 max, we can try to reverse complement the query
    //     auto fwd_score = wf_aligner.getAlignmentScore(); // this is broken
    // compute with
    //int calculate_alignment_score(const char* cigar, const wflign_penalties_t& penalties) {
    //if (wf_aligner.getAlignmentScore() == -2147483648) {
    // Try reverse complement alignment
    std::string rev_comp_query = reverse_complement(std::string(query + j, query_length));
    const int rev_status = wf_aligner.alignEnd2End(target + i, target_length, rev_comp_query.c_str(), query_length);

    //auto rev_score = wf_aligner.getAlignmentScore();
    //rev_aln.ok = (rev_score > fwd_score && rev_status == WF_STATUS_ALG_COMPLETED);
    rev_aln.ok = (rev_status == WF_STATUS_ALG_COMPLETED);

    if (rev_aln.ok) {
        //std::cerr << "reverse complement alignment worked!" << std::endl;
#ifdef VALIDATE_WFA_WFLIGN
        if (!validate_cigar(wf_aligner.cigar, rev_comp_query.c_str(), target + i, query_length,
                            target_length, 0, 0)) {
            std::cerr << "cigar failure at reverse complement alignment " << j << " " << i << std::endl;
            unpack_display_cigar(wf_aligner.cigar, rev_comp_query.c_str(), target + i, query_length,
                                 target_length, 0, 0);
            std::cerr << ">query (reverse complement)" << std::endl
                      << rev_comp_query << std::endl;
            std::cerr << ">target" << std::endl
                      << std::string(target + i, target_length) << std::endl;
            assert(false);
        }
#endif

        wflign_edit_cigar_copy(wf_aligner, &rev_aln.edit_cigar);
        rev_score = calculate_alignment_score(rev_aln.edit_cigar, convex_penalties);
        //std::cerr << "reverse score is " << rev_score << std::endl;
        rev_aln.j = j;
        rev_aln.i = i;
        rev_aln.query_length = query_length;
        rev_aln.target_length = target_length;
    }

    if (rev_aln.ok && rev_score < fwd_score) {
        rev_aln.ok = true;
        aln.ok = false;
        std::cerr << "got better score with reverse complement alignment" << std::endl
              << " query_length " << query_length
              << " target_length " << target_length
              << " chain_gap " << chain_gap
              << " max_patching_score " << max_patching_score
              << " max_score " << max_score
              << " fwd_score " << fwd_score
              << " rev_score " << rev_score
              << std::endl
              << "query " << std::string(query + j, query_length)
              << " target " << std::string(target + i, target_length)
              << std::endl;
    } else {
        rev_aln.ok = false;
    }

    aln.j = j;
    aln.i = i;
    aln.query_length = query_length;
    aln.target_length = target_length;
}

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
        const int& min_wf_length,
        const int& max_dist_threshold,
#ifdef WFA_PNG_TSV_TIMING
        const std::string* prefix_wavefront_plot_in_png,
        const uint64_t& wfplot_max_size,
        const bool& emit_patching_tsv,
        std::ostream* out_patching_tsv,
#endif
        const bool& with_endline) {

    int64_t target_pointer_shift = 0;

    uint64_t target_length_mut = target_length;

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
    //uint64_t total_score = 0;

    // double mash_dist_sum = 0;
    uint64_t ok_alns = 0;

    auto start_time = std::chrono::steady_clock::now();

#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] processing traceback"
      << std::endl;
#endif
    // write trace into single cigar vector
    std::vector<char> tracev;
    std::vector<alignment_t> rev_patch_alns;
    {
        // patch: walk the cigar, patching directly when we have simultaneous
        // gaps in query and ref and adding our results to the final trace as we
        // go

#define MAX_NUM_INDELS_TO_LOOK_AT 2
        auto distance_close_big_enough_indels =	
                [](const uint32_t indel_len, auto iterator,	
                   const std::vector<char> &trace,
                   const uint16_t&max_dist_to_look_at) {	
                    const uint32_t min_indel_len_to_find = indel_len / 3;	

                    // std::cerr << "min_indel_len_to_find " <<	
                    // min_indel_len_to_find << std::endl; std::cerr <<	
                    // "max_dist_to_look_at " << max_dist_to_look_at << std::endl;	

                    auto q = iterator;	

                    uint8_t num_indels_to_find = MAX_NUM_INDELS_TO_LOOK_AT;	
                    uint32_t curr_size_close_indel = 0;	
                    int32_t dist_close_indels = 0;	

                    // std::cerr << "indel_len " << indel_len << std::endl;	
                    // std::cerr << "min_indel_len_to_find " << min_indel_len_to_find << std::endl;	
                    // std::cerr << "max_dist_to_look_at " << max_dist_to_look_at << std::endl;	

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

                    // std::cerr << "dist_close_indels " << dist_close_indels << std::endl;	
                    // std::cerr << "num_indels_found " << MAX_NUM_INDELS_TO_LOOK_AT - num_indels_to_find << std::endl;	

                    return num_indels_to_find < MAX_NUM_INDELS_TO_LOOK_AT	
                           ? dist_close_indels	
                           : -1;	
                };

        auto patching = [&query, &query_name, &query_length, &query_start,
                &query_offset, &target, &target_name,
                &target_length_mut, &target_start, &target_offset,
                &target_total_length, &target_end,
                &target_pointer_shift,
                &wflign_max_len_major,
                &wflign_max_len_minor,
                &distance_close_big_enough_indels, &min_wf_length,
                &max_dist_threshold, &wf_aligner,
                &rev_patch_alns,
                &convex_penalties,
                &chain_gap, &max_patching_score
#ifdef WFA_PNG_TSV_TIMING
                ,&emit_patching_tsv,
                &out_patching_tsv
#endif
        ](std::vector<char> &unpatched,
                                   std::vector<char> &patched,
                                   const uint16_t &min_wfa_head_tail_patch_length,
                                   const uint16_t &min_wfa_patch_length,
                                   const uint16_t &max_dist_to_look_at) {
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

                    // nibble forward if we're below the correct length
                    // this gives a bit of context for the alignment
                    while (q != unpatched.end() &&
                           (query_delta < min_wfa_head_tail_patch_length || target_delta < min_wfa_head_tail_patch_length)) {
                        const auto &c = *q++;
                        switch (c) {
                            case 'M': case 'X':
                                ++query_delta; ++target_delta; break;
                            case 'I': ++query_delta; break;
                            case 'D': ++target_delta; break;
                            default: break;
                        }
                    }

//                    std::cerr << "A HEAD patching in"
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
                    const uint64_t delta_to_ask = query_delta <= target_delta ? 0 : query_delta - target_delta;

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
                    //    std::cerr << "B HEAD patching in "
                    //              << query_name << " "
                    //              << query_offset + query_pos <<
                    //              " - " << query_delta
                    //            << " --- "
                    //            << target_name << " " 
                    //            << target_offset + target_pos << " - "
                    //            << target_delta_x
                    //            << std::endl;

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

                        wfa::WFAlignerGapAffine2Pieces* wf_aligner_heads =
                                new wfa::WFAlignerGapAffine2Pieces(
                                        0,
                                        convex_penalties.mismatch,
                                        convex_penalties.gap_opening1,
                                        convex_penalties.gap_extension1,
                                        convex_penalties.gap_opening2,
                                        convex_penalties.gap_extension2,
                                        wfa::WFAligner::Alignment,
                                        wfa::WFAligner::MemoryMed);
                        wf_aligner_heads->setHeuristicWFmash(min_wf_length,max_dist_threshold);
                        const int status = wf_aligner_heads->alignEndsFree(
                                target_rev.c_str(),target_rev.size(),0,0,
                                query_rev.c_str(),query_rev.size(),0,query_rev.size());
                        //std::cerr << "Head patching status " << status << "\n";
                        if (status == WF_STATUS_ALG_COMPLETED || status == WF_STATUS_ALG_PARTIAL) {
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

                            got_alignment = true;

                            target_pos = target_pos_x;
                            target_delta = target_delta_x;
                            target_pointer_shift = target_pointer_shift_x;

                            target_start = target_start_x;
                            target_length_mut += target_delta_to_shift;

                            char* cigar_ops;
                            int cigar_length;
                            wf_aligner_heads->getAlignment(&cigar_ops,&cigar_length);

                            // Put the missing part of the CIGAR string first and then its aligned part
                            uint64_t missing_query_len = query_rev.size();
                            uint64_t missing_target_len = target_rev.size();
                            for(int xxx = cigar_length - 1; xxx >= 0; --xxx) {
                                // The CIGAR string can be incomplete
                                switch (cigar_ops[xxx]) {
                                    case 'M': case 'X':
                                        --missing_query_len; --missing_target_len; break;
                                    case 'I': --missing_query_len; break;
                                    case 'D': --missing_target_len; break;
                                    default: break;
                                }
                            }
                            // Put the missing part of the CIGAR string
                            for (uint64_t i = 0; i < missing_query_len; ++i) {
                                patched.push_back('I');
                            }
                            for (uint64_t i = 0; i < missing_target_len; ++i) {
                                patched.push_back('D');
                            }

                            for(int xxx = cigar_length - 1; xxx >= 0; --xxx) {
                                //std::cerr << cigar_ops[xxx];
                                patched.push_back(cigar_ops[xxx]);
                            }
                            //std::cerr << "\n";
                        }
                        delete wf_aligner_heads;
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
            while (q != unpatched.end()) {
                // get to the first match
                while (q != unpatched.end() && (*q == 'M' || *q == 'X')) {
                    /*
                std::cerr << "q: " << query[query_pos] << " "
                          << "t: " << target[target_pos - target_pointer_shift]
                          << std::endl;
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
                // unused!

                {
                    got_alignment = false;

                    if ((size_region_to_repatch > 0 ||
                         (query_delta > 0 && target_delta > 0) ||
                         (query_delta > 2 || target_delta > 2) &&
                                 (query_delta < wflign_max_len_major &&
                                  target_delta < wflign_max_len_major) &&
                                 (query_delta < wflign_max_len_minor ||
                                  target_delta < wflign_max_len_minor))) {

                        int32_t distance_close_indels = 
                            (query_delta > 10 || target_delta > 10) ?	
                            distance_close_big_enough_indels(std::max(query_delta, target_delta), q, unpatched, max_dist_to_look_at) :	
                            -1;

                        // Trigger the patching if there is a dropout
                        // (consecutive Is and Ds) or if there is a close and
                        // big enough indel forward
                        if (size_region_to_repatch > 0 ||
                            (query_delta > 0 && target_delta > 0) ||
                            (query_delta > 2 || target_delta > 2) ||
                            distance_close_indels > 0) {
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

                            size_region_to_repatch = 0;
                            {
                                alignment_t patch_aln;
                                alignment_t rev_patch_aln;
                                // WFA is only global
                                do_wfa_patch_alignment(
                                        query, query_pos, query_delta,
                                        target - target_pointer_shift, target_pos, target_delta,
                                        wf_aligner, convex_penalties, patch_aln, rev_patch_aln,
                                        chain_gap, max_patching_score);
                                if (rev_patch_aln.ok) {
                                    // we got a good reverse alignment
                                    rev_patch_alns.push_back(rev_patch_aln);
                                }
                                if (patch_aln.ok) {
                                    // std::cerr << "got an ok patch aln" <<
                                    // std::endl;
                                    got_alignment = true;
                                    const int start_idx =
                                            patch_aln.edit_cigar.begin_offset;
                                    const int end_idx =
                                            patch_aln.edit_cigar.end_offset;
                                    for (int i = start_idx; i < end_idx; i++) {
                                        patched.push_back(patch_aln.edit_cigar.cigar_ops[i]);
                                    }
                                    // std::cerr << "\n";

                                    // Check if there are too many indels in the
                                    // patch
                                    /*
                                    uint32_t size_indel = 0;
                                    for (int i = end_idx - 1; i >= start_idx;
                                         --i) {
                                        // std::cerr <<
                                        // patch_aln.edit_cigar.operations[i];
                                        if (patch_aln.edit_cigar.cigar_ops[i] == 'I' ||
                                            patch_aln.edit_cigar.cigar_ops[i] == 'D') {
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
                                    */
                                    // std::cerr << std::endl;

                                    // Not too big, to avoid repatching
                                    // structural variants boundaries
                                    //std::cerr << "size_region_to_repatch " << size_region_to_repatch << std::endl;
                                    //std::cerr << "end_idx - start_idx " << end_idx - start_idx << std::endl;
                                    /*
                                    if (size_indel > 7 && size_indel <= 4096 &&
                                        size_region_to_repatch <
                                                (end_idx - start_idx)) {
                                        //std::cerr << "REPATCH " << std::endl;
                                    } else {
                                        size_region_to_repatch = 0;
                                    }
                                    */
                                }
#ifdef WFA_PNG_TSV_TIMING
                                if (emit_patching_tsv) {
                                    *out_patching_tsv
                                            << query_name << "\t" << query_pos << "\t" << query_pos + query_delta << "\t"
                                            << target_name << "\t" << (target_pos - target_pointer_shift) << "\t" << (target_pos - target_pointer_shift + target_delta) << "\t"
                                            << patch_aln.ok << std::endl;
                                }
#endif
                            }
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

                    // std::cerr << "query_delta " << query_delta << std::endl;
                    // std::cerr << "target_delta " << target_delta <<
                    // std::endl;
                    query_pos += query_delta;
                    target_pos += target_delta;

                    query_delta = 0;
                    target_delta = 0;
                }
            }

            // Tail patching
            {
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
                    // nibble backward if we're below the correct length
                    // this gives a bit of context for the alignment
                    while (!patched.empty() &&
                           (query_delta < min_wfa_head_tail_patch_length || target_delta < min_wfa_head_tail_patch_length)) {
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

//                    std::cerr << "A TAIL patching in "
//                              << query_name << " " <<
//                              query_offset << " @ " <<
//                              query_pos << " - " <<
//                              query_delta << " "
//                              << target_name << " " <<
//                              target_offset << " @ "
//                              << target_pos - target_pointer_shift << " - "
//                              << target_delta
//                              << std::endl;

                    // min_wfa_patch_length-bps of margins to manage insertions in the query
                    const uint64_t delta_to_ask = query_delta <= target_delta ? 0 : query_delta - target_delta;

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
                    //    std::cerr << "B TAIL patching in "
                    //              << query_name << " " <<
                    //              query_offset + query_pos << " - " <<
                    //              query_delta << " "
                    //              << target_name << " " <<
                    //              target_offset + target_pos - target_pointer_shift << " - "
                    //              << target_delta_x
                    //              << std::endl;

                        wfa::WFAlignerGapAffine2Pieces* wf_aligner_tails =
                                new wfa::WFAlignerGapAffine2Pieces(
                                        0,
                                        convex_penalties.mismatch,
                                        convex_penalties.gap_opening1,
                                        convex_penalties.gap_extension1,
                                        convex_penalties.gap_opening2,
                                        convex_penalties.gap_extension2,
                                        wfa::WFAligner::Alignment,
                                        wfa::WFAligner::MemoryMed);
                        wf_aligner_tails->setHeuristicWFmash(min_wf_length,max_dist_threshold);
                        const int status = wf_aligner_tails->alignEndsFree(
                                target - target_pointer_shift + target_pos, target_delta_x,0,0,
                                query + query_pos, query_delta,0,query_delta);
                        //std::cerr << "Tail patching status " << status << "\n";
                        if (status == WF_STATUS_ALG_COMPLETED || status == WF_STATUS_ALG_PARTIAL) {
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

                            char* cigar_ops;
                            int cigar_length;
                            wf_aligner_tails->getAlignment(&cigar_ops,&cigar_length);
                            uint64_t missing_query_len = query_delta;
                            uint64_t missing_target_len = target_delta;
                            for(int xxx = 0; xxx < cigar_length; ++xxx) {
                                //std::cerr << cigar_ops[xxx];
                                patched.push_back(cigar_ops[xxx]);

                                // The CIGAR string can be incomplete
                                switch (cigar_ops[xxx]) {
                                    case 'M': case 'X':
                                        --missing_query_len; --missing_target_len; break;
                                    case 'I': --missing_query_len; break;
                                    case 'D': --missing_target_len; break;
                                    default: break;
                                }
                            }
                            //std::cerr << "\n";
                            // Put the missing part of the CIGAR string
                            for (uint64_t i = 0; i < missing_query_len; ++i) {
                                patched.push_back('I');
                            }
                            for (uint64_t i = 0; i < missing_target_len; ++i) {
                                patched.push_back('D');
                            }
                        }
                        delete wf_aligner_tails;
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
                // not used
                // query_pos += query_delta;
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
                            const auto &c = aln.edit_cigar.cigar_ops[i];
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

            //std::cerr << "FIRST PATCH ROUND" << std::endl;
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            patching(erodev, pre_tracev, 4096, 32, 512); // In the 1st round, we patch more aggressively

#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_trace(pre_tracev, query,
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

        //std::cerr << "SECOND PATCH ROUND" << std::endl;
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        patching(pre_tracev, tracev, 256, 8, 128); // In the 2nd round, we patch less aggressively
    }

    // normalize the indels
    //sort_indels(tracev);

#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] got full patched traceback: ";
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
    uint64_t begin_offset;
    uint64_t end_offset;
    {
        uint64_t trim_del_first;
        uint64_t trim_del_last;

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

#ifdef WFA_PNG_TSV_TIMING
    bool emit_png = !prefix_wavefront_plot_in_png->empty() && wfplot_max_size > 0;
    if (emit_png) {
        const int pattern_length = (int)query_length;
        const int text_length = (int)target_length;

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
                            uint64_t _v = v;
                            uint64_t _h = h;
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
            const std::string filename = *prefix_wavefront_plot_in_png +
                                         query_name + "_" + std::to_string(query_offset) + "_" + std::to_string(query_offset+query_length) + " _ " + (query_is_rev ? "-" : "+") +
                                         "_" + target_name + "_" + std::to_string(target_offset) + "_" + std::to_string(target_offset+target_length) + ".2.trace.png";
            encodeOneStep(filename.c_str(), bytes, width, height);
        }
    }
#endif

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

#ifdef WFA_PNG_TSV_TIMING
        const long elapsed_time_patching_ms =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - start_time)
                        .count();

        const std::string timings_and_num_alignements =
                "wt:i:" + std::to_string(elapsed_time_wflambda_ms) +
                "\tpt:i:" + std::to_string(elapsed_time_patching_ms) +
                "\taa:i:" + std::to_string(num_alignments) +
                "\tap:i:" + std::to_string(num_alignments_performed);
#endif

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
                << matches + mismatches + inserted_bp + deleted_bp
                << "\t"
                << std::round(float2phred(1.0 - block_identity))
                //<< "\t" << "as:i:" << total_score
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
                << "\t"
                << "md:f:" << mashmap_estimated_identity;

            if (emit_md_tag) {
                out << "\t";

                write_tag_and_md_string(out, cigarv, target_start);
            }

#ifdef WFA_PNG_TSV_TIMING
            out << "\t" << timings_and_num_alignements << "\t"
                << "cg:Z:" << cigarv << "\n";
#else
            out << "\t" << "cg:Z:" << cigarv << "\n";
#endif
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
#ifdef WFA_PNG_TSV_TIMING
            out << "\t" << timings_and_num_alignements << "\n";
#else
            out << "\n";
#endif
        }
    }

    // always clean up
    free(cigarv);

    // write how many reverse complement alignments were found
    //std::cerr << "got " << rev_patch_alns.size() << " rev patch alns" << std::endl;
    for (auto& rev_patch_aln : rev_patch_alns) {
        bool rev_query_is_rev = !query_is_rev;  // Flip the orientation
        write_alignment(
            out,
            rev_patch_aln,
            query_name,
            query_total_length,
            query_offset,
            query_length,
            rev_query_is_rev,  // Use the flipped orientation
            target_name,
            target_total_length,
            target_offset,
            target_length,
            min_identity,
            mashmap_estimated_identity,
            false,  // Don't add an endline after each alignment
            true);  // This is a reverse complement alignment
        // write tag indicating that this is a reverse complement alignment
        out << "\t" << "rc:Z:true" << "\n";
    }
    out << std::flush;
}

void write_alignment(
        std::ostream& out,
        const alignment_t& aln,
        const std::string& query_name,
        const uint64_t& query_total_length,
        const uint64_t& query_offset, // query offset on the forward strand
        const uint64_t& query_length, // used to compute the coordinates for reversed alignments
        const bool& query_is_rev,
        const std::string& target_name,
        const uint64_t& target_total_length,
        const uint64_t& target_offset,
        const uint64_t& target_length, // unused
        const float& min_identity,
        const float& mashmap_estimated_identity,
        const bool& with_endline,
        const bool& is_rev_patch) {

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
            if (query_is_rev && !is_rev_patch) {
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

double float2phred(const double& prob) {
    if (prob == 1)
        return 255; // guards against "-0"
    double p = -10 * (double)log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
}

void sort_indels(std::vector<char>& v) {
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

    } /* namespace wavefront */
} /* namespace wflign */
