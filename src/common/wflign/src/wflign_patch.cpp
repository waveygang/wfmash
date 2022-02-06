#include <stddef.h>
#include <chrono>
#include <cstdlib>
#include <iterator>
#include <string>
#include "rkmh.hpp"
#include "wflign_patch.hpp"

namespace wflign {

namespace wavefront {

//wfa::wavefront_aligner_t* get_wavefront_aligner(
//    const wfa::affine_penalties_t& wfa_affine_penalties,
//    //ToDo: to remove if wavefront_aligner_new will not re-take the seqs' lens in input
//    const uint64_t& target_length,
//    const uint64_t& query_length,
//    const bool& low_memory) {
//    // Configure the attributes of the wf-aligner
//    wfa::wavefront_aligner_attr_t attributes =
//        wfa::wavefront_aligner_attr_default;
//    attributes.distance_metric = wfa::gap_affine;
//    attributes.affine_penalties = wfa_affine_penalties;
//    // attributes.distance_metric = gap_affine2p;
//    // attributes.affine2p_penalties = affine2p_penalties;
//    attributes.reduction.reduction_strategy =
//        wfa::wavefront_reduction_none; // wavefront_reduction_dynamic
//    // attributes.reduction.min_wavefront_length = 10;
//    // attributes.reduction.max_distance_threshold = 50;
//    attributes.alignment_scope =
//        wfa::compute_alignment; // alignment_scope_score
//    attributes.low_memory = low_memory;
//    //wfa::wavefront_aligner_t *const wf_aligner =
//    //return wfa::wavefront_aligner_new(target_length, query_length, &attributes);
//    return wfa::wavefront_aligner_new(&attributes);
//}

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
		const uint64_t& j,
		const std::string& target_name,
		const char* target,
		std::vector<rkmh::hash_t>*& target_sketch,
		const uint64_t& target_length,
		const uint64_t& i,
		const uint16_t& segment_length_q,
		const uint16_t& segment_length_t,
		const uint16_t& step_size,
		const uint64_t& minhash_kmer_size,
		const int& min_wavefront_length,
		const int& max_distance_threshold,
		const float& max_mash_dist,
		const float& mashmap_estimated_identity,
		wfa::WFAlignerGapAffine& wf_aligner,
		const wflign_penalties_t& affine_penalties,
		alignment_t& aln) {

    // first make the sketches if we haven't yet
    if (query_sketch == nullptr) {
        query_sketch = new std::vector<rkmh::hash_t>();
        *query_sketch = rkmh::hash_sequence(
            query + j, segment_length_q, minhash_kmer_size, segment_length_q / 8);
    }
    if (target_sketch == nullptr) {
        target_sketch = new std::vector<rkmh::hash_t>();
        *target_sketch = rkmh::hash_sequence(
            target + i, segment_length_t, minhash_kmer_size, segment_length_t / 8);
    }

    // first check if our mash dist is inbounds
    const float mash_dist =
        rkmh::compare(*query_sketch, *target_sketch, minhash_kmer_size);

    // this threshold is set low enough that we tend to randomly sample wflambda
    // matrix cells for alignment the threshold is adaptive, based on the mash
    // distance of the mapping we are aligning we should obtain enough
    // alignments that we can still patch between them
    if (mash_dist > max_mash_dist) {
        // if it isn't, return false
        return false;
    } else {
        // if it is, we'll align

        const int max_score = std::max(segment_length_q, segment_length_t) * (0.75 + mash_dist);
        wf_aligner.setMaxAlignmentScore(max_score);
        const int status = wf_aligner.alignEnd2End(
        		target + i,segment_length_t,
                query + j,segment_length_q);

        aln.j = j;
        aln.i = i;

        // aln.mash_dist = mash_dist;
        aln.ok = (status == 0);

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

            wflign_edit_cigar_copy(wf_aligner,&aln.edit_cigar);

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
		const int& segment_length,
		const int& min_wavefront_length,
		const int& max_distance_threshold,
		wfa::WFAlignerGapAffine& _wf_aligner,
		const wflign_penalties_t& affine_penalties,
		alignment_t& aln) {
    const long max_seg_len = 3 * segment_length;
    const bool big_wave = (query_length > max_seg_len || target_length > max_seg_len);
    wfa::WFAlignerGapAffine* wf_aligner = &_wf_aligner;
    if (big_wave) {
    	wf_aligner = new wfa::WFAlignerGapAffine(
    			affine_penalties.mismatch,
    			affine_penalties.gap_opening,
    			affine_penalties.gap_extension,
    			wfa::WFAligner::Alignment,
				wfa::WFAligner::MemoryMed);
    }

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
    	wf_aligner->setReductionNone();
    } else {
    	wf_aligner->setReductionAdaptive(min_wavefront_length,max_distance_threshold);
    }
    const int max_score = 2 * std::max(target_length, query_length);
    wf_aligner->setMaxAlignmentScore(max_score);
    const int status = wf_aligner->alignEnd2End(target + i,target_length,query + j,query_length);

    aln.ok = (status == 0);
    if (aln.ok) {
        // No more necessary: correct X/M errors in the cigar
        // hack_cigar(wf_aligner->cigar, query, target, query_length, target_length, j, i);

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

        wflign_edit_cigar_copy(*wf_aligner,&aln.edit_cigar);
    }

    if (big_wave) {
        delete wf_aligner;
    }
}
void write_merged_alignment(
    std::ostream &out,
    const std::vector<alignment_t *> &trace,
    wfa::WFAlignerGapAffine& wf_aligner,
    const wflign_penalties_t& affine_penalties,
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
    const uint16_t& segment_length,
    const float& min_identity,
    const long& elapsed_time_wflambda_ms,
    const uint64_t& num_alignments,
    const uint64_t& num_alignments_performed,
    const float& mashmap_estimated_identity,
    const uint64_t& wflign_max_len_major,
    const uint64_t& wflign_max_len_minor,
    const uint16_t& erode_k,
    const int& min_wf_length,
    const int& max_dist_threshold,
    const bool& with_endline) {

    int64_t target_pointer_shift = 0;

    uint64_t target_length_mut = target_length;

    // patching parameters
    // we will nibble patching back to this length
    const uint64_t min_wfa_patch_length = 128;

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

                    	wfa::WFAlignerGapAffine* wf_aligner_heads =
                    			new wfa::WFAlignerGapAffine(
                    					affine_penalties.mismatch,
                    					affine_penalties.gap_opening,
                    					affine_penalties.gap_extension,
                    					wfa::WFAligner::Alignment,
										wfa::WFAligner::MemoryMed);
                    	wf_aligner_heads->setReductionAdaptive(min_wf_length,max_dist_threshold);
                    	const int status = wf_aligner_heads->alignEndsFree(
                    			target_rev.c_str(),target_rev.size(),0,0,
                    			query_rev.c_str(),query_rev.size(),0,query_rev.size());
                        if (status == 0) { // WF_ALIGN_SUCCESSFUL
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

                            char* cigar_ops;
                            int cigar_length;
                            wf_aligner_heads->getAlignmentCigar(&cigar_ops,&cigar_length);
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

                    if ((size_region_to_repatch > 0 || (query_delta > 0 && target_delta > 0) || (query_delta > 2 || target_delta > 2)) &&
						(query_delta < wflign_max_len_major && target_delta < wflign_max_len_major) &&
						(query_delta < wflign_max_len_minor || target_delta < wflign_max_len_minor)) {

                        int32_t distance_close_indels = (query_delta > 3 || target_delta > 3) ?
                            distance_close_big_enough_indels(std::max(query_delta, target_delta), q, unpatched) :
                            -1;
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
                            if (query_delta >= 10 && target_delta >= 10) {
                                alignment_t patch_aln;
                                // WFA is only global
                                do_wfa_patch_alignment(
                                    query, query_pos, query_delta,
                                    target - target_pointer_shift, target_pos,
                                    target_delta, segment_length,
                                    min_wf_length, max_dist_threshold,
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
                                        patched.push_back(patch_aln.edit_cigar.cigar_ops[i]);
                                    }
                                    // std::cerr << "\n";

                                    // Check if there are too many indels in the
                                    // patch
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

                    	wfa::WFAlignerGapAffine* wf_aligner_tails =
                    			new wfa::WFAlignerGapAffine(
                    					affine_penalties.mismatch,
                    					affine_penalties.gap_opening,
                    					affine_penalties.gap_extension,
                    					wfa::WFAligner::Alignment,
										wfa::WFAligner::MemoryMed);
                    	wf_aligner_tails->setReductionAdaptive(min_wf_length,max_dist_threshold);
                    	const int status = wf_aligner_tails->alignEndsFree(
                    			target - target_pointer_shift + target_pos, target_delta_x,0,0,
                    			query + query_pos, query_delta,0,query_delta);

                        if (status == 0) { // WF_ALIGN_SUCCESSFUL
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

                            char* cigar_ops;
                            int cigar_length;
                            wf_aligner_tails->getAlignmentCigar(&cigar_ops,&cigar_length);
                            for(int xxx = 0; xxx < cigar_length; ++xxx) {
                                //std::cerr << cigar_ops[xxx];
                                patched.push_back(cigar_ops[xxx]);
                            }
                            //std::cerr << "\n";
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

            // std::cerr << "FIRST PATCH ROUND
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            patching(erodev, pre_tracev);

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
            sort_indels(pre_tracev);
        }

        // std::cerr << "SECOND PATCH ROUND
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        patching(pre_tracev, tracev);
    }

    // normalize the indels
    sort_indels(tracev);

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
    const bool& with_endline) {

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
