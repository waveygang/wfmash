#include <chrono>
#include "wflign_wfa.hpp"

namespace wflign {

namespace wavefront {

void wflign_affine_wavefront(
    std::ostream& out,
    const bool& merge_alignments,
    const bool& emit_md_tag,
    const bool& paf_format_else_sam,
    const std::string& query_name,
    const char* query,
    const uint64_t& query_total_length,
    const uint64_t& query_offset, // todo this is broken for reverse alignment reporting!!!!
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
    const int& wflambda_max_distance_threshold) {
    //const int& wfa_min_wavefront_length, // with these set at 0 we do exact WFA for WFA itself
    //const int& wfa_max_distance_threshold) {

    if (query_offset + query_length > query_total_length
        || target_offset + target_length > target_total_length) {
        return;
    }

    // set up our implicit matrix
    const uint64_t steps_per_segment = 2;
    const uint64_t step_size = segment_length / steps_per_segment;

    // Pattern & Text
    const int pattern_length = query_length / step_size;
    const int text_length = target_length / step_size;

    // uncomment to use reduced WFA locally
    // currently not supported due to issues with traceback when applying WF-reduction on small problems
    const int wfa_min_wavefront_length = 0; //segment_length / 16;
    const int wfa_max_distance_threshold = 0; //segment_length / 8;

    // Allocate MM
    wflambda::mm_allocator_t* const wflambda_mm_allocator = wflambda::mm_allocator_new(BUFFER_SIZE_8M);
    // Set penalties
    wflambda::affine_penalties_t wflambda_affine_penalties = {
        .match = 0,
        .mismatch = 7,
        .gap_opening = 11,
        .gap_extension = 1,
    };
    // Init Affine wflambda
    wflambda::affine_wavefronts_t* affine_wavefronts;
    if (wflambda_min_wavefront_length || wflambda_max_distance_threshold) {
        affine_wavefronts = wflambda::affine_wavefronts_new_reduced(
            pattern_length+1, text_length+1, &wflambda_affine_penalties,
            wflambda_min_wavefront_length, wflambda_max_distance_threshold,
            NULL, wflambda_mm_allocator);
    } else {
        affine_wavefronts = wflambda::affine_wavefronts_new_complete(
            pattern_length+1, text_length+1, &wflambda_affine_penalties, NULL, wflambda_mm_allocator);
    }

    // save computed alignments in a pair-indexed patchmap
    whash::patchmap<uint64_t,alignment_t*> alignments;

    // allocate vectors to store our sketches
    std::vector<std::vector<rkmh::hash_t>*> query_sketches(pattern_length, nullptr);
    std::vector<std::vector<rkmh::hash_t>*> target_sketches(text_length, nullptr);

    //std::cerr << "v" << "\t" << "h" << "\t" << "score" << "\t" << "aligned" << std::endl;

    // setup affine WFA
    wfa::mm_allocator_t* const wfa_mm_allocator = wfa::mm_allocator_new(BUFFER_SIZE_8M);
    wfa::affine_penalties_t wfa_affine_penalties = {
        .match = 0,
        .mismatch = 7,
        .gap_opening = 11,
        .gap_extension = 1,
    };
    const uint64_t minhash_kmer_size = 13;
    int v_max = 0;
    int h_max = 0;

    auto extend_match = [&](const int& v, const int& h) {
            bool aligned = false;
            if (v >= 0 && h >= 0 && v < pattern_length && h < text_length) {
                uint64_t k = encode_pair(v, h);
                auto f = alignments.find(k); //TODO: it can be removed using an edit-distance mode as high-level of WF-inception
                if (f != alignments.end()) {
                    aligned = true;
                } else  {
                    uint64_t query_begin = (v < pattern_length-1 ? v * step_size
                                            : query_length - segment_length);
                    uint64_t target_begin = (h < text_length-1 ? h * step_size
                                             : target_length - segment_length);
                    auto* aln = new alignment_t();
                    aligned =
                        do_wfa_segment_alignment(
                            query_name,
                            query,
                            query_sketches[v],
                            query_length,
                            query_begin,
                            target_name,
                            target,
                            target_sketches[h],
                            target_length,
                            target_begin,
                            segment_length,
                            step_size,
                            minhash_kmer_size,
                            wfa_min_wavefront_length,
                            wfa_max_distance_threshold,
                            wfa_mm_allocator,
                            &wfa_affine_penalties,
                            *aln);
                    //std::cerr << v << "\t" << h << "\t" << aln->score << "\t" << aligned << std::endl;
                    if (aligned) {
                        alignments[k] = aln;
                    } else {
                        delete aln;
                    }
                    // cleanup old sketches
                    if (v > v_max) {
                        v_max = v;
                        if (v >= wflambda_max_distance_threshold) {
                            auto& s = query_sketches[v - wflambda_max_distance_threshold];
                            // The C++ language guarantees that delete p will do nothing if p is equal to NULL
                            delete s;
                            s = nullptr;
                        }
                    }
                    if (h > h_max) {
                        h_max = h;
                        if (h >= wflambda_max_distance_threshold) {
                            auto& s = target_sketches[h - wflambda_max_distance_threshold];
                            // The C++ language guarantees that delete p will do nothing if p is equal to NULL
                            delete s;
                            s = nullptr;
                        }
                    }
                }
            } else if (h < 0 || v < 0) { //TODO: it can be removed using an edit-distance mode as high-level of WF-inception
                aligned = true;
            }
            return aligned;
        };

    // accumulate runs of matches in reverse order
    // then trim the cigars of successive mappings
    //
    std::vector<alignment_t*> trace;

    auto trace_match = [&](const int& v, const int& h) {
        if (v < 0 || h < 0 || v > pattern_length || h > text_length) {
            return false;
        } else {
            uint64_t k = encode_pair(v, h);
            trace.push_back(alignments[k]);
            return true;
        }
    };

    auto start_time = std::chrono::steady_clock::now();

    // Align
    wflambda::affine_wavefronts_align(
        affine_wavefronts,
        extend_match,
        trace_match,
        pattern_length,
        text_length);

    long elapsed_time_wflambda_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count();

//#define WFLIGN_DEBUG
#ifdef WFLIGN_DEBUG
    // get alignment score
    const int score = wflambda::edit_cigar_score_gap_affine(
        &affine_wavefronts->edit_cigar, &wflambda_affine_penalties);

    std::cerr << "[wflign::wflign_affine_wavefront] alignment score " << score << " for query: " << query_name << " target: " << target_name << std::endl;
#endif

    // clean up sketches
    // The C++ language guarantees that delete p will do nothing if p is equal to NULL
    for (auto& s : query_sketches) {
        delete s;
        s = nullptr;
    }
    for (auto& s : target_sketches) {
        delete s;
        s = nullptr;
    }

    // todo: implement alignment identifier based on hash of the input, params, and commit
    // annotate each PAF record with it and the full alignment score

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
        for (auto x = trace.rbegin()+1; x != trace.rend(); ++x) {
            auto& curr = **x;
            auto& last = **(x-1);
            /*
            std::cerr << "splicing" << std::endl;
            std::cerr << "curr = ";
            curr.display();
            std::cerr << "last = ";
            last.display();
            */
#ifdef VALIDATE_WFA_WFLIGN
            if (curr.ok && !curr.validate(query, target)) {
                std::cerr << "curr traceback is wrong before trimming @ " << curr.j << " " << curr.i << std::endl;
                curr.display();
                assert(false);
            }
#endif
            trace_pos_t last_pos = { last.j, last.i,
                                     &last.edit_cigar,
                                     last.edit_cigar.begin_offset };
            trace_pos_t curr_pos = { curr.j, curr.i,
                                     &curr.edit_cigar,
                                     curr.edit_cigar.begin_offset };

            // trace the last alignment until we overlap the next
            // to record our match
            trace_pos_t match_pos;

            // walk until they are matched at the query position
            while (!last_pos.at_end() && !curr_pos.at_end()) {
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
            int trim_last=0, trim_curr=0;
            if (match_pos.assigned()) {
                // we'll use our match position to set up the trims
                trim_last = (last.j + last.query_length) - match_pos.j;
                trim_curr = match_pos.j - curr.j;
            } else {
                // we want to remove any possible overlaps in query or target
                // walk back last until we don't overlap in i or j
                // recording the distance walked as an additional trim on last
                bool flip = false;
                while (last_pos.j > curr_pos.j
                       || last_pos.i > curr_pos.i) {
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
            if (trim_last > 0) {
                last.trim_back(trim_last);
#ifdef VALIDATE_WFA_WFLIGN
                if (last.ok && !last.validate(query, target)) {
                    std::cerr << "traceback is wrong after last trimming @ " << last.j << " " << last.i << std::endl;
                    last.display();
                    assert(false);
                }
#endif
            }

            if (trim_curr > 0) {
                curr.trim_front(trim_curr);
#ifdef VALIDATE_WFA_WFLIGN
                if (curr.ok && !curr.validate(query, target)) {
                    std::cerr << "traceback is wrong after curr trimming @ " << curr.j << " " << curr.i << std::endl;
                    curr.display();
                    assert(false);
                }
#endif
            }
            /*
            std::cerr << "spliced" << std::endl;
            std::cerr << "curr = ";
            curr.display();
            std::cerr << "last = ";
            last.display();
            std::cerr << std::endl;
            */
        }

        if (merge_alignments) {
            // write a merged alignment
            write_merged_alignment(out, trace, wfa_mm_allocator, &wfa_affine_penalties,
                                   emit_md_tag, paf_format_else_sam,
                                   query,
                                   query_name, query_total_length, query_offset, query_length,
                                   query_is_rev,
                                   target,
                                   target_name, target_total_length, target_offset, target_length,
                                   min_identity,
                                   segment_length * 128,
                                   elapsed_time_wflambda_ms);
        } else {
            for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
                //std::cerr << "on alignment" << std::endl;
                write_alignment(out, **x,
                                query_name, query_total_length, query_offset, query_length,
                                query_is_rev,
                                target_name, target_total_length, target_offset, target_length,
                                min_identity);
            }
        }
    }

    // clean up our WFA allocator
    wfa::mm_allocator_delete(wfa_mm_allocator);

    // Free
    wflambda::affine_wavefronts_delete(affine_wavefronts);
    wflambda::mm_allocator_delete(wflambda_mm_allocator);
    for (const auto& p : alignments) {
        delete p.second;
    }
}

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
    const uint64_t& segment_length,
    const uint64_t& step_size,
    const uint64_t& minhash_kmer_size,
    const int min_wavefront_length,
    const int max_distance_threshold,
    wfa::mm_allocator_t* const mm_allocator,
    wfa::affine_penalties_t* const affine_penalties,
    alignment_t& aln) {

    aln.query_length = segment_length;
    aln.target_length = segment_length;

    // first make the sketches if we haven't yet
    if (query_sketch == nullptr) {
        query_sketch = new std::vector<rkmh::hash_t>();
        *query_sketch = rkmh::hash_sequence(query+j, segment_length, minhash_kmer_size, segment_length/8);
    }
    if (target_sketch == nullptr) {
        target_sketch = new std::vector<rkmh::hash_t>();
        *target_sketch = rkmh::hash_sequence(target+i, segment_length, minhash_kmer_size, segment_length/8);
    }

    // first check if our mash dist is inbounds
    const double mash_dist = rkmh::compare(*query_sketch, *target_sketch, minhash_kmer_size);

    const int max_score = segment_length;

    // the mash distance generally underestimates the actual divergence
    // but when it's high we are almost certain that it's not a match
    if (mash_dist > 0.618034) {
        // if it isn't, return false
        aln.score = max_score;
        aln.ok = false;
        return false;
    } else {
        // if it is, we'll align
        wfa::affine_wavefronts_t* affine_wavefronts;
        if (min_wavefront_length || max_distance_threshold) {
            // adaptive affine WFA setup
            affine_wavefronts = affine_wavefronts_new_reduced(
                segment_length, segment_length, affine_penalties,
                min_wavefront_length, max_distance_threshold,
                NULL, mm_allocator);
        } else {
            // exact WFA
            affine_wavefronts = affine_wavefronts_new_complete(
                segment_length, segment_length, affine_penalties, NULL, mm_allocator);
        }

        aln.score = wfa::affine_wavefronts_align_bounded(
            affine_wavefronts,
            target+i,
            segment_length,
            query+j,
            segment_length,
            max_score);

        aln.j = j;
        aln.i = i;
        aln.mash_dist = mash_dist;
        // copy our edit cigar if we aligned
        aln.ok = aln.score < max_score;
        if (aln.ok) {
#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(affine_wavefronts->edit_cigar, query, target, segment_length, segment_length, aln.j, aln.i)) {
                std::cerr << "cigar failure at alignment " << aln.j << " " << aln.i << std::endl;
                unpack_display_cigar(affine_wavefronts->edit_cigar, query, target, segment_length, segment_length, aln.j, aln.i);
                std::cerr << ">query" << std::endl << std::string(query+j, segment_length) << std::endl;
                std::cerr << ">target" << std::endl << std::string(target+i, segment_length) << std::endl;
                assert(false);
            }
#endif
            wflign_edit_cigar_copy(&aln.edit_cigar, &affine_wavefronts->edit_cigar);
#ifdef VALIDATE_WFA_WFLIGN
            if (!validate_cigar(affine_wavefronts->edit_cigar, query, target, segment_length, segment_length, aln.j, aln.i)) {
                std::cerr << "cigar failure after cigar copy in alignment " << aln.j << " " << aln.i << std::endl;
                assert(false);
            }
#endif
        }
        // cleanup wavefronts to keep memory low
        affine_wavefronts_delete(affine_wavefronts);

        return aln.ok;
    }
}

void do_wfa_patch_alignment(
    const char* query,
    const uint64_t& j,
    const uint64_t& query_length,
    const char* target,
    const uint64_t& i,
    const uint64_t& target_length,
    const int min_wavefront_length,
    const int max_distance_threshold,
    wfa::mm_allocator_t* const mm_allocator,
    wfa::affine_penalties_t* const affine_penalties,
    alignment_t& aln) {

    //std::cerr << "do_wfa_patch " << j << " " << query_length << " " << i << " " << target_length << std::endl;
    wfa::mm_allocator_clear(mm_allocator);

    wfa::affine_wavefronts_t* affine_wavefronts;
    if (min_wavefront_length || max_distance_threshold) {
        // adaptive affine WFA setup
        affine_wavefronts = affine_wavefronts_new_reduced(
            target_length, query_length, affine_penalties,
            min_wavefront_length, max_distance_threshold,
            NULL, mm_allocator);
    } else {
        // exact WFA
        affine_wavefronts = affine_wavefronts_new_complete(
            target_length, query_length, affine_penalties, NULL, mm_allocator);
    }

    aln.j = j;
    aln.i = i;
    aln.query_length = query_length;
    aln.target_length = target_length;

    const int max_score = (target_length + query_length) * 2;

    aln.score = wfa::affine_wavefronts_align_bounded(
        affine_wavefronts,
        target+i,
        target_length,
        query+j,
        query_length,
        max_score);

    aln.ok = aln.score < max_score;
    if (aln.ok) {
        // correct X/M errors in the cigar
        hack_cigar(affine_wavefronts->edit_cigar, query, target, query_length, target_length, aln.j, aln.i);
#ifdef VALIDATE_WFA_WFLIGN
        if (!validate_cigar(affine_wavefronts->edit_cigar, query, target, query_length, target_length, aln.j, aln.i)) {
            std::cerr << "cigar failure at alignment " << aln.j << " " << aln.i << std::endl;
            unpack_display_cigar(affine_wavefronts->edit_cigar, query, target, query_length, target_length, aln.j, aln.i);
            std::cerr << ">query" << std::endl << std::string(query+j, query_length) << std::endl;
            std::cerr << ">target" << std::endl << std::string(target+i, target_length) << std::endl;
            assert(false);
        }
#endif
        wflign_edit_cigar_copy(&aln.edit_cigar, &affine_wavefronts->edit_cigar);
    }
    // cleanup wavefronts to keep memory low
    affine_wavefronts_delete(affine_wavefronts);

}

EdlibAlignResult do_edlib_patch_alignment(
    const char* query,
    const uint64_t& j,
    const uint64_t& query_length,
    const char* target,
    const uint64_t& i,
    const uint64_t& target_length,
    const EdlibAlignMode align_mode) {

    //std::cerr << "do_edlib_patch " << j << " " << query_length << " " << i << " " << target_length << std::endl;

    auto edlib_config = edlibNewAlignConfig(-1,
                                            align_mode,
                                            EDLIB_TASK_PATH,
                                            NULL, 0);

    return edlibAlign(query+j, query_length,
                      target+i, target_length,
                      edlib_config);

}

bool hack_cigar(
    wfa::edit_cigar_t& cigar,
    const char* query, const char* target,
    const uint64_t& query_aln_len, const uint64_t& target_aln_len,
    uint64_t j, uint64_t i) {
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    uint64_t j_max = j + query_aln_len;
    uint64_t i_max = i + target_aln_len;
    bool ok = true;
    //std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
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
                //std::cerr << "mismatch @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                cigar.operations[c] = 'X';
                ok = false;
            }
            if (j >= j_max) {
                //std::cerr << "query out of bounds @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                cigar.operations[c] = 'D';
                ok = false;
            }
            if (i >= i_max) {
                //std::cerr << "target out of bounds @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                cigar.operations[c] = 'I';
                ok = false;
            }
            ++j; ++i;
            break;
        case 'X':
            if (j < j_max && i < i_max && query[j] == target[i]) {
                cigar.operations[c] = 'M';
                ok = false;
            }
            ++j; ++i;
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

/*
wfa::edit_cigar_t filter_short_matches(
    const wfa::edit_cigar_t& cigar) {
// todo
}
*/

bool validate_cigar(
    const wfa::edit_cigar_t& cigar,
    const char* query, const char* target,
    const uint64_t& query_aln_len, const uint64_t& target_aln_len,
    uint64_t j, uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    uint64_t j_max = j + query_aln_len;
    uint64_t i_max = i + target_aln_len;
    bool ok = true;
    //std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (cigar.operations[c]) {
        case 'M':
            // check that we match
            if (query[j] != target[i]) {
                std::cerr << "mismatch @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            if (j >= j_max) {
                std::cerr << "query out of bounds @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            if (i >= i_max) {
                std::cerr << "target out of bounds @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            ++j; ++i;
            break;
        case 'X':
            ++j; ++i;
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

bool validate_trace(
    const std::vector<char>& tracev,
    const char* query, const char* target,
    const uint64_t& query_aln_len, const uint64_t& target_aln_len,
    uint64_t j, uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = 0;
    const int end_idx = tracev.size();
    uint64_t j_max = j + query_aln_len;
    uint64_t i_max = i + target_aln_len;
    bool ok = true;
    //std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (tracev[c]) {
        case 'M':
            // check that we match
            if (query[j] != target[i]) {
                std::cerr << "mismatch @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            if (j >= j_max) {
                std::cerr << "query out of bounds @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            if (i >= i_max) {
                std::cerr << "target out of bounds @ " << j << " " << i << " " << query[j] << " " << target[i] << std::endl;
                ok = false;
            }
            ++j; ++i;
            break;
        case 'X':
            ++j; ++i;
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

bool unpack_display_cigar(
    const wfa::edit_cigar_t& cigar,
    const char* query, const char* target,
    const uint64_t& query_aln_len, const uint64_t& target_aln_len,
    uint64_t j, uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    uint64_t j_max = j + query_aln_len;
    uint64_t i_max = i + target_aln_len;
    //std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (cigar.operations[c]) {
        case 'M':
            // check that we match
            std::cerr << "M" << " "
                      << j << " " << i << " "
                      << "q:" << query[j] << " "
                      << "t:" << target[i] << " "
                      << (query[j] == target[i] ? " " : "👾")
                      << (j >= j_max ? "⚡" : " ")
                      << (i >= i_max ? "🔥" : " ")
                      << std::endl;
            ++j; ++i;
            break;
        case 'X':
            std::cerr << "X" << " "
                      << j << " " << i << " "
                      << "q:" << query[j] << " "
                      << "t:" << target[i] << " "
                      << (query[j] != target[i] ? " " : "👾")
                      << (j >= j_max ? "⚡" : " ")
                      << (i >= i_max ? "🔥" : " ")
                      << std::endl;
            ++j; ++i;
            break;
        case 'I':
            std::cerr << "I" << " "
                      << j << " " << i << " "
                      << "q:" << query[j] << " "
                      << "t:" << "|" << "  "
                      << (j >= j_max ? "⚡" : " ")
                      << (i >= i_max ? "🔥" : " ")
                      << std::endl;
            ++j;
            break;
        case 'D':
            std::cerr << "D" << " "
                      << j << " " << i << " "
                      << "q:" << "|" << " "
                      << "t:" << target[i] << "  "
                      << (j >= j_max ? "⚡" : " ")
                      << (i >= i_max ? "🔥" : " ")
                      << std::endl;
            ++i;
            break;
        default:
            break;
        }
    }
    return true;
}

void write_merged_alignment(
    std::ostream& out,
    const std::vector<alignment_t*>& trace,
    wfa::mm_allocator_t* const mm_allocator,
    wfa::affine_penalties_t* const affine_penalties,
    const bool& emit_md_tag,
    const bool& paf_format_else_sam,
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
    const uint64_t& dropout_rescue_max,
    const long& elapsed_time_wflambda_ms,
    const bool& with_endline) {

    uint64_t target_length_mut = target_length;

    // patching parameters
    const uint64_t min_wfa_length = 16;
    const uint64_t min_edlib_length = 0;
    const int min_wf_length = 16;
    const int max_dist_threshold = 128;

    // we need to get the start position in the query and target
    // then run through the whole alignment building up the cigar
    // finally emitting it
    // our final cigar
    //
    //std::string cigarstr;
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

    //double mash_dist_sum = 0;
    uint64_t ok_alns = 0;

    const uint64_t erode_k = 13;

    auto start_time = std::chrono::steady_clock::now();

#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] processing traceback" << std::endl;
#endif
    // write trace into single cigar vector
    std::vector<char> tracev;
    {
        std::vector<char> erodev;
        {
            std::vector<char> rawv;
            // copy
#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] copying traceback" << std::endl;
#endif
            for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
                auto& aln = **x;
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
                        const auto& c = aln.edit_cigar.operations[i];
                        switch (c) {
                        case 'M': case 'X':
                            ++query_aligned_length; ++target_aligned_length; break;
                        case 'I': ++query_aligned_length; break;
                        case 'D': ++target_aligned_length; break;
                        default: break;
                        }
                        rawv.push_back(c);
                    }
                    query_end = aln.j + query_aligned_length;
                    target_end = aln.i + target_aligned_length;
                }
            }
#ifdef WFLIGN_DEBUG
            std::cerr << "[wflign::wflign_affine_wavefront] eroding traceback at k=" << erode_k << std::endl;
#endif
            //erode by removing matches < k
            for (uint64_t i = 0; i < rawv.size(); ) {
                if (rawv[i] == 'M') {
                    uint64_t j = i;
                    while (++j < rawv.size() && rawv[j] == 'M') { }
                    if (j-i < erode_k) {
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
        // normalize: sort so that I<D and otherwise leave it as-is
#ifdef WFLIGN_DEBUG
        std::cerr << "[wflign::wflign_affine_wavefront] normalizing eroded traceback" << std::endl;
#endif
        sort_indels(erodev);

#ifdef WFLIGN_DEBUG
        std::cerr << "[wflign::wflign_affine_wavefront] got normalized eroded traceback: ";
        for (auto c : erodev) {
            std::cerr << c;
        }
        std::cerr << std::endl;
#endif
        // patch: walk the cigar, patching directly when we have simultaneous gaps in query and ref
        // and adding our results to the final trace as we go
        auto q = erodev.begin();
        uint64_t query_pos = query_start;
        uint64_t target_pos = target_start;
        int64_t last_match_query = -1;
        int64_t last_match_target = -1;
        // get to the first match ... we'll not yet try to patch the alignment tips
        while (q != erodev.end()) {
            while (q != erodev.end() && (*q == 'M' || *q == 'X')) {
                if (*q == 'M') {
                    /*
                    std::cerr << "q: " << query[query_pos] << " "
                              << "t: " << target[target_pos] << std::endl;
                    */
                    if (query_pos >= query_length || target_pos >= target_length) {
                        std::cerr << "[wflign::wflign_affine_wavefront] corrupted traceback (out of bounds) for "
                                  << query_name << " " << query_offset << " "
                                  << target_name << " " << target_offset << std::endl;
                        exit(1);
                    }
                    if (query[query_pos] != target[target_pos]) {
                        std::cerr << "[wflign::wflign_affine_wavefront] corrupted traceback (match/mismatch confusion) for "
                                  << query_name << " " << query_offset << " "
                                  << target_name << " " << target_offset << std::endl;
                        exit(1);
                    }
                }
                tracev.push_back(*q);
                last_match_query = query_pos;
                last_match_target = target_pos;
                ++query_pos;
                ++target_pos;
                ++q;
            }
            // how long a gap?
            uint64_t query_delta = 0;
            while (q != erodev.end() && *q == 'I') {
                ++query_delta;
                ++q;
            }
            uint64_t target_delta = 0;
            while (q != erodev.end() && *q == 'D') {
                ++target_delta;
                ++q;
            }
            // how long was our last gap?
            // if it's long enough, patch it
            if (q != erodev.end()) {
                bool got_alignment = false;
                if (last_match_query > -1 && last_match_target > -1 &&
                    query_delta > 0 && target_delta > 0 &&
                    (query_delta < dropout_rescue_max || target_delta < dropout_rescue_max)) {

                    uint64_t patch_target_aligned_length = 0;
                    uint64_t patch_query_aligned_length = 0;
                    //int query_delta = aln.j - query_end;
                    //int target_delta = aln.i - target_end;
                    if (query_delta > min_wfa_length && target_delta > min_wfa_length) {
                        alignment_t patch_aln;
#ifdef WFLIGN_DEBUG
                        std::cerr << "do_wfa_patch_alignment" << std::endl;
#endif
#ifdef WFLIGN_SHOW_PATCH
                        std::cerr << "[wflign::wflign_affine_wavefront] patching in "
                                  << query_name << " " << query_offset << " @ " << query_pos
                                  << target_name << " " << target_offset << " @ " << target_pos
                                  << std::endl;
#endif
                        do_wfa_patch_alignment(
                            query, query_pos, query_delta,
                            target, target_pos, target_delta,
                            min_wf_length, max_dist_threshold,
                            mm_allocator, affine_penalties, patch_aln);
                        if (patch_aln.ok) {
                            //std::cerr << "got an ok patch aln" << std::endl;
                            got_alignment = true;
                            const int start_idx = patch_aln.edit_cigar.begin_offset;
                            const int end_idx = patch_aln.edit_cigar.end_offset;
                            for (int i = start_idx; i < end_idx; i++) {
                                tracev.push_back(patch_aln.edit_cigar.operations[i]);
                            }
                        }
                    } else if (query_delta > min_edlib_length && target_delta > min_edlib_length) {
#ifdef WFLIGN_DEBUG
                        std::cerr << "do_edlib_patch_alignment" << std::endl;
#endif
                        auto result = do_edlib_patch_alignment(
                            query, query_pos, query_delta,
                            target, target_pos, target_delta,
                            EDLIB_MODE_NW);
                        if (result.status == EDLIB_STATUS_OK
                            && result.alignmentLength != 0
                            && result.editDistance >= 0) {
                            got_alignment = true;
                            // copy it into the trace
                            char moveCodeToChar[] = {'M', 'I', 'D', 'X'};
                            auto& end_idx = result.alignmentLength;
                            for (int i = 0; i < end_idx; i++) {
                                tracev.push_back(moveCodeToChar[result.alignment[i]]);
                            }
                        }
                        edlibFreeAlignResult(result);
                    }
                }
                // add in stuff if we didn't align
                if (!got_alignment) {
                    for (uint64_t i = 0; i < query_delta; ++i) {
                        tracev.push_back('I');
                    }
                    for (uint64_t i = 0; i < target_delta; ++i) {
                        tracev.push_back('D');
                    }
                }
                //std::cerr << "query_delta " << query_delta << std::endl;
                //std::cerr << "target_delta " << target_delta << std::endl;
                query_pos += query_delta;
                target_pos += target_delta;
            } else {
                // we're at the end

                bool got_alignment = false;

                if (query_delta > 0) {
                    // there is a piece of query
                    auto target_delta_x = target_delta +
                            (target_pos + target_delta + query_delta < target_total_length ?
                                    query_delta :
                                    target_total_length - (target_pos + target_delta));

                    auto result = do_edlib_patch_alignment(
                            query, query_pos, query_delta,
                            target, target_pos, target_delta_x,
                            EDLIB_MODE_SHW);

                    if (result.status == EDLIB_STATUS_OK
                        && result.alignmentLength != 0
                        && result.editDistance >= 0) {
                        got_alignment = true;

                        if (target_pos + target_delta_x > target_length) {
                            target_length_mut = target_pos + target_delta_x;
                        }

                        target_delta = target_delta_x;

                        // copy it into the trace
                        char moveCodeToChar[] = {'M', 'I', 'D', 'X'};

                        for (int i = 0; i < *result.startLocations; ++i) {
                            tracev.push_back('D');
                        }

                        auto& end_idx = result.alignmentLength;
                        for (int i = 0; i < end_idx; ++i) {
                            tracev.push_back(moveCodeToChar[result.alignment[i]]);
                        }

                        for (int i = *result.startLocations + result.alignmentLength; i < target_delta; ++i) {
                            tracev.push_back('D');
                        }
                    }

                    edlibFreeAlignResult(result);
                }

                if(!got_alignment){
                    // add in our tail gap / softclip
                    for (uint64_t i = 0; i < query_delta; ++i) {
                        tracev.push_back('I');
                    }
                    for (uint64_t i = 0; i < target_delta; ++i) {
                        tracev.push_back('D');
                    }
                }
                query_pos += query_delta; // not used
                target_pos += target_delta;
            }
        }
#ifdef WFLIGN_DEBUG
        std::cerr << "[wflign::wflign_affine_wavefront] got unsorted patched traceback: ";
        for (auto c : tracev) {
            std::cerr << c;
        }
        std::cerr << std::endl;
#endif

        //std::cerr << "sorting the indels in tracev" << std::endl;
        // normalize the indels
        sort_indels(tracev);
    }

#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] got full patched traceback: ";
    for (auto c : tracev) {
        std::cerr << c;
    }
    std::cerr << std::endl;
#endif

#ifdef VALIDATE_WFA_WFLIGN
    if (!validate_trace(tracev, query, target, query_length, target_length_mut, query_start, target_start)) {
        std::cerr << "cigar failure at alignment "
                  << "\t" << query_total_length
                  << "\t" << query_offset + (query_is_rev ? query_length - query_end : query_start)
                  << "\t" << query_offset + (query_is_rev ? query_length - query_start : query_end)
                  << "\t" << (query_is_rev ? "-" : "+")
                  << "\t" << target_name
                  << "\t" << target_total_length
                  << "\t" << target_offset + target_start
                  << "\t" << target_offset + target_end << std::endl;
        exit(1);
    }
#endif

    // convert trace to cigar, get correct start and end coordinates
    char* cigarv = alignment_to_cigar(tracev,
                                      total_target_aligned_length,
                                      total_query_aligned_length,
                                      matches,
                                      mismatches,
                                      insertions,
                                      inserted_bp,
                                      deletions,
                                      deleted_bp);

    // gap-compressed identity
    const double gap_compressed_identity = (double)matches / (matches + mismatches + insertions + deletions);

    const uint64_t edit_distance = mismatches + inserted_bp + deleted_bp;

    // identity over the full block
    const double block_identity = (double)matches / (matches + edit_distance);

    auto write_tag_and_md_string = [&](std::ostream& out, const char* c, const int target_start) {
        out << "MD:Z:";

        char last_op = '\0';
        int last_len = 0;
        int t_off = target_start, l_MD = 0;
        int l = 0;
        int x = 0;
        while (c[x] != '\0') {
            while (isdigit(c[x])) ++x;
            char op = c[x];
            int len = 0;
            std::from_chars(c+l, c+x, len);
            l = ++x;
            if (last_len) {
                if (last_op == op) {
                    len += last_len;
                } else {
                    //std::cerr << t_off << "   " << last_len << last_op << std::endl;

                    if (last_op == '=') {
                        l_MD += last_len;
                        t_off += last_len;
                    }else if (last_op == 'X') {
                        for (uint64_t ii = 0; ii < last_len; ++ii) {
                            out << l_MD << target[t_off + ii];
                            l_MD = 0;
                        }

                        t_off += last_len;
                    }else if (last_op == 'D') {
                        out << l_MD << "^";
                        for (uint64_t ii = 0; ii < last_len; ++ii) {
                            out << target[t_off + ii];
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
            //std::cerr << t_off << "   " << last_len << last_op << std::endl;

            if (last_op == '=') {
                out << last_len + l_MD;
            }else if (last_op == 'X') {
                for (uint64_t ii = 0; ii < last_len; ++ii) {
                    out << l_MD << target[t_off + ii];
                    l_MD = 0;
                }
                out << "0";
            }else if (last_op == 'I'){
                out << l_MD;
            }else if (last_op == 'D') {
                out << l_MD << "^";
                for (uint64_t ii = 0; ii < last_len; ++ii) {
                    out << target[t_off + ii];
                }
                out << "0";
            }
        }
    };
    //std::cerr << "target_offset: " << target_offset << " -- target_start: " << target_start << std::endl;
    if (gap_compressed_identity >= min_identity) {
        long elapsed_time_patching_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count();

        const std::string timings = "wt:i:" + std::to_string(elapsed_time_wflambda_ms) + "\tpt:i:" + std::to_string(elapsed_time_patching_ms);

        if (paf_format_else_sam) {
            out << query_name
                << "\t" << query_total_length
                << "\t" << query_offset + (query_is_rev ? query_length - query_end : query_start)
                << "\t" << query_offset + (query_is_rev ? query_length - query_start : query_end)
                << "\t" << (query_is_rev ? "-" : "+")
                << "\t" << target_name
                << "\t" << target_total_length
                << "\t" << target_offset + target_start
                << "\t" << target_offset + target_end
                << "\t" << matches
                << "\t" << std::max(total_target_aligned_length, total_query_aligned_length)
                << "\t" << std::round(float2phred(1.0-block_identity))
                << "\t" << "as:i:" << total_score
                << "\t" << "gi:f:" << gap_compressed_identity
                << "\t" << "bi:f:" << block_identity
                //<< "\t" << "md:f:" << mash_dist_sum / trace.size()
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "ii:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "dd:i:" << deleted_bp
                << "\t" << "cg:Z:" << cigarv;

                if (emit_md_tag) {
                    out << "\t";

                    write_tag_and_md_string(out, cigarv, target_start);
                }

                out << "\t" << timings << "\n";
        } else {
            uint64_t query_start_pos = query_offset + (query_is_rev ? query_length - query_end : query_start);
            uint64_t query_end_pos = query_offset + (query_is_rev ? query_length - query_start : query_end);

            out << query_name                                                   // Query template NAME
                << "\t" << (query_is_rev ? "16" : "0")                          // bitwise FLAG
                << "\t" << target_name                                          // Reference sequence NAME
                << "\t" << target_offset + target_start + 1                     // 1-based leftmost mapping POSition
                << "\t" << std::round(float2phred(1.0-block_identity))          // MAPping Quality
                << "\t";

            // CIGAR
            if (query_is_rev) {
                if (query_length > query_end_pos) {
                    out << (query_length - query_end_pos) << "S";
                }
            } else {
                if (query_start_pos > 0) {
                    out << query_start_pos << "S";
                }
            }
            out << cigarv;
            if (query_is_rev) {
                if (query_start_pos > 0) {
                    out << query_start_pos << "S";
                }
            } else {
                if (query_length > query_end_pos) {
                    out << (query_length - query_end_pos) << "S";
                }
            }

            out << "\t" << "*"                                                  // Reference name of the mate/next read
                << "\t" << "0"                                                  // Position of the mate/next read
                << "\t" << "0"                                                  // observed Template LENgth
                << "\t";

            // segment SEQuence
            for (uint64_t p = 0; p < query_length; ++p) {
                out << query[p];
            }

            out << "\t" << "*"                                                  // ASCII of Phred-scaled base QUALity+33
                << "\t" << "NM:i:" << edit_distance
                << "\t" << "AS:i:" << total_score
                << "\t" << "gi:f:" << gap_compressed_identity
                << "\t" << "bi:f:" << block_identity
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

            out << "\t" << timings << "\n";
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
    const bool& with_endline) {
//    bool aligned = false;
    //Output to file
    //auto& result = aln.result;
    if (aln.ok) {

        /*
        if (result.status == EDLIB_STATUS_OK
        && result.alignmentLength != 0
        && result.editDistance >= 0) {
        */

        uint64_t matches = 0;
        uint64_t mismatches = 0;
        uint64_t insertions = 0;
        uint64_t inserted_bp = 0;
        uint64_t deletions = 0;
        uint64_t deleted_bp = 0;
        uint64_t refAlignedLength = 0;
        uint64_t qAlignedLength = 0;

        char* cigar = wfa_alignment_to_cigar(&aln.edit_cigar,
                                             refAlignedLength,
                                             qAlignedLength,
                                             matches,
                                             mismatches,
                                             insertions,
                                             inserted_bp,
                                             deletions,
                                             deleted_bp);

        size_t alignmentRefPos = aln.i;
        //double identity = (double)matches / (matches + mismatches + insertions + deletions);
        double gap_compressed_identity = (double)matches / (matches + mismatches + insertions + deletions);
        double block_identity = (double)matches / (matches + mismatches + inserted_bp + deleted_bp);
        // convert our coordinates to be relative to forward strand (matching PAF standard)

        if (gap_compressed_identity >= min_identity) {
            uint64_t q_start;
            if (query_is_rev) {
                q_start = query_offset + (query_length - (aln.j + qAlignedLength));
            } else {
                q_start = query_offset + aln.j;
            }
            out << query_name
                << "\t" << query_total_length
                << "\t" << q_start
                << "\t" << q_start + qAlignedLength
                << "\t" << (query_is_rev ? "-" : "+")
                << "\t" << target_name
                << "\t" << target_total_length
                << "\t" << target_offset + alignmentRefPos
                << "\t" << target_offset + alignmentRefPos + refAlignedLength
                << "\t" << matches
                << "\t" << std::max(refAlignedLength, qAlignedLength)
                << "\t" << std::round(float2phred(1.0-block_identity))
                << "\t" << "as:i:" << aln.score
                << "\t" << "gi:f:" << gap_compressed_identity
                << "\t" << "bi:f:" << block_identity
                //<< "\t" << "md:f:" << aln.mash_dist
                //<< "\t" << "ma:i:" << matches
                //<< "\t" << "mm:i:" << mismatches
                //<< "\t" << "ni:i:" << insertions
                //<< "\t" << "bi:i:" << inserted_bp
                //<< "\t" << "nd:i:" << deletions
                //<< "\t" << "bd:i:" << deleted_bp
                << "\t" << "cg:Z:" << cigar;
            if (with_endline) {
                out << std::endl;
            }
        }
        free(cigar);
    }
}

char* alignment_to_cigar(
    const std::vector<char>& edit_cigar,
    uint64_t& target_aligned_length,
    uint64_t& query_aligned_length,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& inserted_bp,
    uint64_t& deletions,
    uint64_t& deleted_bp) {

    // the edit cigar contains a character string of ops
    // here we compress them into the standard cigar representation

    auto* cigar = new std::vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    const int start_idx = 0;
    const int end_idx = edit_cigar.size();

    //std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
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
    cigar->push_back(0);  // Null character termination.

    char* cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}

char* wfa_alignment_to_cigar(
    const wfa::edit_cigar_t* const edit_cigar,
    uint64_t& target_aligned_length,
    uint64_t& query_aligned_length,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& inserted_bp,
    uint64_t& deletions,
    uint64_t& deleted_bp) {

    // the edit cigar contains a character string of ops
    // here we compress them into the standard cigar representation

    auto* cigar = new std::vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    const int start_idx = edit_cigar->begin_offset;
    const int end_idx = edit_cigar->end_offset;

    //std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int i = start_idx; i <= end_idx; i++) {
        // if new sequence of same moves started
        if (i == end_idx || (edit_cigar->operations[i] != lastMove && lastMove != 0)) {
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
    cigar->push_back(0);  // Null character termination.

    char* cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}

char* edlib_alignment_to_cigar(
    const unsigned char* const alignment,
    const int alignment_length,
    uint64_t& target_aligned_length,
    uint64_t& query_aligned_length,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& inserted_bp,
    uint64_t& deletions,
    uint64_t& deleted_bp) {

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'M', 'I', 'D', 'X'};

    std::vector<char>* cigar = new std::vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    auto& end_idx = alignment_length;
    for (int i = 0; i <= end_idx; i++) {
        if (i == end_idx || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
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
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < end_idx) {
                numOfSameMoves = 0;
            }
        }
        if (i < end_idx) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0);  // Null character termination.
    char* cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}

double float2phred(const double& prob) {
    if (prob == 1)
        return 255;  // guards against "-0"
    double p = -10 * (double) log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
}

void wflign_edit_cigar_copy(
    wfa::edit_cigar_t* const edit_cigar_dst,
    wfa::edit_cigar_t* const edit_cigar_src) {
    edit_cigar_dst->max_operations = edit_cigar_src->max_operations;
    edit_cigar_dst->begin_offset = 0;
    edit_cigar_dst->end_offset = edit_cigar_src->end_offset - edit_cigar_src->begin_offset;
    edit_cigar_dst->score = edit_cigar_src->score;
    // alloc our ops
    edit_cigar_dst->operations = (char*)malloc(edit_cigar_dst->end_offset);
    memcpy(edit_cigar_dst->operations,
           edit_cigar_src->operations+edit_cigar_src->begin_offset,
           edit_cigar_dst->end_offset);
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

/*
void edlib_to_wflign_edit_cigar_copy(
    wfa::edit_cigar_t* const edit_cigar_dst,
    char* const edlib_cigar_src,
    const uint64_t& edit_distance,
    const uint64_t& edlib_cigar_len) {
    //edit_cigar_dst->max_operations = edlib_cigar_src->max_operations; //hmmm
    edit_cigar_dst->begin_offset = 0;
    edit_cigar_dst->end_offset = edlib_cigar_len;
    edit_cigar_dst->score = edit_distance;
    edit_cigar_dst->operations = (char*)malloc(edit_cigar_dst->end_offset);
    memcpy(edit_cigar_dst->operations,
           edlib_cigar_src,
           edit_cigar_dst->end_offset);
}
*/

}

}
