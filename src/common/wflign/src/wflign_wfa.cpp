#include "wflign_wfa.hpp"

namespace wflign {

namespace wavefront {

void wflign_affine_wavefront(
    std::ostream& out,
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

    // set up our implicit matrix
    uint64_t steps_per_segment = 2;
    uint64_t step_size = segment_length / steps_per_segment;

    // Pattern & Text
    const int pattern_length = query_length / step_size;
    const int text_length = target_length / step_size;

    // use exact WFA locally
    const int wfa_min_wavefront_length = 0;
    const int wfa_max_distance_threshold = 0;

    // Allocate MM
    wflambda::mm_allocator_t* const wflambda_mm_allocator = wflambda::mm_allocator_new(BUFFER_SIZE_8M);
    // Set penalties
    wflambda::affine_penalties_t wflambda_affine_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening = 6,
        .gap_extension = 2,
    };
    // Init Affine wflambda
    wflambda::affine_wavefronts_t* affine_wavefronts;
    if (wflambda_min_wavefront_length || wflambda_max_distance_threshold) {
        affine_wavefronts = wflambda::affine_wavefronts_new_reduced(
            pattern_length, text_length, &wflambda_affine_penalties,
            wflambda_min_wavefront_length, wflambda_max_distance_threshold,
            NULL, wflambda_mm_allocator);
    } else {
        affine_wavefronts = wflambda::affine_wavefronts_new_complete(
            pattern_length, text_length, &wflambda_affine_penalties, NULL, wflambda_mm_allocator);
    }

    // save computed alignments in a pair-indexed patchmap
    whash::patchmap<uint64_t,alignment_t*> alignments;

    // allocate vectors to store our sketches
    //whash::patchmap<uint64_t, std::vector<rkmh::hash_t>*> query_sketches;
    //whash::patchmap<uint64_t, std::vector<rkmh::hash_t>*> target_sketches;
    std::vector<std::vector<rkmh::hash_t>*> query_sketches(pattern_length, nullptr);
    std::vector<std::vector<rkmh::hash_t>*> target_sketches(text_length, nullptr);

    //std::cerr << "v" << "\t" << "h" << "\t" << "score" << "\t" << "aligned" << std::endl;

    // setup affine WFA
    wfa::mm_allocator_t* const wfa_mm_allocator = wfa::mm_allocator_new(BUFFER_SIZE_8M);
    wfa::affine_penalties_t wfa_affine_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening = 6,
        .gap_extension = 1,
    };
    const uint64_t minhash_kmer_size = 19;
    int v_max = 0;
    int h_max = 0;

    auto extend_match =
        [&](const int& v,
            const int& h) {
            bool aligned = false;
            uint64_t k = encode_pair(v, h);
            if (v >= 0 && h >= 0
                && v < pattern_length
                && h < text_length) {
                auto f = alignments.find(k);
                if (f != alignments.end()) {
                    aligned = true;
                } else  {
                    uint64_t query_begin = (v < pattern_length-1 ? v * step_size
                                            : query_length - segment_length);
                    uint64_t target_begin = (h < text_length-1 ? h * step_size
                                             : target_length - segment_length);
                    alignment_t* aln = new alignment_t();
                    aligned =
                        do_alignment(
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
                            min_identity,
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
                            if (s != nullptr) {
                                delete s;
                                s = nullptr;
                            }
                        }
                    }
                    if (h > h_max) {
                        h_max = h;
                        if (h >= wflambda_max_distance_threshold) {
                            auto& s = target_sketches[h - wflambda_max_distance_threshold];
                            if (s != nullptr) {
                                delete s;
                                s = nullptr;
                            }
                        }
                    }
                }
            }
            return aligned;
        };

    // accumulate runs of matches in reverse order
    // then trim the cigars of successive mappings
    //
    std::vector<alignment_t*> trace;
    
    auto trace_match =
        [&](const int& v, const int& h) {
            uint64_t k = encode_pair(v, h);
            auto f = alignments.find(k);
            if (f != alignments.end()) {
                trace.push_back(f->second);
                return true;
            } else {
                return false;
            }
        };

    // Align
    wflambda::affine_wavefronts_align(
        affine_wavefronts,
        extend_match,
        trace_match,
        pattern_length,
        text_length);

#ifdef WFLIGN_DEBUG
    // get alignment score
    const int score = wflambda::edit_cigar_score_gap_affine(
        &affine_wavefronts->edit_cigar, &wflambda_affine_penalties);

    std::cerr << "[wflign::wflign_affine_wavefront] alignment score " << score << " for query: " << query_name << " target: " << target_name << std::endl;
#endif

    // clean up sketches
    for (auto& s : query_sketches) {
        if (s != nullptr) {
            delete s;
            s = nullptr;
        }
    }
    for (auto& s : target_sketches) {
        if (s != nullptr) {
            delete s;
            s = nullptr;
        }
    }

    // clean up our WFA allocator
    wfa::mm_allocator_delete(wfa_mm_allocator);

    // todo: implement alignment identifier based on hash of the input, params, and commit
    // annotate each PAF record with it and the full alignment score

    // Trim alignments that overlap in the query
    if (!trace.empty()) {
        //int last_i = 0;
        //int last_j = 0;
        for (auto & x : trace) {
            (*x).keep_query_length = segment_length;
        }
        for (auto x = trace.rbegin()+1; x != trace.rend(); ++x) {
            auto& curr = **x;
            auto& last = **(x-1);
            int ovlp_j = last.j + segment_length - curr.j;
            if (ovlp_j > 0) {
                int trim_last = ovlp_j / 2;
                int trim_curr = ovlp_j - trim_last;
                last.keep_query_length -= trim_last;
                curr.keep_query_length -= trim_curr;
                curr.skip_query_start += trim_curr;
            }
        }
        for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
            write_alignment(out, **x,
                            query_name, query_total_length, query_offset, query_length,
                            query_is_rev,
                            target_name, target_total_length, target_offset, target_length,
                            min_identity);
        }
    }

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
bool do_alignment(
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
    const float& min_identity,
    alignment_t& aln) {

    // first make the sketches if we haven't yet
    if (query_sketch == nullptr) {
        query_sketch = new std::vector<rkmh::hash_t>();
        *query_sketch = rkmh::hash_sequence(query+j, segment_length, minhash_kmer_size, segment_length/10);
    }
    if (target_sketch == nullptr) {
        target_sketch = new std::vector<rkmh::hash_t>();
        *target_sketch = rkmh::hash_sequence(target+i, segment_length, minhash_kmer_size, segment_length/10);
    }

    // first check if our mash dist is inbounds
    double mash_dist = rkmh::compare(*query_sketch, *target_sketch, minhash_kmer_size);
    //std::cerr << "mash_dist = " << mash_dist << std::endl;

    int max_score = segment_length * 0.8;

    // the mash distance generally underestimates the actual divergence
    // but when it's high we are almost certain that it's not a match
    if (mash_dist > 0.1) {
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

        int match_count = 0;
        for (int q = affine_wavefronts->edit_cigar.begin_offset;
             q != affine_wavefronts->edit_cigar.end_offset; ++q) {
            match_count += affine_wavefronts->edit_cigar.operations[q] == 'M';
        }

        //aln.identity = std::min(1.0, (double)match_count / ((double)segment_length / 2));
        aln.identity = (double)match_count / ((double)segment_length);
        //std::cerr << "identity is " << aln.identity << " and score " << aln.score << std::endl;

        aln.j = j;
        aln.i = i;
        aln.mash_dist = mash_dist;
        // copy our edit cigar if we aligned
        aln.ok = aln.score < max_score && aln.identity > min_identity;
        if (aln.ok) {
            wflign_edit_cigar_copy(&aln.edit_cigar, &affine_wavefronts->edit_cigar);
        }
        // cleanup wavefronts to keep memory low
        affine_wavefronts_delete(affine_wavefronts);

        return aln.ok;
    }
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
        uint64_t deletions = 0;
        uint64_t softclips = 0;
        uint64_t refAlignedLength = 0;
        uint64_t qAlignedLength = 0;
        int skipped_target_start = 0;
        int kept_target_length = 0;

        char* cigar = alignmentToCigar(&aln.edit_cigar,
                                       aln.skip_query_start,
                                       aln.keep_query_length,
                                       skipped_target_start,
                                       kept_target_length,
                                       refAlignedLength,
                                       qAlignedLength,
                                       matches,
                                       mismatches,
                                       insertions,
                                       deletions,
                                       softclips);

        size_t alignmentRefPos = aln.i;
        double total = refAlignedLength + (qAlignedLength - softclips);
        double identity = (double)(total - mismatches * 2 - insertions - deletions) / total;
        // convert our coordinates to be relative to forward strand (matching PAF standard)
        uint64_t q_start;
        if (query_is_rev) {
            q_start = query_offset + (query_length - (aln.j + aln.skip_query_start + qAlignedLength));
        } else {
            q_start = query_offset + aln.j + aln.skip_query_start;
        }
        if (identity >= min_identity) {
            out << query_name
                << "\t" << query_total_length
                << "\t" << q_start
                << "\t" << q_start + qAlignedLength
                << "\t" << (query_is_rev ? "-" : "+")
                << "\t" << target_name
                << "\t" << target_total_length
                << "\t" << target_offset + alignmentRefPos + skipped_target_start
                << "\t" << target_offset + alignmentRefPos + skipped_target_start + refAlignedLength
                << "\t" << matches
                << "\t" << std::max(refAlignedLength, qAlignedLength)
                << "\t" << std::round(float2phred(1.0-identity))
                << "\t" << "as:i:" << aln.score
                << "\t" << "id:f:" << identity
                << "\t" << "md:f:" << aln.mash_dist
                << "\t" << "ma:i:" << matches
                << "\t" << "mm:i:" << mismatches
                << "\t" << "ni:i:" << insertions
                << "\t" << "nd:i:" << deletions
                << "\t" << "ns:i:" << softclips
                << "\t" << "cg:Z:" << cigar;
            if (with_endline) {
                out << std::endl;
            }
        }
        free(cigar);
    }
}

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
    uint64_t& softclips) {

    // the edit cigar contains a character string of ops
    // here we compress them into the standard cigar representation

    std::vector<char>* cigar = new std::vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    int seen_query = 0;
    int start_idx = edit_cigar->begin_offset;

    while (start_idx < edit_cigar->end_offset
           && seen_query < skip_query_start) {
        switch (edit_cigar->operations[start_idx++]) {
        case 'M':
        case 'X':
            ++skipped_target_start;
            ++seen_query;
            break;
        case 'I':
            ++seen_query;
            break;
        case 'D':
            ++skipped_target_start;
            break;
        default:
            break;
        }
    }
    int end_idx = start_idx;
    seen_query = 0;
    while (end_idx < edit_cigar->end_offset
        && seen_query < keep_query_length) {
        switch (edit_cigar->operations[end_idx++]) {
        case 'M':
        case 'X':
            ++kept_target_length;
            ++seen_query;
            break;
        case 'I':
            ++seen_query;
            break;
        case 'D':
            ++kept_target_length;
            break;
        default:
            break;
        }
    }
    if (end_idx == start_idx) {
        end_idx = edit_cigar->end_offset;
    }
    for (int i = start_idx; i <= end_idx; i++) {
        // if new sequence of same moves started
        if (i == end_idx || (edit_cigar->operations[i] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
            case 'M':
                matches += numOfSameMoves;
                qAlignedLength += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            case 'X':
                mismatches += numOfSameMoves;
                qAlignedLength += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            case 'I':
                // assume that starting and ending insertions are softclips
                if (i == end_idx || cigar->empty()) {
                    softclips += numOfSameMoves;
                } else {
                    insertions += numOfSameMoves;
                }
                qAlignedLength += numOfSameMoves;
                break;
            case 'D':
                deletions += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
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

}

}
