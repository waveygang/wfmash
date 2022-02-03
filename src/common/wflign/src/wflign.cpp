#include <stddef.h>
#include <stdlib.h>
#include <cassert>
#include <chrono>
#include <iterator>
#include <string>

#include "src/wflign.hpp"

// Namespaces
namespace wflign {
namespace wavefront {

/*
 * Configuration
 */
#define MAX_LEN_FOR_PURE_WFA    50000
#define MIN_WF_LENGTH           256
#define MAX_DIST_THRESHOLD      4096

/*
 * Utils
 */
void wflign_edit_cigar_copy(
		wfa::WFAligner& wf_aligner,
		wflign_cigar_t* const cigar_dst) {
	// Retrieve WFA CIGAR
	// Retrieve CIGAR
	char* cigar_ops;
	int cigar_length;
	wf_aligner.getAlignmentCigar(&cigar_ops,&cigar_length);
	// Allocate
	cigar_dst->cigar_ops = (char*)malloc(cigar_length);
	// Copy
	cigar_dst->cigar_length = cigar_length;
	memcpy(cigar_dst->cigar_ops,cigar_ops,cigar_length);
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
		const uint16_t erode_k) {
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
	this->merge_alignments = false;
	this->emit_md_tag = false;
	this->paf_format_else_sam = false;
	this->no_seq_in_sam = false;
	// Internals
    num_alignments = 0;
    num_alignments_performed = 0;
}
/*
 * Output configuration
 */
void WFlign::set_output(
	std::ostream* const out,
	const bool emit_tsv,
	std::ostream* const out_tsv,
	const bool merge_alignments,
	const bool emit_md_tag,
	const bool paf_format_else_sam,
	const bool no_seq_in_sam) {
	this->out = out;
	this->emit_tsv = emit_tsv;
	this->out_tsv = out_tsv;
	this->merge_alignments = merge_alignments;
	this->emit_md_tag = emit_md_tag;
	this->paf_format_else_sam = paf_format_else_sam;
	this->no_seq_in_sam = no_seq_in_sam;
}
/*
 * WFling align
 */
void WFlign::wflign_affine_wavefront(
		std::string* const query_name,
		char* const query,
		const uint64_t query_total_length,
		const uint64_t query_offset,
		const uint64_t query_length,
		const bool query_is_rev,
		std::string* const target_name,
		char* const target,
		const uint64_t target_total_length,
		const uint64_t target_offset,
		const uint64_t target_length) {
	// Set query
	this->query_name = query_name;
	this->query = query;
	this->query_total_length = query_total_length;
	this->query_offset = query_offset;
	this->query_length = query_length;
	this->query_is_rev = query_is_rev;
	// Set target
	this->target_name = target_name;
	this->target = target;
	this->target_total_length = target_total_length;
	this->target_offset = target_offset;
	this->target_length = target_length;

    if (query_offset + query_length > query_total_length ||
        target_offset + target_length > target_total_length) {
        return;
    }

    auto minhash_kmer_size = _minhash_kmer_size;

    // Set penalties
    wflign_penalties_t wfa_affine_penalties;
    if (wfa_mismatch_score > 0 && wfa_gap_opening_score > 0 && wfa_gap_extension_score > 0){
    	wfa_affine_penalties.match = 0;
    	wfa_affine_penalties.mismatch = wfa_mismatch_score;
    	wfa_affine_penalties.gap_opening = wfa_gap_opening_score;
    	wfa_affine_penalties.gap_extension = wfa_gap_extension_score;
        minhash_kmer_size = 17;
    } else {
        if (mashmap_estimated_identity >= 0.99999) {
        	wfa_affine_penalties.match = 0;
        	wfa_affine_penalties.mismatch = 15;
        	wfa_affine_penalties.gap_opening = 25;
        	wfa_affine_penalties.gap_extension = 1;
            minhash_kmer_size = 17;
        } else if (mashmap_estimated_identity >= 0.97) {
        	wfa_affine_penalties.match = 0;
        	wfa_affine_penalties.mismatch = 9;
        	wfa_affine_penalties.gap_opening = 13;
        	wfa_affine_penalties.gap_extension = 1;
            minhash_kmer_size = 17;
        } else if (mashmap_estimated_identity >= 0.9) {
        	wfa_affine_penalties.match = 0;
        	wfa_affine_penalties.mismatch = 7;
        	wfa_affine_penalties.gap_opening = 11;
        	wfa_affine_penalties.gap_extension = 1;
            minhash_kmer_size = 16;
        } else if (mashmap_estimated_identity >= 0.8) {
        	wfa_affine_penalties.match = 0;
        	wfa_affine_penalties.mismatch = 3;
        	wfa_affine_penalties.gap_opening = 5;
        	wfa_affine_penalties.gap_extension = 1;
            minhash_kmer_size = 15;
        } else {
        	wfa_affine_penalties.match = 0;
        	wfa_affine_penalties.mismatch = 2;
        	wfa_affine_penalties.gap_opening = 4;
        	wfa_affine_penalties.gap_extension = 1;
            minhash_kmer_size = 13;
        }
    }

    // accumulate runs of matches in reverse order
    // then trim the cigars of successive mappings
    std::vector<alignment_t *> trace;

    uint64_t num_alignments = 0;
    uint64_t num_alignments_performed = 0;
    const auto start_time = std::chrono::steady_clock::now();

    if (query_length <= MAX_LEN_FOR_PURE_WFA && target_length <= MAX_LEN_FOR_PURE_WFA) {

    	wfa::WFAlignerGapAffine& wf_aligner =
    			new wfa::WFAlignerGapAffine(
    					wfa_affine_penalties.mismatch,
    					wfa_affine_penalties.gap_opening,
    					wfa_affine_penalties.gap_extension,
    					false,
    					wfa::WFAligner::WavefrontMemoryHigh);
    	wf_aligner.setReductionAdaptive(MIN_WF_LENGTH,MAX_DIST_THRESHOLD);
    	const int status = wf_aligner.alignEnd2End(target,target_length,query,query_length);

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

            wflign_edit_cigar_copy(wf_aligner,&aln->edit_cigar);

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
                MAX_LEN_FOR_PURE_WFA, min_identity,
                elapsed_time_wflambda_ms, num_alignments,
                num_alignments_performed, mashmap_estimated_identity,
                wflign_max_len_major, wflign_max_len_minor,
                erode_k,
                MIN_WF_LENGTH, MAX_DIST_THRESHOLD);

        // Free
        delete wf_aligner;
    } else {
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

        wflign_penalties_t wflambda_affine_penalties;
        if (wflign_mismatch_score > 0 && wflign_gap_opening_score > 0 && wflign_gap_extension_score > 0){
        	wflambda_affine_penalties.match = 0;
        	wflambda_affine_penalties.mismatch = wflign_mismatch_score;
        	wflambda_affine_penalties.gap_opening = wflign_gap_opening_score;
        	wflambda_affine_penalties.gap_extension = wflign_gap_extension_score;
        } else {
            if (mashmap_estimated_identity >= 0.99999) {
            	wflambda_affine_penalties.match = 0;
            	wflambda_affine_penalties.mismatch = 17;
            	wflambda_affine_penalties.gap_opening = 27;
            	wflambda_affine_penalties.gap_extension = 1;
            } else if (mashmap_estimated_identity >= 0.97) {
            	wflambda_affine_penalties.match = 0;
            	wflambda_affine_penalties.mismatch = 13;
            	wflambda_affine_penalties.gap_opening = 21;
            	wflambda_affine_penalties.gap_extension = 1;
            } else if (mashmap_estimated_identity >= 0.9) {
            	wflambda_affine_penalties.match = 0;
            	wflambda_affine_penalties.mismatch = 9;
            	wflambda_affine_penalties.gap_opening = 14;
            	wflambda_affine_penalties.gap_extension = 1;
            } else if (mashmap_estimated_identity >= 0.8) {
            	wflambda_affine_penalties.match = 0;
            	wflambda_affine_penalties.mismatch = 7;
            	wflambda_affine_penalties.gap_opening = 11;
            	wflambda_affine_penalties.gap_extension = 1;
            } else {
            	wflambda_affine_penalties.match = 0;
            	wflambda_affine_penalties.mismatch = 4;
            	wflambda_affine_penalties.gap_opening = 6;
            	wflambda_affine_penalties.gap_extension = 1;
            }
        }

        float max_mash_dist_to_evaluate;
        if (wflign_max_mash_dist > 0) {
            max_mash_dist_to_evaluate = wflign_max_mash_dist;
        } else {
            // heuristic bound on the max mash dist, adaptive based on estimated
            // identity the goal here is to sparsify the set of alignments in the
            // wflambda layer we then patch up the gaps between them
            if (mashmap_estimated_identity >= 0.97) {
                max_mash_dist_to_evaluate = 0.05;
            } else if (mashmap_estimated_identity >= 0.9) {
                max_mash_dist_to_evaluate = 0.2;
            } else if (mashmap_estimated_identity >= 0.8) {
                max_mash_dist_to_evaluate = 0.7;
            } else if (mashmap_estimated_identity >= 0.7) {
                max_mash_dist_to_evaluate = 0.9;
            }
        }

        //std::cerr << "wfa_affine_penalties.mismatch " << wfa_affine_penalties.mismatch << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_opening " << wfa_affine_penalties.gap_opening << std::endl;
        //std::cerr << "wfa_affine_penalties.gap_extension " << wfa_affine_penalties.gap_extension << std::endl;
        //std::cerr << "wflambda_affine_penalties.mismatch " << wflambda_affine_penalties.mismatch << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_opening " << wflambda_affine_penalties.gap_opening << std::endl;
        //std::cerr << "wflambda_affine_penalties.gap_extension " << wflambda_affine_penalties.gap_extension << std::endl;
        //std::cerr << "max_mash_dist_to_evaluate " << max_mash_dist_to_evaluate << std::endl;

        // Configure the attributes of the wflambda-aligner
    	wfa::WFAlignerGapAffine& wflambda_aligner =
    			new wfa::WFAlignerGapAffine(
    					wfa_affine_penalties.mismatch,
    					wfa_affine_penalties.gap_opening,
    					wfa_affine_penalties.gap_extension,
    					false,
    					wfa::WFAligner::WavefrontMemoryHigh);
        if (wflambda_min_wavefront_length || wflambda_max_distance_threshold) {
        	wflambda_aligner.setReductionAdaptive(wflambda_min_wavefront_length,wflambda_max_distance_threshold);
        } else {
        	wflambda_aligner.setReductionNone();
        }

        wflambda_aligner.setMatchFunct(
              int (*matchFunct)(int,int,void*),
              void* matchFunctArguments);


        // save computed alignments in a pair-indexed map
        robin_hood::unordered_flat_map<uint64_t, alignment_t *> alignments;

        // allocate vectors to store our sketches
        std::vector<std::vector<rkmh::hash_t> *> query_sketches(pattern_length,
                                                                nullptr);
        std::vector<std::vector<rkmh::hash_t> *> target_sketches(text_length,
                                                                 nullptr);

    	wfa::WFAlignerGapAffine& wf_aligner =
    			new wfa::WFAlignerGapAffine(
    					wfa_affine_penalties.mismatch,
    					wfa_affine_penalties.gap_opening,
    					wfa_affine_penalties.gap_extension,
    					false,
    					wfa::WFAligner::WavefrontMemoryFull);
    	wflambda_aligner.setReductionNone();

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
                            wfa_max_distance_threshold, max_mash_dist_to_evaluate, mashmap_estimated_identity,
                            wf_aligner, &wfa_affine_penalties, *aln);
                    if (emit_tsv) {
                        // 0) Mis-match, alignment skipped
                        // 1) Mis-match, alignment performed
                        // 2) Match, alignment performed
                        out_tsv << v << "\t" << h << "\t" << (alignment_performed ? (aln->ok ? 2 : 1) : 0) << std::endl;
                    }
                    ++num_alignments;
                    if (alignment_performed) {
                        ++num_alignments_performed;
                        if (aln->ok){
                            is_a_match = true;
                            alignments[k] = aln;
                        } else {
                            alignments[k] = nullptr;
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
                        segment_length_to_use, min_identity,
                        elapsed_time_wflambda_ms, num_alignments,
                        num_alignments_performed, mashmap_estimated_identity,
                        wflign_max_len_major, wflign_max_len_minor,
                        erode_k,
                        MIN_WF_LENGTH, MAX_DIST_THRESHOLD);
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
        delete wf_aligner;
    }
}






















































} // namespace wavefront
} // namespace wflign
