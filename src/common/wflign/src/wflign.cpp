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
 * Setup
 */
WFlign::WFlign(
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
























































} // namespace wavefront
} // namespace wflign
