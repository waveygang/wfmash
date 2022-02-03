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

#include "dna.hpp"
#include "rkmh.hpp"
#include "wflign_alignment.hpp"
#include "robin-hood-hashing/robin_hood.h"
#include "WFA/bindings/cpp/WFAligner.hpp"

/*
 * Configuration
 */
//#define WFLIGN_DEBUG true // for debugging messages
//#define VALIDATE_WFA_WFLIGN

// Namespaces
namespace wflign {
namespace wavefront {

class WFlign {
public:
	// WFlambda parameters
	uint16_t segment_length;
	float min_identity;
	int _minhash_kmer_size;
	int wfa_mismatch_score;
	int wfa_gap_opening_score;
	int wfa_gap_extension_score;
	int wflambda_min_wavefront_length; // with these set at 0 we do exact WFA for wflambda
	int wflambda_max_distance_threshold;
	float mashmap_estimated_identity;
	// WFlign parameters
	int wflign_mismatch_score;
	int wflign_gap_opening_score;
	int wflign_gap_extension_score;
	float wflign_max_mash_dist;
	uint64_t wflign_max_len_major;
	uint64_t wflign_max_len_minor;
	uint16_t erode_k;
	// Query
	std::string* query_name;
	char* query;
	uint64_t query_total_length;
	uint64_t query_offset;
	uint64_t query_length;
	bool query_is_rev;
	// Target
	std::string* target_name;
	char* target;
	uint64_t target_total_length;
	uint64_t target_offset;
	uint64_t target_length;
	// Output
	std::ostream* out;
	const bool emit_tsv;
	std::ostream* out_tsv;
	bool merge_alignments;
	bool emit_md_tag;
	bool paf_format_else_sam;
	bool no_seq_in_sam;
	// Setup
	WFlign(
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
		const uint16_t erode_k);
	// Set output configuration
	void set_output(
		std::ostream* const out,
		const bool emit_tsv,
		std::ostream* const out_tsv,
		const bool merge_alignments,
		const bool emit_md_tag,
		const bool paf_format_else_sam,
		const bool no_seq_in_sam);
	// WFling affine
	void wflign_affine_wavefront(
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
		const uint64_t target_length);
};

} // namespace wavefront
} // namespace wflign

#endif /* WFLIGN_HPP_ */
