#ifndef WFLIGN_ALIGNMENT_HPP_
#define WFLIGN_ALIGNMENT_HPP_

#include <vector>
#include <cstdint>
#include "WFA2-lib/bindings/cpp/WFAligner.hpp"

/*
 * Cigar
 */
typedef struct {
    char* cigar_ops;
    int begin_offset;
    int end_offset;
} wflign_cigar_t;

/*
 * Penalties
 */
typedef struct {
    int match;
    int mismatch;
    int gap_opening1;
    int gap_extension1;
    int gap_opening2;
    int gap_extension2;
} wflign_penalties_t;

/*
 * Wflign Alignment
 */
class alignment_t {
public:
    int j;
    int i;
    int query_length;
    int target_length;
    bool ok;
    bool keep;
    bool is_rev;
    wflign_cigar_t edit_cigar;
    // Default constructor
    alignment_t();
    // Destructor
    ~alignment_t();
    // Copy constructor
    alignment_t(const alignment_t& other);
    // Move constructor
    alignment_t(alignment_t&& other) noexcept;
    // Copy assignment operator
    alignment_t& operator=(const alignment_t& other);
    // Move assignment operator
    alignment_t& operator=(alignment_t&& other) noexcept;
    // Trim functions
    void trim_front(int query_trim);
    void trim_back(int query_trim);
    // query_begin, query_end, target_begin, target_end
    // Accessors
    int query_begin();
    int query_end();
    int target_begin();
    int target_end();
};


/*
 * Wflign Trace-Pos: Links a position in a traceback matrix to its edit
 */
class trace_pos_t {
public:
    int j = 0;
    int i = 0;
    wflign_cigar_t* edit_cigar = nullptr;
    int offset = 0;
    // Setup
    trace_pos_t(
            const int j,
            const int i,
            wflign_cigar_t* const edit_cigar,
            const int offset);
    trace_pos_t();
    // Accessors
    bool incr();
    bool decr();
    bool at_end();
    char curr();
    bool equal(trace_pos_t& other);
    bool assigned();
};
/*
 * Validate
 */
bool validate_cigar(
        const wflign_cigar_t& cigar,
        const char* query,
        const char* target,
        const uint64_t& query_aln_len,
        const uint64_t& target_aln_len,
        uint64_t j,
        uint64_t i);
bool validate_trace(
        std::vector<char>& tracev,
        const char* query,
        const char* target,
        const uint64_t& query_aln_len,
        const uint64_t& target_aln_len,
        uint64_t j,
        uint64_t i);
/*
 * Alignment-CIGAR Adaptors
 */
char* alignment_to_cigar(
        const std::vector<char>& edit_cigar,
        const uint64_t& start_idx,
        const uint64_t& end_idx,
        uint64_t& target_aligned_length,
        uint64_t& query_aligned_length,
        uint64_t& matches,
        uint64_t& mismatches,
        uint64_t& insertions,
        uint64_t& inserted_bp,
        uint64_t& deletions,
        uint64_t& deleted_bp);
char* wfa_alignment_to_cigar(
        const wflign_cigar_t* const edit_cigar,
        uint64_t& target_aligned_length,
        uint64_t& query_aligned_length,
        uint64_t& matches,
        uint64_t& mismatches,
        uint64_t& insertions,
        uint64_t& inserted_bp,
        uint64_t& deletions,
        uint64_t& deleted_bp);
/*
 * Utils
 */
bool unpack_display_cigar(
        const wflign_cigar_t& cigar,
        const char* query,
        const char* target,
        const uint64_t query_aln_len,
        const uint64_t target_aln_len,
        uint64_t j,
        uint64_t i);
void wflign_edit_cigar_copy(
        wfa::WFAligner& wf_aligner,
        wflign_cigar_t* const cigar_dst);

int calculate_alignment_score(const wflign_cigar_t& cigar, const wflign_penalties_t& penalties);

#endif /* WFLIGN_ALIGNMENT_HPP_ */
