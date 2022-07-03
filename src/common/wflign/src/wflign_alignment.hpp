#ifndef WFLIGN_ALIGNMENT_HPP_
#define WFLIGN_ALIGNMENT_HPP_

#include <vector>
#include <cstdint>
#include "WFA/bindings/cpp/WFAligner.hpp"

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
    int gap_opening;
    int gap_extension;
} wflign_penalties_t;

/*
 * Wflign Alignment
 */
class alignment_t {
public:
    int j;
    int i;
    uint16_t query_length;
    uint16_t target_length;
    bool ok = false;
    bool keep = false;
    wflign_cigar_t edit_cigar;
    // Setup
    alignment_t();
    ~alignment_t();
    // Utils
    bool validate(
            const char* query,
            const char* target);
    void trim_front(int query_trim);
    void trim_back(int query_trim);
    // Display
    void display(void);
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

#endif /* WFLIGN_ALIGNMENT_HPP_ */