#include "wflign_alignment.hpp"

#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <iterator>
#include <cassert>
#include <vector>

/*
 * Wflign Alignment
 */

// Default constructor
alignment_t::alignment_t()
    : j(0), i(0), query_length(0), target_length(0), score(std::numeric_limits<int>::max()), ok(false), keep(false), is_rev(false) {
    edit_cigar = {nullptr, 0, 0};
}

// Destructor
alignment_t::~alignment_t() {
        free(edit_cigar.cigar_ops);
    }

// Copy constructor
alignment_t::alignment_t(const alignment_t& other)
    : j(other.j), i(other.i), query_length(other.query_length), is_rev(other.is_rev),
      target_length(other.target_length), ok(other.ok), keep(other.keep) {
    if (other.edit_cigar.cigar_ops) {
        edit_cigar.cigar_ops = (char*)malloc((other.edit_cigar.end_offset - other.edit_cigar.begin_offset) * sizeof(char));
        memcpy(edit_cigar.cigar_ops, other.edit_cigar.cigar_ops + other.edit_cigar.begin_offset, 
               (other.edit_cigar.end_offset - other.edit_cigar.begin_offset) * sizeof(char));
    } else {
        edit_cigar.cigar_ops = nullptr;
    }
    edit_cigar.begin_offset = 0;
    edit_cigar.end_offset = other.edit_cigar.end_offset - other.edit_cigar.begin_offset;
}

// Move constructor
alignment_t::alignment_t(alignment_t&& other) noexcept
    : j(other.j), i(other.i), query_length(other.query_length), is_rev(other.is_rev),
      target_length(other.target_length), ok(other.ok), keep(other.keep),
      edit_cigar(other.edit_cigar) {
    other.edit_cigar = {nullptr, 0, 0};
}

// Copy assignment operator
alignment_t& alignment_t::operator=(const alignment_t& other) {
    if (this != &other) {
        j = other.j;
        i = other.i;
        query_length = other.query_length;
        target_length = other.target_length;
        ok = other.ok;
        keep = other.keep;
        is_rev = other.is_rev;

        free(edit_cigar.cigar_ops);
        if (other.edit_cigar.cigar_ops) {
            edit_cigar.cigar_ops = (char*)malloc((other.edit_cigar.end_offset - other.edit_cigar.begin_offset) * sizeof(char));
            memcpy(edit_cigar.cigar_ops, other.edit_cigar.cigar_ops + other.edit_cigar.begin_offset, 
                   (other.edit_cigar.end_offset - other.edit_cigar.begin_offset) * sizeof(char));
        } else {
            edit_cigar.cigar_ops = nullptr;
        }
        edit_cigar.begin_offset = 0;
        edit_cigar.end_offset = other.edit_cigar.end_offset - other.edit_cigar.begin_offset;
    }
    return *this;
}

// Move assignment operator
alignment_t& alignment_t::operator=(alignment_t&& other) noexcept {
    if (this != &other) {
        j = other.j;
        i = other.i;
        query_length = other.query_length;
        target_length = other.target_length;
        ok = other.ok;
        keep = other.keep;
        is_rev = other.is_rev;

        free(edit_cigar.cigar_ops);
        edit_cigar = other.edit_cigar;
        other.edit_cigar = {nullptr, 0, 0};
    }
    return *this;
}

int alignment_t::query_begin() const {
    return j;
}

int alignment_t::query_end() const {
    return j + query_length;
}

int alignment_t::target_begin() const {
    return i;
}

int alignment_t::target_end() const {
    return i + target_length;
}


//void alignment_t::display(void) {
//    std::cerr << j << " " << i << " " << query_length << " "
//              << target_length << " " << ok << std::endl;
//    for (int x = edit_cigar.begin_offset; x < edit_cigar.end_offset; ++x) {
//        std::cerr << edit_cigar.cigar_ops[x];
//    }
//    std::cerr << std::endl;
//}
//bool alignment_t::validate(
//        const char *query,
//        const char *target) {
//    return validate_cigar(edit_cigar,query,target,query_length,target_length,j,i);
//}
void alignment_t::trim_front(int query_trim) {
    // this kills the alignment
    if (query_trim >= query_length) {
        ok = false;
        return;
    }
    // increment j and i appropriately
    int trim_to_j = j + query_trim;
    int x = edit_cigar.begin_offset;
    while (x < edit_cigar.end_offset && j < trim_to_j) {
        switch (edit_cigar.cigar_ops[x++]) {
            case 'M':
            case 'X':
                --query_length;
                --target_length;
                ++j;
                ++i;
                break;
            case 'I':
                --query_length;
                ++j;
                break;
            case 'D':
                --target_length;
                ++i;
                break;
            default:
                break;
        }
        if (target_length <= 0 || query_length <= 0) {
            ok = false;
            return;
        }
    }
    while (x < edit_cigar.end_offset && edit_cigar.cigar_ops[x] == 'D') {
        ++x;
        --target_length;
        ++i;
    }
    if (x == edit_cigar.end_offset)
        ok = false;
    edit_cigar.begin_offset = x;
}
void alignment_t::trim_back(int query_trim) {
    if (query_trim >= query_length) {
        ok = false;
        return;
    }
    int x = edit_cigar.end_offset;
    int q = 0;
    while (x > edit_cigar.begin_offset && q < query_trim) {
        switch (edit_cigar.cigar_ops[--x]) {
            case 'M':
            case 'X':
                --query_length;
                --target_length;
                ++q;
                break;
            case 'I':
                --query_length;
                ++q;
                break;
            case 'D':
                --target_length;
                break;
            default:
                break;
        }
        if (target_length <= 0 || query_length <= 0) {
            ok = false;
            return;
        }
    }
    while (x >= edit_cigar.begin_offset &&
           edit_cigar.cigar_ops[x - 1] == 'D') {
        --x;
        --target_length;
    }
    if (x == edit_cigar.begin_offset) ok = false;
    edit_cigar.end_offset = x;
}
/*
 * Wflign Trace-Pos: Links a position in a traceback matrix to its edit
 */
trace_pos_t::trace_pos_t(
        const int j,
        const int i,
        wflign_cigar_t* const edit_cigar,
        const int offset) {
    this->j = j;
    this->i = i;
    this->edit_cigar = edit_cigar;
    this->offset = offset;
}
trace_pos_t::trace_pos_t() {
    this->j = 0;
    this->i = 0;
    this->edit_cigar = nullptr;
    this->offset = 0;
}
bool trace_pos_t::incr() {
    // FIXME ANDREA Shouldn't it be
    // if (offset < edit_cigar->end_offset-1)
    if (offset < edit_cigar->end_offset) {
        switch (curr()) {
            case 'M':
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
        ++offset;
        return true;
    } else {
        return false;
    }
}
bool trace_pos_t::decr() {
    if (offset > 0) {
        --offset;
        switch (curr()) {
            case 'M':
            case 'X':
                --j;
                --i;
                break;
            case 'I':
                --j;
                break;
            case 'D':
                --i;
                break;
            default:
                break;
        }
        return true;
    } else {
        return false;
    }
}
bool trace_pos_t::at_end() {
    return offset == edit_cigar->end_offset;
}
char trace_pos_t::curr() {
    assert(!at_end());
    return edit_cigar->cigar_ops[offset];
}
bool trace_pos_t::equal(trace_pos_t& other) {
    return j == other.j &&
           i == other.i &&
           curr() == 'M' &&
           curr() == other.curr();
}
bool trace_pos_t::assigned() {
    return edit_cigar != nullptr;
}
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
        uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    bool ok = true;
    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (cigar.cigar_ops[c]) {
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
bool validate_trace(
        std::vector<char>& tracev,
        const char* query,
        const char* target,
        const uint64_t& query_aln_len,
        const uint64_t& target_aln_len,
        uint64_t j,
        uint64_t i) {
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
        uint64_t& deleted_bp) {
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
char* wfa_alignment_to_cigar(
        const wflign_cigar_t* const edit_cigar,
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

    auto *cigar = new std::vector<char>();
    char lastMove = 0; // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    const int start_idx = edit_cigar->begin_offset;
    const int end_idx = edit_cigar->end_offset;

    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int i = start_idx; i <= end_idx; i++) {
        // if new sequence of same moves started
        if (i == end_idx ||
            (edit_cigar->cigar_ops[i] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
                case '=':
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
            lastMove = edit_cigar->cigar_ops[i];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0); // Null character termination.

    char *cigar_ = (char *)malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}
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
        uint64_t i) {
    // check that our cigar matches where it claims it does
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = start_idx; c < end_idx; c++) {
        // if new sequence of same moves started
        switch (cigar.cigar_ops[c]) {
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
void wflign_edit_cigar_copy(
        wfa::WFAligner& wf_aligner,
        wflign_cigar_t* const cigar_dst) {
    // Retrieve CIGAR
    char* cigar_ops;
    int cigar_length;
    wf_aligner.getAlignment(&cigar_ops,&cigar_length);
    // Allocate
    cigar_dst->cigar_ops = (char*)malloc(cigar_length);
    // Copy
    cigar_dst->begin_offset = 0;
    cigar_dst->end_offset = cigar_length;
    memcpy(cigar_dst->cigar_ops,cigar_ops,cigar_length);
}

int calculate_alignment_score(const wflign_cigar_t& cigar, const wflign_penalties_t& penalties) {
    int score = 0;
    char prev_op = '\0';
    int gap_length = 0;

    auto process_gap = [&](char op, int length) {
        switch (op) {
            case 'M':
                // Match is free (best case)
                break;
            case 'X':
                score += length * penalties.mismatch;
                break;
            case 'I':
            case 'D':
                score += penalties.gap_opening1 + penalties.gap_extension1;
                if (length > 1) {
                    score += std::min(
                        penalties.gap_extension1 * (length - 1),
                        penalties.gap_opening2 + penalties.gap_extension2 * (length - 1)
                    );
                }
                break;
        }
    };

    for (int i = cigar.begin_offset; i <= cigar.end_offset; ++i) {
        char op = (i < cigar.end_offset) ? cigar.cigar_ops[i] : '\0';  // '\0' to process last gap
        
        if (op != prev_op || i == cigar.end_offset) {
            if (gap_length > 0) {
                process_gap(prev_op, gap_length);
            }
            gap_length = 1;
        } else {
            ++gap_length;
        }
        
        prev_op = op;
    }

    return score;
}

std::string cigar_to_string(const wflign_cigar_t& cigar) {
    std::stringstream ss;
    const int start_idx = cigar.begin_offset;
    const int end_idx = cigar.end_offset;
    for (int c = start_idx; c < end_idx; c++) {
        ss << cigar.cigar_ops[c];
    }
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const alignment_t& aln) {
    return os << "Alignment: "
              << "Query(" << aln.query_begin() << "-" << aln.query_end() << "/" << aln.query_length << ") "
              << "Target(" << aln.target_begin() << "-" << aln.target_end() << "/" << aln.target_length << ") "
              << "Score=" << aln.score << " "
              << "Rev=" << (aln.is_rev ? "Yes" : "No") << " "
              << "Status=" << (aln.ok ? "OK" : "NotOK") << " "
              << "Keep=" << (aln.keep ? "Yes" : "No") << " "
              << "CIGAR=" << cigar_to_string(aln.edit_cigar) << " "
              << "Indices(i,j)=(" << aln.i << "," << aln.j << ")";
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
