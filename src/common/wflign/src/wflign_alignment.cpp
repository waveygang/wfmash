#include "wflign_alignment.hpp"

#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <iterator>
#include <vector>

/*
 * Wflign Alignment
 */
void alignment_t::display(void) {
	std::cerr << j << " " << i << " " << query_length << " "
			  << target_length << " " << ok << std::endl;
	for (int x = edit_cigar.begin_offset; x < edit_cigar.end_offset; ++x) {
		std::cerr << edit_cigar.operations[x];
	}
	std::cerr << std::endl;
}
bool alignment_t::validate(
		const char *query,
		const char *target) {
    return validate_cigar(edit_cigar, query, target, query_length,
                          target_length, j, i);
}
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
        switch (edit_cigar.operations[x++]) {
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
    while (x < edit_cigar.end_offset && edit_cigar.operations[x] == 'D') {
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
        switch (edit_cigar.operations[--x]) {
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
           edit_cigar.operations[x - 1] == 'D') {
        --x;
        --target_length;
    }
    if (x == edit_cigar.begin_offset)
        ok = false;
    edit_cigar.end_offset = x;
}
alignment_t::~alignment_t() {
	free(edit_cigar.operations);
}
/*
 * Wflign Trace-Pos: Links a position in a traceback matrix to its edit
 */
bool trace_pos_t::incr() {
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
    if (offset > edit_cigar->begin_offset) {
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
    return edit_cigar->operations[offset];
}
bool trace_pos_t::equal(const trace_pos_t &other) {
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
    const int cigar_length = cigar.cigar_length;
    const char* const cigar_ops = cigar.cigar_ops;
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    bool ok = true;
    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = 0; c < cigar_length; c++) {
        // if new sequence of same moves started
        switch (cigar_ops[c]) {
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
		const std::vector<char>& tracev,
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
		const std::vector<char> &edit_cigar,
		const uint64_t &start_idx,
		const uint64_t &end_idx,
		uint64_t &target_aligned_length,
		uint64_t &query_aligned_length,
		uint64_t &matches,
		uint64_t &mismatches,
		uint64_t &insertions,
		uint64_t &inserted_bp,
		uint64_t &deletions,
		uint64_t &deleted_bp) {

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
		uint64_t &target_aligned_length,
		uint64_t &query_aligned_length,
		uint64_t &matches,
		uint64_t &mismatches,
		uint64_t &insertions,
		uint64_t &inserted_bp,
		uint64_t &deletions,
		uint64_t &deleted_bp) {

    // the edit cigar contains a character string of ops
    // here we compress them into the standard cigar representation

    auto *cigar = new std::vector<char>();
    char lastMove = 0; // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    const int cigar_length = edit_cigar->cigar_length;
    const char* const cigar_ops = edit_cigar->cigar_ops;

    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int i = 0; i < cigar_length; i++) {
        // if new sequence of same moves started
        if (i == cigar_length || (cigar_ops[i] != lastMove && lastMove != 0)) {
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
            if (i < cigar_length) {
                numOfSameMoves = 0;
            }
        }
        if (i < cigar_length) {
            lastMove = cigar_ops[i];
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
    const int cigar_length = cigar.cigar_length;
    const char* const cigar_ops = cigar.cigar_ops;
    const uint64_t j_max = j + query_aln_len;
    const uint64_t i_max = i + target_aln_len;
    // std::cerr << "start to end " << start_idx << " " << end_idx << std::endl;
    for (int c = 0; c < cigar_length; c++) {
        // if new sequence of same moves started
        switch (cigar_ops[c]) {
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
