#include "wflign_swizzle.hpp"
#include <iostream>

namespace wflign {

// Implementation copied from cigar_swap.cpp and wrapped in namespace
static std::string merge_cigar_ops(const std::string &cigar) {
    std::string merged;
    int current_count = 0;
    char current_op = '\0';

    size_t i = 0;
    while (i < cigar.size()) {
        int val = 0;
        while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
            val = val * 10 + (cigar[i] - '0');
            i++;
        }
        if (i >= cigar.size()) break;
        char op = cigar[i++];
        if (op == current_op) {
            current_count += val;
        } else {
            if (current_op != '\0') {
                merged += std::to_string(current_count);
                merged += current_op;
            }
            current_op = op;
            current_count = val;
        }
    }
    if (current_op != '\0') {
        merged += std::to_string(current_count);
        merged += current_op;
    }
    return merged;
}

static bool sequences_match(
    const std::string &query_seq,
    const std::string &target_seq,
    int64_t query_start,
    int64_t target_start,
    int N,
    bool debug) {
    if (query_start < 0 || target_start < 0) {
        if (debug) {
            std::cerr << "[swizzle-debug] sequences_match failed: negative start position"
                      << " query_start=" << query_start 
                      << " target_start=" << target_start << std::endl;
        }
        return false;
    }
    if (query_start + N > (int64_t)query_seq.size()) {
        if (debug) {
            std::cerr << "[swizzle-debug] sequences_match failed: query bounds"
                      << " query_start=" << query_start
                      << " N=" << N 
                      << " query_size=" << query_seq.size() << std::endl;
        }
        return false;
    }
    if (target_start + N > (int64_t)target_seq.size()) {
        if (debug) {
            std::cerr << "[swizzle-debug] sequences_match failed: target bounds"
                      << " target_start=" << target_start
                      << " N=" << N
                      << " target_size=" << target_seq.size() << std::endl;
        }
        return false;
    }
    for (int i = 0; i < N; i++) {
        if (query_seq[query_start + i] != target_seq[target_start + i]) {
            if (debug) {
                std::cerr << "[swizzle-debug] sequences_match failed: mismatch at position " << i
                          << " query=" << query_seq[query_start + i]
                          << " target=" << target_seq[target_start + i] << std::endl;
            }
            return false;
        }
    }
    return true;
}

static bool verify_cigar_alignment(
    const std::string &cigar,
    const std::string &query_seq,
    const std::string &target_seq,
    int64_t query_start,
    int64_t target_start,
    bool debug) {
    if (debug) {
        std::cerr << "[swizzle-debug] verify_cigar_alignment:"
                  << " cigar=" << cigar
                  << " query_start=" << query_start
                  << " target_start=" << target_start
                  << " query_size=" << query_seq.size()
                  << " target_size=" << target_seq.size() << std::endl;
    }

    int64_t qPos = query_start;
    int64_t tPos = target_start;

    size_t i = 0;
    while (i < cigar.size()) {
        int val = 0;
        while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
            val = val * 10 + (cigar[i] - '0');
            i++;
        }
        if (i >= cigar.size()) break;
        char op = cigar[i++];
        
        if (op == '=') {
            if (qPos < 0 || tPos < 0 ||
                qPos + val > (int64_t)query_seq.size() ||
                tPos + val > (int64_t)target_seq.size()) {
                if (debug) {
                    std::cerr << "[swizzle-debug] verify_cigar_alignment failed: bounds check"
                              << " op=" << op
                              << " val=" << val
                              << " qPos=" << qPos
                              << " tPos=" << tPos << std::endl;
                }
                return false;
            }
            for (int k = 0; k < val; k++) {
                if (query_seq[qPos + k] != target_seq[tPos + k]) {
                    if (debug) {
                        std::cerr << "[swizzle-debug] verify_cigar_alignment failed: mismatch"
                                  << " pos=" << k
                                  << " query=" << query_seq[qPos + k]
                                  << " target=" << target_seq[tPos + k] << std::endl;
                    }
                    return false;
                }
            }
            qPos += val;
            tPos += val;
        } else if (op == 'D') {
            if (tPos + val > (int64_t)target_seq.size()) {
                if (debug) {
                    std::cerr << "[swizzle-debug] verify_cigar_alignment failed: deletion bounds"
                              << " val=" << val
                              << " tPos=" << tPos
                              << " target_size=" << target_seq.size() << std::endl;
                }
                return false;
            }
            tPos += val;
        } else {
            if (debug) {
                std::cerr << "[swizzle-debug] verify_cigar_alignment failed: invalid op "
                          << op << std::endl;
            }
            return false;
        }
    }
    return true;
}

static bool extract_first_two_ops(
    const std::string &cigar,
    int &count1, char &op1,
    int &count2, char &op2,
    size_t &second_op_start,
    size_t &second_op_end) {
    size_t i = 0;
    {
        int val = 0;
        while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
            val = val * 10 + (cigar[i] - '0');
            i++;
        }
        if (i >= cigar.size()) return false;
        op1 = cigar[i++];
        count1 = val;
    }
    second_op_start = i;
    {
        int val = 0;
        while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
            val = val * 10 + (cigar[i] - '0');
            i++;
        }
        if (i >= cigar.size()) return false;
        op2 = cigar[i++];
        count2 = val;
        second_op_end = i;
    }
    return true;
}

static bool extract_last_two_ops(
    const std::string &cigar,
    int &count2, char &op2,
    int &count1, char &op1,
    size_t &second_last_op_start,
    size_t &last_op_start) {
    if (cigar.empty()) return false;

    size_t pos = cigar.size();
    while (pos > 0 && std::isdigit(static_cast<unsigned char>(cigar[pos - 1]))) {
        pos--;
    }
    if (pos == 0) return false;
    op2 = cigar[pos - 1];
    {
        size_t digit_end = pos - 1;
        int val = 0;
        int factor = 1;
        size_t j = digit_end;
        while (j > 0 && std::isdigit(static_cast<unsigned char>(cigar[j - 1]))) {
            val += (cigar[j - 1] - '0') * factor;
            factor *= 10;
            j--;
        }
        count2 = val;
        last_op_start = j;
    }

    if (last_op_start == 0) return false;
    {
        size_t pos2 = last_op_start;
        while (pos2 > 0 && std::isdigit(static_cast<unsigned char>(cigar[pos2 - 1]))) {
            pos2--;
        }
        if (pos2 == 0) return false;
        op1 = cigar[pos2 - 1];
        {
            size_t digit_end = pos2 - 1;
            int val = 0;
            int factor = 1;
            size_t j = digit_end;
            while (j > 0 && std::isdigit(static_cast<unsigned char>(cigar[j - 1]))) {
                val += (cigar[j - 1] - '0') * factor;
                factor *= 10;
                j--;
            }
            count1 = val;
            second_last_op_start = j;
        }
    }
    return true;
}

std::pair<int64_t,int64_t> alignment_end_coords(
    const std::string &cigar,
    int64_t query_start,
    int64_t target_start) {
    int64_t qPos = query_start;
    int64_t tPos = target_start;
    size_t i = 0;
    while (i < cigar.size()) {
        int val = 0;
        while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
            val = val * 10 + (cigar[i] - '0');
            i++;
        }
        if (i >= cigar.size()) break;
        char op = cigar[i++];
        if (op == '=') {
            qPos += val;
            tPos += val;
        } else if (op == 'D') {
            tPos += val;
        }
    }
    return {qPos, tPos};
}

std::string try_swap_start_pattern(
    const std::string &cigar,
    const std::string &query_seq,
    const std::string &target_seq,
    int64_t query_start,
    int64_t target_start) {
    if (!verify_cigar_alignment(cigar, query_seq, target_seq, query_start, target_start)) {
        return cigar;
    }

    int N, Dlen;
    char op1, op2;
    size_t second_op_start, second_op_end;
    if (!extract_first_two_ops(cigar, N, op1, Dlen, op2, second_op_start, second_op_end)) {
        return cigar;
    }

    if (op1 == '=' && op2 == 'D') {
        // Always enable debug logging
        const bool debug = true;
        std::cerr << "\n[swizzle-debug] Starting start pattern swizzle attempt..." << std::endl;
        std::cerr << "[swizzle-debug] Found pattern: " << N << op1 << " followed by " << Dlen << op2 << std::endl;

        if (debug) {
            std::cerr << "[swizzle-debug] try_swap_start_pattern:"
                      << " cigar=" << cigar
                      << " N=" << N
                      << " Dlen=" << Dlen
                      << " query_start=" << query_start
                      << " target_start=" << target_start
                      << " query_size=" << query_seq.size()
                      << " target_size=" << target_seq.size() << std::endl;
            
            // Print just the relevant portions (100 chars) around where we're trying to match
            std::cerr << "[swizzle-debug] query region ["
                      << query_start << ".." << (query_start + std::min(N, 100)) << "]: "
                      << query_seq.substr(query_start, std::min(N, 100)) << std::endl;
            std::cerr << "[swizzle-debug] target region ["
                      << (target_start + Dlen) << ".." << (target_start + Dlen + std::min(N, 100)) << "]: "
                      << target_seq.substr(target_start + Dlen, std::min(N, 100)) << std::endl;
        }

        // Check if we can match the sequences after shifting by Dlen
        if (sequences_match(query_seq, target_seq, query_start, target_start + Dlen, N, debug)) {
            std::string remainder = cigar.substr(second_op_end);
            std::string swapped = std::to_string(Dlen) + "D" +
                                std::to_string(N) + "=" +
                                remainder;
            swapped = merge_cigar_ops(swapped);

            // Verify the entire new alignment
            if (!verify_cigar_alignment(swapped, query_seq, target_seq, query_start, target_start, debug)) {
                if (debug) {
                    std::cerr << "[swizzle-debug] try_swap_start_pattern: verification failed for swapped CIGAR"
                              << " swapped=" << swapped << std::endl;
                }
                return cigar;
            }
            return swapped;
        }
    }
    return cigar;
}

std::string try_swap_end_pattern(
    const std::string &cigar,
    const std::string &query_seq,
    const std::string &target_seq,
    int64_t query_start,
    int64_t target_start) {
    if (!verify_cigar_alignment(cigar, query_seq, target_seq, query_start, target_start)) {
        return cigar;
    }

    int count2, count1;
    char op2, op1;
    size_t second_last_op_start, last_op_start;
    if (!extract_last_two_ops(cigar, count2, op2, count1, op1, second_last_op_start, last_op_start)) {
        return cigar;
    }

    if (op1 == 'D' && op2 == '=') {
        std::cerr << "\n[swizzle-debug] Starting end pattern swizzle attempt..." << std::endl;
        std::cerr << "[swizzle-debug] Found pattern: " << count1 << op1 << " followed by " << count2 << op2 << std::endl;
        int N = count2;
        int Dlen = count1;
        auto endCoords = alignment_end_coords(cigar, query_start, target_start);
        int64_t endQ = endCoords.first;
        int64_t endT = endCoords.second;

        if (sequences_match(query_seq, target_seq, endQ - N, endT - N - Dlen, N)) {
            std::string before = cigar.substr(0, second_last_op_start);
            std::string swapped = before +
                                std::to_string(N) + "=" +
                                std::to_string(Dlen) + "D";
            swapped = merge_cigar_ops(swapped);
            if (!verify_cigar_alignment(swapped, query_seq, target_seq, query_start, target_start)) {
                return cigar;
            }
            return swapped;
        }
    }
    return cigar;
}

static std::string drop_leading_trailing_deletions_if_all_eq(
    const std::string &cigar) {
    std::vector<std::pair<int,int>> ops;
    {
        size_t i = 0;
        while (i < cigar.size()) {
            int val = 0;
            while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
                val = val * 10 + (cigar[i] - '0');
                i++;
            }
            if (i >= cigar.size()) break;
            char op = cigar[i++];
            ops.push_back(std::make_pair(val, (int)op));
        }
    }
    if (ops.empty()) return cigar;

    for (auto &p : ops) {
        char op = (char)p.second;
        if (op != 'D' && op != '=') {
            return cigar;
        }
    }

    bool found_eq = false;
    for (size_t i = 0; i < ops.size(); i++) {
        char op = (char)ops[i].second;
        if (op == '=') {
            found_eq = true;
            break;
        }
    }
    if (!found_eq) {
        return cigar;
    }

    char firstOp = (char)ops[0].second;
    char lastOp  = (char)ops[ops.size()-1].second;

    if (firstOp != 'D' || lastOp != 'D') {
        return cigar;
    }

    ops.erase(ops.begin());
    ops.pop_back();

    std::string out;
    for (size_t i = 0; i < ops.size(); i++) {
        out += std::to_string(ops[i].first);
        out += (char)ops[i].second;
    }

    out = merge_cigar_ops(out);
    return out;
}

// Convert edit_cigar_t to string CIGAR format
std::string wfa_edit_cigar_to_string(const wflign_cigar_t& edit_cigar) {
    std::string cigar;
    int count = 0;
    char last_op = '\0';
    
    for (int i = edit_cigar.begin_offset; i < edit_cigar.end_offset; ++i) {
        char op = edit_cigar.cigar_ops[i];
        if (op == 'M') {
            op = '=';  // Convert any M to = since these are exact matches
        }
        if (op == last_op) {
            count++;
        } else {
            if (last_op != '\0') {
                cigar += std::to_string(count) + last_op;
            }
            last_op = op;
            count = 1;
        }
    }
    if (last_op != '\0') {
        cigar += std::to_string(count) + last_op;
    }
    return cigar;
}

// Convert string CIGAR format to edit_cigar_t
void wfa_string_to_edit_cigar(const std::string& cigar_str, wflign_cigar_t* edit_cigar) {
    // Clear existing cigar
    edit_cigar->begin_offset = 0;
    edit_cigar->end_offset = 0;
    
    std::vector<char> ops;
    size_t i = 0;
    while (i < cigar_str.size()) {
        // Parse count
        int count = 0;
        while (i < cigar_str.size() && std::isdigit(static_cast<unsigned char>(cigar_str[i]))) {
            count = count * 10 + (cigar_str[i] - '0');
            i++;
        }
        if (i >= cigar_str.size()) break;
        
        // Get operation and convert M to = since we know these are exact matches
        char op = cigar_str[i++];
        if (op == 'M') {
            op = '=';  // Convert M to = since we're dealing with exact matches
        }
        
        // Add operation count times
        for (int j = 0; j < count; j++) {
            edit_cigar->cigar_ops[edit_cigar->end_offset++] = op;
        }
    }
}

} // namespace wflign

