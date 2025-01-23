#ifndef WFLIGN_CIGAR_UTILS_HPP
#define WFLIGN_CIGAR_UTILS_HPP

#include <string>
#include <sstream>
#include "edit_cigar.hpp"

namespace wflign {

// Convert WFA edit_cigar_t to SAM-style CIGAR string
inline std::string wfa_edit_cigar_to_string(const edit_cigar_t& edit_cigar) {
    std::stringstream cigar;
    int count = 0;
    char last_op = '\0';
    
    for (int i = edit_cigar.begin_offset; i < edit_cigar.end_offset; ++i) {
        const char op = edit_cigar.operations[i];
        if (op == last_op) {
            ++count;
        } else {
            if (count > 0) {
                cigar << count << last_op;
            }
            count = 1;
            last_op = op;
        }
    }
    
    if (count > 0) {
        cigar << count << last_op;
    }
    
    return cigar.str();
}

// Convert SAM-style CIGAR string back to WFA edit_cigar_t
inline void wfa_string_to_edit_cigar(const std::string& cigar_str, edit_cigar_t* edit_cigar) {
    edit_cigar_clear(edit_cigar);
    
    int num = 0;
    for (char c : cigar_str) {
        if (std::isdigit(c)) {
            num = num * 10 + (c - '0');
        } else {
            for (int i = 0; i < num; ++i) {
                edit_cigar_add_operation(edit_cigar, c);
            }
            num = 0;
        }
    }
}

} // namespace wflign

#endif // WFLIGN_CIGAR_UTILS_HPP
