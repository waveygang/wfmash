#pragma once

#include <string>
#include <cctype>
#include <utility>
#include <vector>
#include "../deps/WFA2-lib/bindings/cpp/WFAligner.hpp"

namespace wflign {

// Forward declare WFAlignerEdit functions we'll use
static std::string wfa_edit_cigar_to_string(const wfa::WFAlignerEdit& wf_aligner);
static void wfa_string_to_edit_cigar(const std::string& cigar_str, wfa::WFAlignerEdit* wf_aligner);

// Try to swap start pattern: "N= DlenD" => "DlenD N="
static std::string try_swap_start_pattern(
    const std::string &cigar,
    const std::string &query_seq,
    const std::string &target_seq,
    int64_t query_start,
    int64_t target_start);

// Try to swap end pattern: "DlenD N=" => "N= DlenD"  
static std::string try_swap_end_pattern(
    const std::string &cigar,
    const std::string &target_seq,
    const std::string &query_seq,
    int64_t query_start,
    int64_t target_start);

// Drop leading/trailing deletions if the rest is purely '='
static std::string drop_leading_trailing_deletions_if_all_eq(
    const std::string &cigar);

// Helper functions for CIGAR format conversion
static std::string wfa_edit_cigar_to_string(const wfa::WFAlignerEdit& wf_aligner);
static void wfa_string_to_edit_cigar(const std::string& cigar_str, wfa::WFAlignerEdit* wf_aligner);

// Helper functions
static std::string merge_cigar_ops(const std::string &cigar);
static bool sequences_match(
    const std::string &query_seq,
    const std::string &target_seq,
    int64_t query_start,
    int64_t target_start,
    int N);
static bool extract_first_two_ops(
    const std::string &cigar,
    int &count1, char &op1,
    int &count2, char &op2, 
    size_t &second_op_start,
    size_t &second_op_end);
static bool extract_last_two_ops(
    const std::string &cigar,
    int &count2, char &op2,
    int &count1, char &op1,
    size_t &second_last_op_start,
    size_t &last_op_start);
static std::pair<int64_t,int64_t> alignment_end_coords(
    const std::string &cigar,
    int64_t query_start,
    int64_t target_start);
static bool verify_cigar_alignment(
    const std::string &cigar,
    const std::string &query_seq,
    const std::string &target_seq,
    int64_t query_start,
    int64_t target_start);

} // namespace wflign
