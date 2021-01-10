/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Edit cigar data-structure (match/mismatch/insertion/deletion)
 */

#ifndef EDIT_CIGAR_H_
#define EDIT_CIGAR_H_

#include "utils/commons.hpp"
#include "system/mm_allocator.hpp"
#include "gap_lineal/lineal_penalties.hpp"
#include "gap_affine/affine_penalties.hpp"

namespace wflambda {

/*
 * CIGAR
 */
typedef struct {
  char* operations;
  int max_operations;
  int begin_offset;
  int end_offset;
  int score;
} edit_cigar_t;

/*
 * Distance metrics
 */
typedef enum {
  edit,
  gap_lineal,
  gap_affine
} distance_metric_t;

/*
 * Setup
 */
void edit_cigar_allocate(
    edit_cigar_t* const edit_cigar,
    const int pattern_length,
    const int text_length,
    mm_allocator_t* const mm_allocator);
void edit_cigar_clear(
    edit_cigar_t* const edit_cigar);
void edit_cigar_free(
    edit_cigar_t* const edit_cigar,
    mm_allocator_t* const mm_allocator);

/*
 * Accessors
 */
int edit_cigar_get_matches(
    edit_cigar_t* const edit_cigar);
void edit_cigar_add_mismatches(
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    edit_cigar_t* const edit_cigar);

/*
 * Score
 */
int edit_cigar_score_edit(
    edit_cigar_t* const edit_cigar);
int edit_cigar_score_gap_lineal(
    edit_cigar_t* const edit_cigar,
    lineal_penalties_t* const penalties);
int edit_cigar_score_gap_affine(
    edit_cigar_t* const edit_cigar,
    affine_penalties_t* const penalties);

/*
 * Utils
 */
int edit_cigar_cmp(
    edit_cigar_t* const edit_cigar_a,
    edit_cigar_t* const edit_cigar_b);
void edit_cigar_copy(
    edit_cigar_t* const edit_cigar_dst,
    edit_cigar_t* const edit_cigar_src);
bool edit_cigar_check_alignment(
    FILE* const stream,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    edit_cigar_t* const edit_cigar,
    const bool verbose);

/*
 * Display
 */
void edit_cigar_print(
    FILE* const stream,
    edit_cigar_t* const edit_cigar);
void edit_cigar_print_pretty(
    FILE* const stream,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    edit_cigar_t* const edit_cigar,
    mm_allocator_t* const mm_allocator);

}

#endif /* EDIT_CIGAR_H_ */
