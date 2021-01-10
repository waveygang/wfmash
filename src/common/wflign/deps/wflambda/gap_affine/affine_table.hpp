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
 * DESCRIPTION: Gap-affine DP table
 */

#ifndef AFFINE_LAMBDA_TABLE_H_
#define AFFINE_LAMBDA_TABLE_H_

#include "utils/commons.hpp"
#include "edit/edit_cigar.hpp"

namespace wflambda {

typedef struct {
  int M; // Alignment matching/mismatching
  int I; // Alignment ends with a gap in the reference (insertion)
  int D; // Alignment ends with a gap in the read (deletion)
} affine_cell_t;
typedef struct {
  // DP Table
  affine_cell_t** columns;
  int num_rows;
  int num_columns;
  // CIGAR
  edit_cigar_t edit_cigar;
} affine_table_t;

/*
 * Gap-Affine Table Setup
 */
void affine_table_allocate(
    affine_table_t* const table,
    const int pattern_length,
    const int text_length,
    mm_allocator_t* const mm_allocator);
void affine_table_free(
    affine_table_t* const table,
    mm_allocator_t* const mm_allocator);

/*
 * Display
 */
void affine_table_print(
    FILE* const stream,
    const affine_table_t* const table,
    const char* const pattern,
    const char* const text);
void affine_table_print_extended(
    FILE* const stream,
    const affine_table_t* const table,
    const char* const pattern,
    const char* const text);

}

#endif /* AFFINE_LAMBDA_TABLE_H_ */
