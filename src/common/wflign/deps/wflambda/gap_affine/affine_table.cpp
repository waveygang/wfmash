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

#include "edit/edit_table.hpp"
#include "gap_affine/affine_table.hpp"

namespace wflambda {

/*
 * Gap-Affine Table Setup
 */
void affine_table_allocate(
    affine_table_t* const table,
    const int pattern_length,
    const int text_length,
    mm_allocator_t* const mm_allocator) {
  // Allocate DP table
  int h;
  table->num_rows = pattern_length+1;
  const int num_columns = text_length+1;
  table->num_columns = num_columns;
  table->columns = (affine_cell_t**)mm_allocator_malloc(mm_allocator,(text_length+1)*sizeof(affine_cell_t*)); // Columns
  for (h=0;h<num_columns;++h) {
    table->columns[h] = mm_allocator_calloc(mm_allocator,pattern_length+1,affine_cell_t,false); // Rows
  }
  // CIGAR
  edit_cigar_allocate(&table->edit_cigar,pattern_length,text_length,mm_allocator);
}
void affine_table_free(
    affine_table_t* const table,
    mm_allocator_t* const mm_allocator) {
  const int num_columns = table->num_columns;
  int h;
  for (h=0;h<num_columns;++h) {
    mm_allocator_free(mm_allocator,table->columns[h]);
  }
  mm_allocator_free(mm_allocator,table->columns);
  // CIGAR
  edit_cigar_free(&table->edit_cigar,mm_allocator);
}
/*
 * Display
 */
#define AFFINE_LAMBDA_TABLE_PRINT_VALUE(value)  \
  if (value >= 0 && value < SCORE_MAX) { \
    fprintf(stream,"%2d",value); \
  } else { \
    fprintf(stream," *"); \
  }
void affine_table_print(
    FILE* const stream,
    const affine_table_t* const table,
    const char* const pattern,
    const char* const text) {
  // Parameters
  affine_cell_t** const dp = table->columns;
  int i, j;
  // Print Header
  fprintf(stream,"     ");
  for (i=0;i<table->num_columns-1;++i) {
    fprintf(stream," %c ",text[i]);
  }
  fprintf(stream,"\n ");
  for (i=0;i<table->num_columns;++i) {
    fprintf(stream," ");
    AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[i][0].M);
  }
  fprintf(stream,"\n");
  // Print Rows
  for (i=1;i<table->num_rows;++i) {
    fprintf(stream,"%c",pattern[i-1]);
    for (j=0;j<table->num_columns;++j) {
      fprintf(stream," ");
      AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[j][i].M);
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}
void affine_table_print_extended(
    FILE* const stream,
    const affine_table_t* const table,
    const char* const pattern,
    const char* const text) {
  // Parameters
  affine_cell_t** const dp = table->columns;
  int i, j;
  // Print Header
  fprintf(stream,"         ");
  for (i=0;i<table->num_columns-1;++i) {
    fprintf(stream,"     %c     ",text[i]);
  }
  fprintf(stream,"\n ");
  for (i=0;i<table->num_columns;++i) {
    fprintf(stream," ");
    AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[i][0].M);
    fprintf(stream,"{");
    AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[i][0].I);
    fprintf(stream,",");
    AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[i][0].D);
    fprintf(stream,"} ");
  }
  fprintf(stream,"\n");
  // Print Rows
  for (i=1;i<table->num_rows;++i) {
    fprintf(stream,"%c",pattern[i-1]);
    for (j=0;j<table->num_columns;++j) {
      fprintf(stream," ");
      AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[j][i].M);
      fprintf(stream,"{");
      AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[j][i].I);
      fprintf(stream,",");
      AFFINE_LAMBDA_TABLE_PRINT_VALUE(dp[j][i].D);
      fprintf(stream,"} ");
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}

}
