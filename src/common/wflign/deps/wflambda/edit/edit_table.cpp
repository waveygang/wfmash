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
 * DESCRIPTION: DP-table to compute edit alignment
 */

#include "edit/edit_table.hpp"
#include "edit/edit_cigar.hpp"

namespace wflambda {

/*
 * DP-Table Setup
 */
void edit_table_allocate(
    edit_table_t* const edit_table,
    const int pattern_length,
    const int text_length,
    mm_allocator_t* const mm_allocator) {
  // Allocate DP table
  int h;
  edit_table->num_rows = pattern_length + 1;
  const int num_columns = text_length + 1;
  edit_table->num_columns = num_columns;
  edit_table->columns = (int**)mm_allocator_malloc(mm_allocator,(text_length+1)*sizeof(int*)); // Columns
  for (h=0;h<num_columns;++h) {
    edit_table->columns[h] = mm_allocator_calloc(mm_allocator,pattern_length+1,int,false); // Rows
  }
  // CIGAR
  edit_cigar_allocate(&edit_table->edit_cigar,pattern_length,text_length,mm_allocator);
}
void edit_table_free(
    edit_table_t* const edit_table,
    mm_allocator_t* const mm_allocator) {
  const int num_columns = edit_table->num_columns;
  int h;
  for (h=0;h<num_columns;++h) {
    mm_allocator_free(mm_allocator,edit_table->columns[h]);
  }
  mm_allocator_free(mm_allocator,edit_table->columns);
  // CIGAR
  edit_cigar_free(&edit_table->edit_cigar,mm_allocator);
}
/*
 * Display
 */
void edit_table_print(
    FILE* const stream,
    const edit_table_t* const edit_table,
    const char* const pattern,
    const char* const text) {
  // Parameters
  int** const dp = edit_table->columns;
  int i, j;
  // Print Header
  fprintf(stream,"       ");
  for (i=0;i<edit_table->num_columns-1;++i) {
    fprintf(stream,"  %c  ",text[i]);
  }
  fprintf(stream,"\n ");
  for (i=0;i<edit_table->num_columns;++i) {
    fprintf(stream," %3d ",i);
  }
  fprintf(stream,"\n ");
  for (i=0;i<edit_table->num_columns;++i) {
    if (-1 < dp[i][0] && dp[i][0] < 10000) {
      fprintf(stream," %3d ",dp[i][0]);
    } else {
      fprintf(stream,"   * ");
    }
  }
  fprintf(stream,"\n");
  // Print Rows
  for (i=1;i<edit_table->num_rows;++i) {
    fprintf(stream,"%c",pattern[i-1]);
    for (j=0;j<edit_table->num_columns;++j) {
      if (-1 < dp[j][i] && dp[j][i] < 10000) {
        fprintf(stream," %3d ",dp[j][i]);
      } else {
        fprintf(stream,"   * ");
      }
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}

}
