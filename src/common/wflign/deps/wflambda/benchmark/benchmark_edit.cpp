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
 * DESCRIPTION: Edit-distance alignment algorithms wrapper
 */

#include "benchmark/benchmark_utils.hpp"
#include "edit/edit_dp.hpp"

namespace wflambda {

/*
 * Benchmark Edit
 */
void benchmark_edit_dp(
    align_input_t* const align_input) {
  // Parameters
  edit_table_t edit_table;
  // Allocate
  edit_table_allocate(
      &edit_table,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  edit_dp_compute(&edit_table,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length);
  timer_stop(&align_input->timer);
  // Debug alignment
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&edit_table.edit_cigar);
  }
  // Free
  edit_table_free(&edit_table,align_input->mm_allocator);
}
void benchmark_edit_dp_banded(
    align_input_t* const align_input,
    const int bandwidth) {
  // Parameters
  edit_table_t edit_table;
  // Allocate
  edit_table_allocate(
      &edit_table,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  edit_dp_compute_banded(&edit_table,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,bandwidth);
  timer_stop(&align_input->timer);
  // Debug alignment
  if (align_input->debug_flags) {
    benchmark_check_alignment(align_input,&edit_table.edit_cigar);
  }
  // Free
  edit_table_free(&edit_table,align_input->mm_allocator);
}

}
