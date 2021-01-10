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
 * DESCRIPTION: Gap-lineal alignment algorithms wrapper
 */

#include "benchmark/benchmark_gap_lineal.h"
#include "gap_lineal/nw.h"

/*
 * Benchmark NW
 */
void benchmark_gap_lineal_nw(
    align_input_t* const align_input,
    lineal_penalties_t* const penalties) {
  // Allocate
  edit_table_t edit_table;
  edit_table_allocate(
      &edit_table,align_input->pattern_length,
      align_input->text_length,align_input->mm_allocator);
  // Align
  timer_start(&align_input->timer);
  nw_compute(&edit_table,
      align_input->pattern,align_input->pattern_length,
      align_input->text,align_input->text_length,
      penalties);
  timer_stop(&align_input->timer);
  // Free
  edit_table_free(&edit_table,align_input->mm_allocator);
}
