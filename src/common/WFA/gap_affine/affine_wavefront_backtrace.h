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
 * DESCRIPTION: WFA extend backtrace component
 */

#ifndef AFFINE_WAVEFRONT_BACKTRACE_H_
#define AFFINE_WAVEFRONT_BACKTRACE_H_

#include "gap_affine/affine_wavefront.h"

/*
 * Sequences DTO
 */
typedef struct {
  char* pattern;
  int pattern_length;
  char* text;
  int text_length;
} alignment_sequences_t;

/*
 * WF type
 */
typedef enum {
  backtrace_wavefront_M = 0,
  backtrace_wavefront_I = 1,
  backtrace_wavefront_D = 2
} backtrace_wavefront_type;

/*
 * Backtrace
 */
void affine_wavefronts_backtrace(
    affine_wavefronts_t* const affine_wavefronts,
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int alignment_score);

#endif /* AFFINE_WAVEFRONT_BACKTRACE_H_ */
