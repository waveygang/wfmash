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
 * DESCRIPTION: WFA support functions for wavefront-reduction
 *   strategies (adaptive or banded strategies)
 */

#include "gap_affine/affine_wavefront_reduction.h"

/*
 * Setup
 */
void affine_wavefronts_reduction_set_none(
    affine_wavefronts_reduction_t* const wavefronts_reduction) {
  wavefronts_reduction->reduction_strategy = wavefronts_reduction_none;
}
void affine_wavefronts_reduction_set_dynamic(
    affine_wavefronts_reduction_t* const wavefronts_reduction,
    const int min_wavefront_length,
    const int max_distance_threshold) {
  wavefronts_reduction->reduction_strategy = wavefronts_reduction_dynamic;
  wavefronts_reduction->min_wavefront_length = min_wavefront_length;
  wavefronts_reduction->max_distance_threshold = max_distance_threshold;
}
