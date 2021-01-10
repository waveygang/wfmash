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
 * DESCRIPTION: WFA support utilities
 */

#ifndef AFFINE_WAVEFRONT_UTILS_H_
#define AFFINE_WAVEFRONT_UTILS_H_

#include "gap_affine/affine_wavefront.h"

/*
 * Accessors
 */
affine_wavefront_t* affine_wavefronts_get_source_mwavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int score);
affine_wavefront_t* affine_wavefronts_get_source_iwavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int score);
affine_wavefront_t* affine_wavefronts_get_source_dwavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int score);
int affine_wavefronts_diagonal_length(
    affine_wavefronts_t* const affine_wavefronts,
    const int k);
int affine_wavefronts_compute_distance(
    const int pattern_length,
    const int text_length,
    const awf_offset_t offset,
    const int k);

/*
 * Initial Conditions and finalization
 */
void affine_wavefront_initialize(
    affine_wavefronts_t* const affine_wavefronts);
bool affine_wavefront_end_reached(
    affine_wavefronts_t* const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int score);

#endif /* AFFINE_WAVEFRONT_UTILS_H_ */
