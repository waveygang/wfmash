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
 * DESCRIPTION: WFA support functions for handling penalties scores
 */

#ifndef AFFINE_LAMBDA_WAVEFRONT_PENALTIES_H_
#define AFFINE_LAMBDA_WAVEFRONT_PENALTIES_H_

#include "utils/commons.hpp"
#include "gap_affine/affine_penalties.hpp"

namespace wflambda {

/*
 * Wavefront Strategy
 */
typedef enum {
  wavefronts_penalties_match_zero,
  wavefronts_penalties_force_zero_match,
  wavefronts_penalties_shifted_penalties,
  wavefronts_penalties_odd_pair_penalties
} wavefronts_penalties_strategy;

/*
 * Wavefront Penalties
 */
typedef struct {
  affine_penalties_t base_penalties;                // Input base Gap-Affine penalties
  affine_penalties_t wavefront_penalties;           // Wavefront Gap-Affine penalties
  wavefronts_penalties_strategy penalties_strategy; // Penalties adaptation strategy
} affine_wavefronts_penalties_t;

/*
 * Setup
 */
void affine_wavefronts_penalties_init(
    affine_wavefronts_penalties_t* const wavefronts_penalties,
    affine_penalties_t* const penalties,
    const wavefronts_penalties_strategy penalties_strategy);

/*
 * Score Adjustment
 */
void affine_penalties_mzero(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties);
void affine_penalties_shift(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties,
    const bool pair_odd_heuristic);

}

#endif /* AFFINE_LAMBDA_WAVEFRONT_PENALTIES_H_ */
