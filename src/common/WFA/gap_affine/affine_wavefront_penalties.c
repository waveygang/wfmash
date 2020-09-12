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

#include "gap_affine/affine_wavefront_penalties.h"

/*
 * Setup
 */
void affine_wavefronts_penalties_init(
    affine_wavefronts_penalties_t* const wavefronts_penalties,
    affine_penalties_t* const penalties,
    const wavefronts_penalties_strategy penalties_strategy) {
  wavefronts_penalties->base_penalties = *penalties;
  wavefronts_penalties->penalties_strategy =
      (penalties->match==0) ? wavefronts_penalties_match_zero : penalties_strategy;
  switch (wavefronts_penalties->penalties_strategy) {
    case wavefronts_penalties_match_zero:
    case wavefronts_penalties_force_zero_match:
      affine_penalties_mzero(penalties,&(wavefronts_penalties->wavefront_penalties));
      break;
    case wavefronts_penalties_shifted_penalties:
      affine_penalties_shift(penalties,&(wavefronts_penalties->wavefront_penalties),false);
      break;
    case wavefronts_penalties_odd_pair_penalties:
      affine_penalties_shift(penalties,&(wavefronts_penalties->wavefront_penalties),true);
      break;
    default:
      break;
  }
}
/*
 * Score Adjustment
 */
void affine_penalties_mzero(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties) {
  // Check base penalties
  if (base_penalties->match > 0) {
    fprintf(stderr,"Match score must be negative or zero (M=%d)\n",base_penalties->match);
    exit(1);
  }
  if (base_penalties->mismatch <= 0 ||
      base_penalties->gap_opening <= 0 ||
      base_penalties->gap_extension <= 0) {
    fprintf(stderr,"Mismatch/Gap scores must be strictly positive (X=%d,O=%d,E=%d)\n",
        base_penalties->mismatch,base_penalties->gap_opening,base_penalties->gap_extension);
    exit(1);
  }
  // Copy base penalties
  *shifted_penalties = *base_penalties;
  // Zero match score
  shifted_penalties->match = 0;
}
void affine_penalties_shift(
    affine_penalties_t* const base_penalties,
    affine_penalties_t* const shifted_penalties,
    const bool pair_odd_heuristic) {
  // Check base penalties
  if (base_penalties->match > 0) {
    fprintf(stderr,"Match score must be negative (M=%d)\n",base_penalties->match);
    exit(1);
  }
  if (base_penalties->mismatch <= 0 ||
      base_penalties->gap_opening <= 0 ||
      base_penalties->gap_extension <= 0) {
    fprintf(stderr,"Mismatch/Gap scores must be strictly positive (X=%d,O=%d,E=%d)\n",
        base_penalties->mismatch,base_penalties->gap_opening,base_penalties->gap_extension);
    exit(1);
  }
  // Copy base penalties
  *shifted_penalties = *base_penalties;
  // Shift to zero match score
  shifted_penalties->match = 0;
  shifted_penalties->mismatch -= base_penalties->match;
  shifted_penalties->gap_opening -= base_penalties->match;
  shifted_penalties->gap_extension -= base_penalties->match;
  // Odd/Pair shift heuristic
  if (pair_odd_heuristic) {
    const bool is_mismatch_pair = ((shifted_penalties->mismatch%2)==0);
    const bool is_gap_opening_pair = ((shifted_penalties->gap_opening%2)==0);
    const bool is_gap_extension_pair = ((shifted_penalties->gap_extension%2)==0);
    const int total_odd = !is_mismatch_pair + !is_gap_opening_pair + !is_gap_extension_pair;
    const int total_pair = is_mismatch_pair + is_gap_opening_pair + is_gap_extension_pair;
    if (total_odd > total_pair) {
      // Shift all to odd
      if (is_mismatch_pair) ++(shifted_penalties->mismatch);
      if (is_gap_opening_pair) ++(shifted_penalties->gap_opening);
      if (is_gap_extension_pair) ++(shifted_penalties->gap_extension);
    } else {
      // Shift all to pair
      if (!is_mismatch_pair) ++(shifted_penalties->mismatch);
      if (!is_gap_opening_pair) ++(shifted_penalties->gap_opening);
      if (!is_gap_extension_pair) ++(shifted_penalties->gap_extension);
    }
  }
}

