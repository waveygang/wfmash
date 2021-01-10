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
 * DESCRIPTION: WFA extend exact-matches component
 */

#include "gap_affine/affine_wavefront_display.h"
#include "gap_affine/affine_wavefront_extend.h"
#include "gap_affine/affine_wavefront_reduction.h"
#include "gap_affine/affine_wavefront_utils.h"
#include "utils/string_padded.h"

/*
 * Reduce wavefront
 */
void affine_wavefronts_reduce_wavefront_offsets(
    affine_wavefronts_t* const affine_wavefronts,
    affine_wavefront_t* const wavefront,
    const int pattern_length,
    const int text_length,
    const int min_distance,
    const int max_distance_threshold,
    const int alignment_k) {
  // Parameters
  const awf_offset_t* const offsets = wavefront->offsets;
  int k;
  // Reduce from bottom
  const int top_limit = MIN(alignment_k-1,wavefront->hi);
  for (k=wavefront->lo;k<top_limit;++k) {
    const int distance = affine_wavefronts_compute_distance(pattern_length,text_length,offsets[k],k);
    if (distance - min_distance  <= max_distance_threshold) break;
    ++(wavefront->lo);
  }
  // Reduce from top
  const int botton_limit = MAX(alignment_k+1,wavefront->lo);
  for (k=wavefront->hi;k>botton_limit;--k) {
    const int distance = affine_wavefronts_compute_distance(pattern_length,text_length,offsets[k],k);
    if (distance - min_distance <= max_distance_threshold) break;
    --(wavefront->hi);
  }
  // Check hi/lo range
  if (wavefront->lo > wavefront->hi) {
    wavefront->null = true;
  }
  // STATS
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_reduced_cells,
      (wavefront->hi_base-wavefront->hi)+(wavefront->lo-wavefront->lo_base));
}
void affine_wavefronts_reduce_wavefronts(
    affine_wavefronts_t* const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int score) {
  // Parameters
  const int min_wavefront_length = affine_wavefronts->reduction.min_wavefront_length;
  const int max_distance_threshold = affine_wavefronts->reduction.max_distance_threshold;
  const int alignment_k = AFFINE_WAVEFRONT_DIAGONAL(text_length,pattern_length);
  // Fetch m-wavefront
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront==NULL) return;
  if ((mwavefront->hi - mwavefront->lo + 1) < min_wavefront_length) return;
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_reduction,1); // STATS
  // Compute min-distance
  const awf_offset_t* const offsets = mwavefront->offsets;
  int min_distance = MAX(pattern_length,text_length);
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    const int distance = affine_wavefronts_compute_distance(pattern_length,text_length,offsets[k],k);
    min_distance = MIN(min_distance,distance);
  }
  // Reduce m-wavefront
  affine_wavefronts_reduce_wavefront_offsets(
      affine_wavefronts,mwavefront,pattern_length,text_length,
      min_distance,max_distance_threshold,alignment_k);
  // Reduce i-wavefront
  affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[score];
  if (iwavefront!=NULL) {
    if (mwavefront->lo > iwavefront->lo) iwavefront->lo = mwavefront->lo;
    if (mwavefront->hi < iwavefront->hi) iwavefront->hi = mwavefront->hi;
    if (iwavefront->lo > iwavefront->hi) iwavefront->null = true;
  }
  // Reduce d-wavefront
  affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[score];
  if (dwavefront!=NULL) {
    if (mwavefront->lo > dwavefront->lo) dwavefront->lo = mwavefront->lo;
    if (mwavefront->hi < dwavefront->hi) dwavefront->hi = mwavefront->hi;
    if (dwavefront->lo > dwavefront->hi) dwavefront->null = true;
  }
}
/*
 * Wavefront offset extension comparing characters
 */
void affine_wavefronts_extend_mwavefront_epiloge(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int pattern_length,
    const int text_length) {
  // STATS
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_extensions,
      affine_wavefronts->mwavefronts[score]->hi-affine_wavefronts->mwavefronts[score]->lo+1);
  // DEBUG
#ifdef AFFINE_WAVEFRONT_DEBUG
  // Parameters
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    if (mwavefront->offsets[k] >= 0) {
      awf_offset_t offset;
      for (offset=mwavefront->offsets_base[k];offset<=mwavefront->offsets[k];++offset) {
        affine_wavefronts_set_edit_table(affine_wavefronts,pattern_length,text_length,k,offset,score);
      }
    }
  }
#endif
}
void affine_wavefronts_extend_mwavefront_compute_packed(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Fetch m-wavefront
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront==NULL) return;
  // Extend diagonally each wavefront point
  awf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    // Fetch positions & offset
    const awf_offset_t offset = offsets[k];
    const int v = AFFINE_WAVEFRONT_V(k,offset);
    const int h = AFFINE_WAVEFRONT_H(k,offset);
    // Fetch pattern/text blocks
    uint64_t* pattern_blocks = (uint64_t*)(pattern+v);
    uint64_t* text_blocks = (uint64_t*)(text+h);
    uint64_t pattern_block = *pattern_blocks;
    uint64_t text_block = *text_blocks;
    // Compare 64-bits blocks
    uint64_t cmp = pattern_block ^ text_block;
    while (__builtin_expect(!cmp,0)) {
      // Increment offset (full block)
      offsets[k] += 8;
      // Next blocks
      ++pattern_blocks;
      ++text_blocks;
      // Fetch
      pattern_block = *pattern_blocks;
      text_block = *text_blocks;
      // Compare
      cmp = pattern_block ^ text_block;
      WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_extend_inner_loop,1); // STATS
    }
    // Count equal characters
    const int equal_right_bits = __builtin_ctzl(cmp);
    const int equal_chars = DIV_FLOOR(equal_right_bits,8);
    // Increment offset
    offsets[k] += equal_chars;
  }
  // DEBUG
  affine_wavefronts_extend_mwavefront_epiloge(
      affine_wavefronts,score,pattern_length,text_length);
}
void affine_wavefronts_extend_mwavefront_compute(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Fetch m-wavefront
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront==NULL) return;
  // Extend diagonally each wavefront point
  awf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    // Exact extend
    const awf_offset_t offset = offsets[k];
    int v = AFFINE_WAVEFRONT_V(k,offset);
    int h = AFFINE_WAVEFRONT_H(k,offset);
    while (pattern[v++]==text[h++]) {
      ++(offsets[k]);
    }
  }
  // DEBUG
  affine_wavefronts_extend_mwavefront_epiloge(
      affine_wavefronts,score,pattern_length,text_length);
}
/*
 * Gap-Affine Wavefront exact extension
 */
void affine_wavefronts_extend_wavefront_packed(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Extend wavefront
  affine_wavefronts_extend_mwavefront_compute_packed(
      affine_wavefronts,pattern,pattern_length,
      text,text_length,score);
  // Reduce wavefront dynamically
  if (affine_wavefronts->reduction.reduction_strategy == wavefronts_reduction_dynamic) {
    affine_wavefronts_reduce_wavefronts(
        affine_wavefronts,pattern_length,
        text_length,score);
  }
}
