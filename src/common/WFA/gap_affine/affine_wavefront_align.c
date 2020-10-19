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
 * DESCRIPTION: WFA main algorithm
 */

#include "affine_wavefront_align.h"
#include "gap_affine/affine_wavefront_backtrace.h"
#include "gap_affine/affine_wavefront_display.h"
#include "gap_affine/affine_wavefront_extend.h"
#include "gap_affine/affine_wavefront_utils.h"
#include "utils/string_padded.h"

/*
 * Fetch & allocate wavefronts
 */
void affine_wavefronts_fetch_wavefronts(
    affine_wavefronts_t* const affine_wavefronts,
    affine_wavefront_set* const wavefront_set,
    const int score) {
  // Compute scores
  const affine_penalties_t* const wavefront_penalties = &(affine_wavefronts->penalties.wavefront_penalties);
  const int mismatch_score = score - wavefront_penalties->mismatch;
  const int gap_open_score = score - wavefront_penalties->gap_opening - wavefront_penalties->gap_extension;
  const int gap_extend_score = score - wavefront_penalties->gap_extension;
  // Fetch wavefronts
  wavefront_set->in_mwavefront_sub = affine_wavefronts_get_source_mwavefront(affine_wavefronts,mismatch_score);
  wavefront_set->in_mwavefront_gap = affine_wavefronts_get_source_mwavefront(affine_wavefronts,gap_open_score);
  wavefront_set->in_iwavefront_ext = affine_wavefronts_get_source_iwavefront(affine_wavefronts,gap_extend_score);
  wavefront_set->in_dwavefront_ext = affine_wavefronts_get_source_dwavefront(affine_wavefronts,gap_extend_score);
}
void affine_wavefronts_allocate_wavefronts(
    affine_wavefronts_t* const affine_wavefronts,
    affine_wavefront_set* const wavefront_set,
    const int score,
    const int lo_effective,
    const int hi_effective) {
  // Allocate M-Wavefront
  wavefront_set->out_mwavefront =
      affine_wavefronts_allocate_wavefront(affine_wavefronts,lo_effective,hi_effective);
  affine_wavefronts->mwavefronts[score] = wavefront_set->out_mwavefront;
  // Allocate I-Wavefront
  if (!wavefront_set->in_mwavefront_gap->null || !wavefront_set->in_iwavefront_ext->null) {
    wavefront_set->out_iwavefront =
        affine_wavefronts_allocate_wavefront(affine_wavefronts,lo_effective,hi_effective);
    affine_wavefronts->iwavefronts[score] = wavefront_set->out_iwavefront;
  } else {
    wavefront_set->out_iwavefront = NULL;
  }
  // Allocate D-Wavefront
  if (!wavefront_set->in_mwavefront_gap->null || !wavefront_set->in_dwavefront_ext->null) {
    wavefront_set->out_dwavefront =
        affine_wavefronts_allocate_wavefront(affine_wavefronts,lo_effective,hi_effective);
    affine_wavefronts->dwavefronts[score] = wavefront_set->out_dwavefront;
  } else {
    wavefront_set->out_dwavefront = NULL;
  }
}
void affine_wavefronts_compute_limits(
    affine_wavefronts_t* const affine_wavefronts,
    const affine_wavefront_set* const wavefront_set,
    const int score,
    int* const lo_effective,
    int* const hi_effective) {
  // Set limits (min_lo)
  int lo = wavefront_set->in_mwavefront_sub->lo;
  if (lo > wavefront_set->in_mwavefront_gap->lo) lo = wavefront_set->in_mwavefront_gap->lo;
  if (lo > wavefront_set->in_iwavefront_ext->lo) lo = wavefront_set->in_iwavefront_ext->lo;
  if (lo > wavefront_set->in_dwavefront_ext->lo) lo = wavefront_set->in_dwavefront_ext->lo;
  --lo;
  // Set limits (max_hi)
  int hi = wavefront_set->in_mwavefront_sub->hi;
  if (hi < wavefront_set->in_mwavefront_gap->hi) hi = wavefront_set->in_mwavefront_gap->hi;
  if (hi < wavefront_set->in_iwavefront_ext->hi) hi = wavefront_set->in_iwavefront_ext->hi;
  if (hi < wavefront_set->in_dwavefront_ext->hi) hi = wavefront_set->in_dwavefront_ext->hi;
  ++hi;
  // Set effective limits values
  *hi_effective = hi;
  *lo_effective = lo;
}
/*
 * Compute wavefront offsets
 */
#define AFFINE_WAVEFRONT_DECLARE(wavefront,prefix) \
  const awf_offset_t* const prefix ## _offsets = wavefront->offsets; \
  const int prefix ## _hi = wavefront->hi; \
  const int prefix ## _lo = wavefront->lo
#define AFFINE_WAVEFRONT_COND_FETCH(prefix,index,value) \
  (prefix ## _lo <= (index) && (index) <= prefix ## _hi) ? (value) : AFFINE_WAVEFRONT_OFFSET_NULL
/*
 * Compute wavefront offsets
 */
void affine_wavefronts_compute_offsets_idm(
    affine_wavefronts_t* const affine_wavefronts,
    const affine_wavefront_set* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap,m_gap);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_iwavefront_ext,i_ext);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_dwavefront_ext,d_ext);
  awf_offset_t* const out_ioffsets = wavefront_set->out_iwavefront->offsets;
  awf_offset_t* const out_doffsets = wavefront_set->out_dwavefront->offsets;
  awf_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  // Compute loop peeling offset (min_hi)
  int min_hi = wavefront_set->in_mwavefront_sub->hi;
  if (!wavefront_set->in_mwavefront_gap->null && min_hi > wavefront_set->in_mwavefront_gap->hi-1) min_hi = wavefront_set->in_mwavefront_gap->hi-1;
  if (!wavefront_set->in_iwavefront_ext->null && min_hi > wavefront_set->in_iwavefront_ext->hi+1) min_hi = wavefront_set->in_iwavefront_ext->hi+1;
  if (!wavefront_set->in_dwavefront_ext->null && min_hi > wavefront_set->in_dwavefront_ext->hi-1) min_hi = wavefront_set->in_dwavefront_ext->hi-1;
  // Compute loop peeling offset (max_lo)
  int max_lo = wavefront_set->in_mwavefront_sub->lo;
  if (!wavefront_set->in_mwavefront_gap->null && max_lo < wavefront_set->in_mwavefront_gap->lo+1) max_lo = wavefront_set->in_mwavefront_gap->lo+1;
  if (!wavefront_set->in_iwavefront_ext->null && max_lo < wavefront_set->in_iwavefront_ext->lo+1) max_lo = wavefront_set->in_iwavefront_ext->lo+1;
  if (!wavefront_set->in_dwavefront_ext->null && max_lo < wavefront_set->in_dwavefront_ext->lo-1) max_lo = wavefront_set->in_dwavefront_ext->lo-1;
  // Compute score wavefronts (prologue)
  int k;
  for (k=lo;k<max_lo;++k) {
    // Update I
    const awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap,k-1,m_gap_offsets[k-1]);
    const awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(i_ext,k-1,i_ext_offsets[k-1]);
    const awf_offset_t ins = MAX(ins_g,ins_i) + 1;
    out_ioffsets[k] = ins;
    // Update D
    const awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap,k+1,m_gap_offsets[k+1]);
    const awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(d_ext,k+1,d_ext_offsets[k+1]);
    const awf_offset_t del = MAX(del_g,del_d);
    out_doffsets[k] = del;
    // Update M
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(del,MAX(sub,ins));
  }
  // Compute score wavefronts (core)
#if defined(__clang__)
  #pragma clang loop vectorize(enable)
#elif defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=max_lo;k<=min_hi;++k) {
    // Update I
    const awf_offset_t m_gapi_value = m_gap_offsets[k-1];
    const awf_offset_t i_ext_value = i_ext_offsets[k-1];
    const awf_offset_t ins = MAX(m_gapi_value,i_ext_value) + 1;
    out_ioffsets[k] = ins;
    // Update D
    const awf_offset_t m_gapd_value = m_gap_offsets[k+1];
    const awf_offset_t d_ext_value = d_ext_offsets[k+1];
    const awf_offset_t del = MAX(m_gapd_value,d_ext_value);
    out_doffsets[k] = del;
    // Update M
    const awf_offset_t sub = m_sub_offsets[k] + 1;
    out_moffsets[k] = MAX(del,MAX(sub,ins));
  }
  // Compute score wavefronts (epilogue)
  for (k=min_hi+1;k<=hi;++k) {
    // Update I
    const awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap,k-1,m_gap_offsets[k-1]);
    const awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(i_ext,k-1,i_ext_offsets[k-1]);
    const awf_offset_t ins = MAX(ins_g,ins_i) + 1;
    out_ioffsets[k] = ins;
    // Update D
    const awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap,k+1,m_gap_offsets[k+1]);
    const awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(d_ext,k+1,d_ext_offsets[k+1]);
    const awf_offset_t del = MAX(del_g,del_d);
    out_doffsets[k] = del;
    // Update M
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(del,MAX(sub,ins));
  }
}
void affine_wavefronts_compute_offsets_im(
    affine_wavefronts_t* const affine_wavefronts,
    const affine_wavefront_set* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap,m_gap);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_iwavefront_ext,i_ext);
  awf_offset_t* const out_ioffsets = wavefront_set->out_iwavefront->offsets;
  awf_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  // Compute score wavefronts
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=lo;k<=hi;++k) {
    // Update I
    const awf_offset_t ins_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap,k-1,m_gap_offsets[k-1]);
    const awf_offset_t ins_i = AFFINE_WAVEFRONT_COND_FETCH(i_ext,k-1,i_ext_offsets[k-1]);
    const awf_offset_t ins = MAX(ins_g,ins_i) + 1;
    out_ioffsets[k] = ins;
    // Update M
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(ins,sub);
  }
}
void affine_wavefronts_compute_offsets_dm(
    affine_wavefronts_t* const affine_wavefronts,
    const affine_wavefront_set* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_gap,m_gap);
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_dwavefront_ext,d_ext);
  awf_offset_t* const out_doffsets = wavefront_set->out_dwavefront->offsets;
  awf_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  // Compute score wavefronts
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=lo;k<=hi;++k) {
    // Update D
    const awf_offset_t del_g = AFFINE_WAVEFRONT_COND_FETCH(m_gap,k+1,m_gap_offsets[k+1]);
    const awf_offset_t del_d = AFFINE_WAVEFRONT_COND_FETCH(d_ext,k+1,d_ext_offsets[k+1]);
    const awf_offset_t del = MAX(del_g,del_d);
    out_doffsets[k] = del;
    // Update M
    const awf_offset_t sub = AFFINE_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
    out_moffsets[k] = MAX(del,sub);
  }
}
void affine_wavefronts_compute_offsets_m(
    affine_wavefronts_t* const affine_wavefronts,
    const affine_wavefront_set* const wavefront_set,
    const int lo,
    const int hi) {
  // Parameters
  AFFINE_WAVEFRONT_DECLARE(wavefront_set->in_mwavefront_sub,m_sub);
  awf_offset_t* const out_moffsets = wavefront_set->out_mwavefront->offsets;
  // Compute score wavefronts
  int k;
#if defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC ivdep
#else
  #pragma ivdep
#endif
  for (k=lo;k<=hi;++k) {
    // Update M
    out_moffsets[k] = AFFINE_WAVEFRONT_COND_FETCH(m_sub,k,m_sub_offsets[k]+1);
  }
}
/*
 * Compute wavefront
 */
void affine_wavefronts_compute_wavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int score) {
  // Select wavefronts
  affine_wavefront_set wavefront_set;
  affine_wavefronts_fetch_wavefronts(affine_wavefronts,&wavefront_set,score);
  // Check null wavefronts
  if (wavefront_set.in_mwavefront_sub->null &&
      wavefront_set.in_mwavefront_gap->null &&
      wavefront_set.in_iwavefront_ext->null &&
      wavefront_set.in_dwavefront_ext->null) {
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_steps_null,1);
    return;
  }
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_mwavefront_sub->null?1:0));
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_mwavefront_gap->null?1:0));
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_iwavefront_ext->null?1:0));
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_null_used,(wavefront_set.in_dwavefront_ext->null?1:0));
  // Set limits
  int hi, lo;
  affine_wavefronts_compute_limits(affine_wavefronts,&wavefront_set,score,&lo,&hi);
  // Allocate score-wavefronts
  affine_wavefronts_allocate_wavefronts(affine_wavefronts,&wavefront_set,score,lo,hi);
  // Compute WF
  const int kernel = ((wavefront_set.out_iwavefront!=NULL) << 1) | (wavefront_set.out_dwavefront!=NULL);
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_compute_kernel[kernel],1);
  switch (kernel) {
    case 3: // 11b
      affine_wavefronts_compute_offsets_idm(affine_wavefronts,&wavefront_set,lo,hi);
      break;
    case 2: // 10b
      affine_wavefronts_compute_offsets_im(affine_wavefronts,&wavefront_set,lo,hi);
      break;
    case 1: // 01b
      affine_wavefronts_compute_offsets_dm(affine_wavefronts,&wavefront_set,lo,hi);
      break;
    case 0: // 00b
      affine_wavefronts_compute_offsets_m(affine_wavefronts,&wavefront_set,lo,hi);
      break;
  }
  // Account for WF operations performed
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_operations,hi-lo+1);
  // DEBUG
#ifdef AFFINE_WAVEFRONT_DEBUG
  // Copy offsets base before extension (for display purposes)
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront!=NULL) {
    int k;
    for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
      mwavefront->offsets_base[k] = mwavefront->offsets[k];
    }
  }
#endif
}
/*
 * Computation using Wavefronts
 */
void affine_wavefronts_align(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Init padded strings
  strings_padded_t* const strings_padded =
      strings_padded_new_rhomb(
          pattern,pattern_length,text,text_length,
          AFFINE_WAVEFRONT_PADDING,affine_wavefronts->mm_allocator);
  // Initialize wavefront
  affine_wavefront_initialize(affine_wavefronts);
  // Compute wavefronts for increasing score
  int score = 0;
  while (true) {
    // Exact extend s-wavefront
    affine_wavefronts_extend_wavefront_packed(
        affine_wavefronts,strings_padded->pattern_padded,pattern_length,
        strings_padded->text_padded,text_length,score);
    // Exit condition
    if (affine_wavefront_end_reached(affine_wavefronts,pattern_length,text_length,score)) {
      // Backtrace & check alignment reached
      affine_wavefronts_backtrace(
          affine_wavefronts,strings_padded->pattern_padded,pattern_length,
          strings_padded->text_padded,text_length,score);
      break;
    }
    // Update all wavefronts
    ++score; // Increase score
    affine_wavefronts_compute_wavefront(
        affine_wavefronts,strings_padded->pattern_padded,pattern_length,
        strings_padded->text_padded,text_length,score);
    // DEBUG
    //affine_wavefronts_debug_step(affine_wavefronts,pattern,text,score);
    WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_steps,1);
  }
  // DEBUG
  //affine_wavefronts_debug_step(affine_wavefronts,pattern,text,score);
  WAVEFRONT_STATS_COUNTER_ADD(affine_wavefronts,wf_score,score); // STATS
  // Free
  strings_padded_delete(strings_padded);
}

