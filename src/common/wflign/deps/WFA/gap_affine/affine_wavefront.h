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
 * DESCRIPTION: Wavefront alignment algorithm for pairwise gap-affine
 *   alignment (Main module)
 */

#ifndef AFFINE_WAVEFRONT_H_
#define AFFINE_WAVEFRONT_H_


#include "utils/commons.h"
#include "system/profiler_counter.h"
#include "system/profiler_timer.h"

#include "gap_affine/affine_table.h"
#include "gap_affine/affine_wavefront_penalties.h"
#include "gap_affine/affine_wavefront_reduction.h"
#include "gap_affine/wavefront_stats.h"

/*
 * Debug
 */
//#define AFFINE_WAVEFRONT_DEBUG
//#define AFFINE_WAVEFRONT_STATS

/*
 * Constants
 */
#define AFFINE_WAVEFRONT_OFFSET_NULL (-10)
#define AFFINE_WAVEFRONT_K_NULL      (INT_MAX/2)

/*
 * Translate k and offset to coordinates h,v
 */
#define AFFINE_WAVEFRONT_V(k,offset) ((offset)-(k))
#define AFFINE_WAVEFRONT_H(k,offset) (offset)

#define AFFINE_WAVEFRONT_DIAGONAL(h,v) ((h)-(v))
#define AFFINE_WAVEFRONT_OFFSET(h,v)   (h)

/*
 * Offset size
 */
//#define AFFINE_WAVEFRONT_W8
//#define AFFINE_WAVEFRONT_W16
#define AFFINE_WAVEFRONT_W32

#ifdef AFFINE_WAVEFRONT_W8
  typedef int8_t awf_offset_t;
#else
  #ifdef AFFINE_WAVEFRONT_W16
    typedef int16_t awf_offset_t;
  #else // AFFINE_WAVEFRONT_W32
    typedef int32_t awf_offset_t;
  #endif
#endif

/*
 * Wavefront
 */
typedef struct {
  // Range
  bool null;                  // Is null interval?
  int lo;                     // Effective lowest diagonal (inclusive)
  int hi;                     // Effective highest diagonal (inclusive)
  int lo_base;                // Lowest diagonal before reduction (inclusive)
  int hi_base;                // Highest diagonal before reduction (inclusive)
  // Offsets
  awf_offset_t* offsets;      // Offsets
#ifdef AFFINE_WAVEFRONT_DEBUG
  awf_offset_t* offsets_base; // Offsets increment
#endif
} affine_wavefront_t;

/*
 * Gap-Affine Wavefronts
 */
typedef struct {
  // Dimensions
  int pattern_length;                          // Pattern length
  int text_length;                             // Text length
  int num_wavefronts;                          // Total number of allocatable wavefronts
  // Limits
  int max_penalty;                             // MAX(mismatch_penalty,single_gap_penalty)
  int max_k;                                   // Maximum diagonal k (used for null-wf, display, and banding)
  int min_k;                                   // Maximum diagonal k (used for null-wf, display, and banding)
  // Wavefronts
  affine_wavefront_t** mwavefronts;            // M-wavefronts
  affine_wavefront_t** iwavefronts;            // I-wavefronts
  affine_wavefront_t** dwavefronts;            // D-wavefronts
  affine_wavefront_t wavefront_null;           // Null wavefront (used to gain orthogonality)
  // Reduction
  affine_wavefronts_reduction_t reduction;     // Reduction parameters
  // Penalties
  affine_wavefronts_penalties_t penalties;     // Penalties parameters
  // CIGAR
  edit_cigar_t edit_cigar;                     // Alignment CIGAR
  // MM
  mm_allocator_t* mm_allocator;                // MM-Allocator
  affine_wavefront_t* wavefronts_mem;          // MM-Slab for affine_wavefront_t (base)
  affine_wavefront_t* wavefronts_current;      // MM-Slab for affine_wavefront_t (next)
  // STATS
  wavefronts_stats_t* wavefronts_stats; // Stats
  // DEBUG
#ifdef AFFINE_WAVEFRONT_DEBUG
  affine_table_t gap_affine_table;             // DP-Table encoded by the wavefronts
#endif
} affine_wavefronts_t;

/*
 * SWF Wavefront Computation Set
 */
typedef struct {
  /* In Wavefronts*/
  affine_wavefront_t* in_mwavefront_sub;
  affine_wavefront_t* in_mwavefront_gap;
  affine_wavefront_t* in_iwavefront_ext;
  affine_wavefront_t* in_dwavefront_ext;
  /* Out Wavefronts */
  affine_wavefront_t* out_mwavefront;
  affine_wavefront_t* out_iwavefront;
  affine_wavefront_t* out_dwavefront;
} affine_wavefront_set;

/*
 * Setup
 */
void affine_wavefronts_clear(
    affine_wavefronts_t* const affine_wavefronts);
void affine_wavefronts_delete(
    affine_wavefronts_t* const affine_wavefronts);

/*
 * Setup WF-modes
 */
affine_wavefronts_t* affine_wavefronts_new_complete(
    const int pattern_length,
    const int text_length,
    affine_penalties_t* const penalties,
    wavefronts_stats_t* const wavefronts_stats,
    mm_allocator_t* const mm_allocator);
affine_wavefronts_t* affine_wavefronts_new_reduced(
    const int pattern_length,
    const int text_length,
    affine_penalties_t* const penalties,
    const int min_wavefront_length,
    const int max_distance_threshold,
    wavefronts_stats_t* const wavefronts_stats,
    mm_allocator_t* const mm_allocator);

/*
 * Allocate individual wavefront
 */
affine_wavefront_t* affine_wavefronts_allocate_wavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int lo_base,
    const int hi_base);

#endif /* AFFINE_WAVEFRONT_H_ */
