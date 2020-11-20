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
 * DESCRIPTION: Module to track relevant stats from the WFA
 */

#ifndef WAVEFRONT_STATS_H_
#define WAVEFRONT_STATS_H_

#include "utils/commons.hpp"
#include "system/profiler_counter.hpp"
#include "system/profiler_timer.hpp"

namespace wflambda {

/*
 * Stats
 */
#ifdef AFFINE_LAMBDA_WAVEFRONT_STATS
  #define WAVEFRONT_STATS_COUNTER_ADD(wavefronts,counter,amount) \
    if (wavefronts->wf_stats!=NULL) counter_add(&(wavefronts->wavefronts_stats->counter),(amount))
  #define WAVEFRONT_STATS_TIMER_START(wavefronts,timer) \
    if (wavefronts->wf_stats!=NULL) timer_start(&(wavefronts->wavefronts_stats->timer))
  #define WAVEFRONT_STATS_TIMER_STOP(wavefronts,timer) \
    if (wavefronts->wf_stats!=NULL) timer_stop(&(wavefronts->wavefronts_stats->timer))
#else
  #define WAVEFRONT_STATS_COUNTER_ADD(wf,counter,amount)
  #define WAVEFRONT_STATS_TIMER_START(wf,timer)
  #define WAVEFRONT_STATS_TIMER_STOP(wf,timer)
#endif

/*
 * Wavefront Stats
 */
typedef struct {
  profiler_counter_t wf_score;              // Score reached by WF-alignment
  profiler_counter_t wf_steps;              // Step performed by WF-alignment
  profiler_counter_t wf_steps_null;         // Avoided WF-alignment steps due to null wavefronts
  profiler_counter_t wf_steps_extra;        // Extra steps performed by WF-alignment searching for a better solution
  profiler_counter_t wf_operations;         // Single cell WF-operations performed
  profiler_counter_t wf_extensions;         // Single cell WF-extensions performed
  profiler_counter_t wf_reduction;          // Calls to reduce wavefront
  profiler_counter_t wf_reduced_cells;      // Total cells reduced
  profiler_counter_t wf_null_used;          // Total times a null-WF was used as padding
  profiler_counter_t wf_extend_inner_loop;  // Total times SIMD-extension had to re-iterate
  profiler_counter_t wf_compute_kernel[4];  // Specialized WF computation kernel used
  profiler_timer_t wf_time_backtrace;       // Time spent doing backtrace
  profiler_counter_t wf_backtrace_paths;    // Total paths explored by the backtrace
  profiler_counter_t wf_backtrace_alg;      // Total alignments explored by the backtrace
} wavefronts_stats_t;

/*
 * Setup
 */
void wavefronts_stats_clear(wavefronts_stats_t* const wavefronts_stats);

/*
 * Display
 */
void wavefronts_stats_print(
    FILE* const stream,
    wavefronts_stats_t* const wavefronts_stats);

}

#endif /* WAVEFRONT_STATS_H_ */
