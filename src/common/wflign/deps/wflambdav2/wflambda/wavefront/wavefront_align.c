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
 * DESCRIPTION: WaveFront alignment module for sequence pairwise alignment
 */

#include "wflambda/utils/string_padded.h"
#include "wflambda/wavefront/wavefront_align.h"
#include "wflambda/wavefront/wavefront_extend.h"
#include "wflambda/wavefront/wavefront_compute.h"
#include "wflambda/wavefront/wavefront_compute_affine.h"
#include "wflambda/wavefront/wavefront_compute_affine2p.h"
#include "wflambda/wavefront/wavefront_backtrace.h"

#ifdef WFLAMBDA_NAMESPACE
namespace wflambda {
#endif

/*
 * Alignment end reached
 */
bool wavefront_align_global_terminate(
    wavefront_aligner_t* const wf_aligner,
    const std::function<bool(const int&, const int&)>& traceback_lambda,
    const int pattern_length,
    const int text_length,
    const int score_final) {
  // Parameters
  const int alignment_k = WAVEFRONT_DIAGONAL(text_length,pattern_length);
  const int alignment_offset = WAVEFRONT_OFFSET(text_length,pattern_length);
  int score = score_final;
  // Check wavefront
  if (wf_aligner->memory_modular) score = score % wf_aligner->max_score_scope;
  wavefront_t* const mwavefront = wf_aligner->mwavefronts[score];
  if (mwavefront==NULL) return false;
  // Check limits
  wf_offset_t* const offsets = mwavefront->offsets;
  if (mwavefront->lo > alignment_k || alignment_k > mwavefront->hi) return false;
  // Check offset
  const wf_offset_t offset = offsets[alignment_k];
  if (offset < alignment_offset) return false; // Global termination condition
  // Retrieve alignment
  if (wf_aligner->alignment_scope == alignment_scope_score) {
    wf_aligner->cigar.begin_offset = 0;
    wf_aligner->cigar.end_offset = 0;
    wf_aligner->cigar.score = -score_final;
  } else {
    if (wf_aligner->bt_piggyback) {
      // Fetch backtrace from buffer and recover alignment
      //ToDo
      wf_backtrace_buffer_recover_cigar(
          wf_aligner->bt_buffer,
          mwavefront->bt_pcigar[alignment_k],
          mwavefront->bt_prev[alignment_k],
          traceback_lambda,
          pattern_length,text_length,&wf_aligner->cigar);
    } else {
      // Backtrace alignment
      wavefront_backtrace_affine(wf_aligner,
                                 traceback_lambda,
                                 pattern_length,text_length,
                                 score);
    }
  }
  // Terminate
  return true;
}
/*
 * Initial Conditions
 */
void wavefront_align_global_initialize(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->distance_metric;
  // Init wavefronts
  wf_aligner->mwavefronts[0] = wavefront_slab_allocate(wf_aligner->wavefront_slab,0,0);
  wf_aligner->mwavefronts[0]->offsets[0] = 0;
  if (wf_aligner->bt_piggyback) {
    wf_aligner->mwavefronts[0]->bt_pcigar[0] = 0;
    wf_aligner->mwavefronts[0]->bt_prev[0] = 0; // TODO Put to MAX
  }
  if (distance_metric==edit || distance_metric==gap_lineal) return;
  wf_aligner->d1wavefronts[0] = NULL;
  wf_aligner->i1wavefronts[0] = NULL;
  if (distance_metric==gap_affine) return;
  wf_aligner->d2wavefronts[0] = NULL;
  wf_aligner->i2wavefronts[0] = NULL;
}
/*
 * Global Alignment
 */
int wavefront_align_global(
    wavefront_aligner_t* const wf_aligner,
    const std::function<bool(const int&, const int&)>& match_lambda,
    const std::function<bool(const int&, const int&)>& traceback_lambda,
    const int pattern_length,
    const int text_length) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->distance_metric;
  // Initialize wavefront
  wavefront_align_global_initialize(wf_aligner);
  // Compute wavefronts of increasing score
  int score = 0;
  while (true) {
    // Exact extend s-wavefront
    wavefront_extend(wf_aligner, match_lambda,
                     pattern_length, text_length,
                     score);
    // Exit condition
    if (wavefront_align_global_terminate(
        wf_aligner,traceback_lambda,pattern_length,text_length,score)) break;
    // Compute (s+1)-wavefront
    ++score;
    switch (distance_metric) {
      case gap_affine:
        wavefront_compute_affine(wf_aligner,
                                 pattern_length, text_length,
                                 score);
        break;
      case gap_affine_2p:
        wavefront_compute_affine2p(wf_aligner,
                                   pattern_length,text_length,
                                   score);
        break;
      default:
        fprintf(stderr,"Distance function not yet implemented\n"); exit(1);
        break;
    }
    // DEBUG
    //wavefront_aligner_print(stderr,wf_aligner,score,score,2,16);
  }
  return score;
}
/*
 * Bounded Global Alignment
 */
int wavefront_align_global_bounded(
    wavefront_aligner_t* const wf_aligner,
    const std::function<bool(const int&, const int&)>& match_lambda,
    const std::function<bool(const int&, const int&)>& traceback_lambda,
    const int pattern_length,
    const int text_length,
    const int max_score) {
  // Parameters
  const distance_metric_t distance_metric = wf_aligner->distance_metric;
  // Initialize wavefront
  wavefront_align_global_initialize(wf_aligner);
  // Compute wavefronts of increasing score
  int score = 0;
  while (true) {
    // Exact extend s-wavefront
    wavefront_extend(wf_aligner,match_lambda,
                     pattern_length,text_length,
                     score);
    // Exit condition
    if (wavefront_align_global_terminate(
        wf_aligner,traceback_lambda,pattern_length,text_length,score)) break;
    // Compute (s+1)-wavefront
    ++score;
    // Bound the alignment at our max_score
    if (score > max_score) {
        break; // todo... signal failure somehow
    }
    switch (distance_metric) {
      case gap_affine:
        wavefront_compute_affine(wf_aligner,
                                 pattern_length,text_length,
                                 score);
        break;
      case gap_affine_2p:
        wavefront_compute_affine2p(wf_aligner,
                                   pattern_length,text_length,
                                   score);
        break;
      default:
        fprintf(stderr,"Distance function not yet implemented\n"); exit(1);
        break;
    }
    // DEBUG
    //wavefront_aligner_print(stderr,wf_aligner,score,score,2,16);
  }
  return score;
}
/*
 * Semi-Global Alignment
 */
void wavefront_align_semiglobal(
    wavefront_aligner_t* const wf_aligner,
    strings_padded_t* const sequences,
    const int pattern_length,
    const int text_length) {
  // TODO Incorporate HERE the end-condition
}
/*
 * Wavefront Alignment
 */
int wavefront_align(
    wavefront_aligner_t* const wf_aligner,
    const std::function<bool(const int&, const int&)>& match_lambda,
    const std::function<bool(const int&, const int&)>& traceback_lambda,
    const int pattern_length,
    const int text_length) {

  // Alignment computing wavefronts
  int score = wavefront_align_global(wf_aligner,match_lambda,traceback_lambda,pattern_length,text_length);

  // Return our resulting score
  return score;
}
/*
 * Bounded Wavefront Alignment
 */
int wavefront_align_bounded(
    wavefront_aligner_t* const wf_aligner,
    const std::function<bool(const int&, const int&)>& match_lambda,
    const std::function<bool(const int&, const int&)>& traceback_lambda,
    const int pattern_length,
    const int text_length,
    const int max_score) {

  // Alignment computing wavefronts
  int score = wavefront_align_global_bounded(wf_aligner,match_lambda,traceback_lambda,pattern_length,text_length,max_score);

  // Return our resulting score
  return score;
}

#ifdef WFLAMBDA_NAMESPACE
}
#endif
