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
 * DESCRIPTION: WaveFront-Alignment module for the "extension" of exact matches
 */

#include <functional>
#include "wflambda/utils/string_padded.h"
#include "wflambda/wavefront/wavefront_extend.h"
#include "wflambda/wavefront/wavefront_compute.h"
#include "wflambda/wavefront/wavefront_reduction.h"

#ifdef WFLAMBDA_NAMESPACE
namespace wflambda {
#endif

/*
 * Wavefront offset extension comparing characters
 */
void wavefront_extend_packed(
    wavefront_aligner_t* const wf_aligner,
    const std::function<bool(const int&, const int&)>& match_lambda,
    const int pattern_length,
    const int text_length,
    const int score) {
  // Fetch m-wavefront
  wavefront_t* const mwavefront = wf_aligner->mwavefronts[score];
  if (mwavefront==NULL) return;
  // Extend diagonally each wavefront point
  wf_offset_t* const offsets = mwavefront->offsets;
  int k;
  for (k=mwavefront->lo;k<=mwavefront->hi;++k) {
    // Fetch offset & positions
    wf_offset_t offset = offsets[k];
    uint32_t h = WAVEFRONT_H(k,offset); // Make unsigned to avoid checking negative
    uint32_t v = WAVEFRONT_V(k,offset); // Make unsigned to avoid checking negative

    while (h < text_length&& v < pattern_length && match_lambda(v++,h++)) {
        ++(offsets[k]);
    }
  }
}
/*
 * Wavefront exact "extension"
 */
void wavefront_extend(
    wavefront_aligner_t* const wf_aligner,
    const std::function<bool(const int&, const int&)>& match_lambda,
    const int pattern_length,
    const int text_length,
    int score) {
  // Modular wavefront
  if (wf_aligner->memory_modular) score = score % wf_aligner->max_score_scope;
  // Extend wavefront
  wavefront_extend_packed(
      wf_aligner,match_lambda,pattern_length,text_length,score);
  // Reduce wavefront dynamically
  if (wf_aligner->reduction.reduction_strategy == wavefront_reduction_dynamic) {
    wavefront_reduce(wf_aligner,pattern_length,text_length,score);
  }
}

#ifdef WFLAMBDA_NAMESPACE
}
#endif
