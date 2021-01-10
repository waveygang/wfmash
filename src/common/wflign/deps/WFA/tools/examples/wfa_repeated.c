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
 * DESCRIPTION: WFA Sample-Code
 */

#include "gap_affine/affine_wavefront_align.h"

int main(int argc,char* argv[]) {
  // Patter & Text
  char* pattern = "TCTTTACTCGCGCGTTTCTTACTCGCGCGTTGGAGAAATACAATAGTGGAGAAATACAATAGTTTTTTTTTTTT";
  char* text    = "TTTTTTCTATACTGCGCGTTTTCTATACTCGCGCGTTGGAGAAATACAATAGTGGAGAAATAAAATAGT";
  // Allocate MM
  mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  // Set penalties
  affine_penalties_t affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  };
  // Init Affine-WFA
  affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(
      strlen(pattern),strlen(text),&affine_penalties,NULL,mm_allocator);
  // Repeat alignment (for the sake of it)
  int i;
  for (i=0;i<100000;++i) {
    // Clean
    affine_wavefronts_clear(affine_wavefronts);
    // Align
    affine_wavefronts_align(
        affine_wavefronts,pattern,strlen(pattern),text,strlen(text));
    // Report
    if ((i%1000) == 0) {
      const int score = edit_cigar_score_gap_affine(
          &affine_wavefronts->edit_cigar,&affine_penalties);
      fprintf(stderr,"Done %d, score is %d\n",i,score);
    }
  }
  // Free
  affine_wavefronts_delete(affine_wavefronts);
  mm_allocator_delete(mm_allocator);
}
