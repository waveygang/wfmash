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
 * DESCRIPTION: WFA display functions
 */

#include "gap_affine/affine_wavefront_display.h"
#include "gap_affine/affine_wavefront_utils.h"

/*
 * Accessors
 */
void affine_wavefronts_set_edit_table(
    affine_wavefronts_t* const affine_wavefronts,
    const int pattern_length,
    const int text_length,
    const int k,
    const awf_offset_t offset,
    const int score) {
#ifdef AFFINE_WAVEFRONT_DEBUG
  if (offset < 0) return;
  const int h = AFFINE_WAVEFRONT_H(k,offset);
  const int v = AFFINE_WAVEFRONT_V(k,offset);
  if (0 <= v && v <= pattern_length &&
      0 <= h && h <= text_length) {
    if (affine_wavefronts->gap_affine_table.columns[h][v].M == -1) {
      affine_wavefronts->gap_affine_table.columns[h][v].M = score;
    }
//    // Overwritte
//    else {
//      fprintf(stderr,"Attempt to overwrite (score %d over %d)\n",
//          score,affine_wavefronts->gap_affine_table.columns[h][v].M);
//    }
//    // Out negative
  }
//  else if (v < 0 || h < 0) {
//    fprintf(stderr,"Out-of-table. Negative values (h,v)=(%d,%d)\n",h,v);
//  }
//    // Out
//  else {
//    fprintf(stderr,"Out-of-table. Beyond table limits\n");
//  }
#endif
}
/*
 * Display
 */
#define AFFINE_WAVEFRONTS_PRINT_ELEMENT(wavefront,k) \
  /* Check limits */ \
  if (wavefront!=NULL && wavefront->lo <= k && k <= wavefront->hi) { \
    if (wavefront->offsets[k] >= 0) { \
      fprintf(stream,"[%2d]",(int)wavefront->offsets[k]); \
    } else { \
      fprintf(stream,"[  ]"); \
    } \
  } else { \
    fprintf(stream,"    "); \
  }
void affine_wavefronts_print_wavefront(
    FILE* const stream,
    affine_wavefronts_t* const affine_wavefronts,
    const int current_score) {
  // Print header
  fprintf(stream,">[SCORE=%3d]\n",current_score);
  if (affine_wavefronts->iwavefronts != NULL) {
    fprintf(stream," [M] [I] [D] ");
  } else {
    fprintf(stream," [M] ");
  }
  fprintf(stream,"\n");
  // Compute min/max k
  int k, max_k=0, min_k=0;
  for (k=affine_wavefronts->max_k;k>=affine_wavefronts->min_k;k--) {
    // Fetch wavefront
    affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[current_score];
    if (mwavefront != NULL) {
      max_k = MAX(max_k,mwavefront->hi);
      min_k = MIN(min_k,mwavefront->lo);
    }
    if (affine_wavefronts->iwavefronts != NULL) {
      affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[current_score];
      affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[current_score];
      if (iwavefront != NULL) {
        max_k = MAX(max_k,iwavefront->hi);
        min_k = MIN(min_k,iwavefront->lo);
      }
      if (dwavefront != NULL) {
        max_k = MAX(max_k,dwavefront->hi);
        min_k = MIN(min_k,dwavefront->lo);
      }
    }
  }
  // Print wavefront values
  for (k=max_k;k>=min_k;k--) {
    fprintf(stream,"[k=%3d] ",k);
    // Fetch wavefront
    affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[current_score];
    if (affine_wavefronts->iwavefronts != NULL) {
      affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[current_score];
      affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[current_score];
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(mwavefront,k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(iwavefront,k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(dwavefront,k);
    } else {
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(mwavefront,k);
    }
    fprintf(stream,"\n");
  }
}
void affine_wavefronts_print_wavefronts(
    FILE* const stream,
    affine_wavefronts_t* const affine_wavefronts,
    const int current_score) {
  // Headers
  int k, s;
  fprintf(stream,">[SCORE=%3d]\n",current_score);
  fprintf(stream,"        ");
  for (s=0;s<=current_score;++s) {
    if (affine_wavefronts->iwavefronts != NULL) {
      fprintf(stream," [M] [I] [D] ");
    } else {
      fprintf(stream," [M] ");
    }
  }
  fprintf(stream,"\n");
  // Compute min/max k
  int max_k=0, min_k=0;
  for (k=affine_wavefronts->max_k;k>=affine_wavefronts->min_k;k--) {
    for (s=0;s<=current_score;++s) {
      // Fetch wavefront
      affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[s];
      if (mwavefront != NULL) {
        max_k = MAX(max_k,mwavefront->hi);
        min_k = MIN(min_k,mwavefront->lo);
      }
      affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[s];
      if (iwavefront != NULL) {
        max_k = MAX(max_k,iwavefront->hi);
        min_k = MIN(min_k,iwavefront->lo);
      }
      affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[s];
      if (dwavefront != NULL) {
        max_k = MAX(max_k,dwavefront->hi);
        min_k = MIN(min_k,dwavefront->lo);
      }
    }
  }
  // Traverse all diagonals
  for (k=max_k;k>=min_k;k--) {
    fprintf(stream,"[k=%3d] ",k);
    // Traverse all scores
    for (s=0;s<=current_score;++s) {
      // Fetch wavefront
      affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[s];
      affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[s];
      affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[s];
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(mwavefront,k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(iwavefront,k);
      AFFINE_WAVEFRONTS_PRINT_ELEMENT(dwavefront,k);
      fprintf(stream," ");
    }
    fprintf(stream,"\n");
  }
  // Headers
  fprintf(stream,"SCORE   ");
  for (s=0;s<=current_score;++s) {
    fprintf(stream,"     %2d      ",s);
  }
  fprintf(stream,"\n");
}
void affine_wavefronts_print_wavefronts_pretty(
    FILE* const stream,
    affine_wavefronts_t* const affine_wavefronts,
    const int current_score) {
  // Header
  fprintf(stream,">[SCORE=%3d]\n",current_score);
  int k, s;
  awf_offset_t offset;
  // Traverse all diagonals
  for (k=affine_wavefronts->max_k;k>=affine_wavefronts->min_k;k--) {
    // Header
    fprintf(stream,"[k=%3d] ",k);
    // Padding
    awf_offset_t max_offset = affine_wavefronts_diagonal_length(affine_wavefronts,k);
    if (max_offset <= 0) {
      fprintf(stream,"\n");
      continue;
    }
    offset = 0; // No table frame (initial conditions)
    if (k > 0) {
      max_offset += k;
      for (;offset<k;++offset) fprintf(stream,"   ");
    }
    // Fetch offsets from 0 - score
    int last_offset = -1;
    for (s=0;s<=current_score;++s) {
      // Fetch wavefront
      affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[s];
      if (mwavefront == NULL) continue;
      if (k < mwavefront->lo || mwavefront->hi < k) continue;
      // Fetch offset
      awf_offset_t offsets_sub = mwavefront->offsets[k];
      awf_offset_t offsets_sub_base = offsets_sub;
#ifdef AFFINE_WAVEFRONT_DEBUG
      offsets_sub_base = mwavefront->offsets_base[k];
#endif
      if (offsets_sub < 0) continue;
      // Print offsets
      if (last_offset == -1) {
        while (offset < offsets_sub_base-1) {
          if (offset == max_offset) fprintf(stream," | ");
          fprintf(stream," * "); ++offset;
        }
      }
      while (offset < offsets_sub) {
        if (offset == max_offset) fprintf(stream," | ");
        fprintf(stream,"%2d ",s); ++offset;
      }
      // Update last
      last_offset = offsets_sub;
    }
    // Print unknowns
    for (;offset<=max_offset;++offset) {
      if (offset == max_offset) fprintf(stream," | ");
      fprintf(stream," * ");
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");
}
/*
 * Debug
 */
void affine_wavefronts_debug_step(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const char* const text,
    const int score) {
#ifdef AFFINE_WAVEFRONT_DEBUG
//  affine_table_print(stderr,&affine_wavefronts->gap_affine_table,pattern,text);
//  affine_wavefronts_print_wavefronts_pretty(stderr,affine_wavefronts,score);
//  affine_wavefronts_print_wavefronts(stderr,affine_wavefronts,score);
//  affine_wavefronts_print_wavefront(stderr,affine_wavefronts,score);
//  affine_wavefronts_print_computation_stats(stderr,affine_wavefronts);
#endif
}


