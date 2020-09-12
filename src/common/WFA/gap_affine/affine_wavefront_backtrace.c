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
 * DESCRIPTION: WFA extend backtrace component
 */

#include "gap_affine/affine_wavefront_backtrace.h"

/*
 * Backtrace Detect Limits
 */
bool affine_wavefronts_valid_location(
    const int k,
    const awf_offset_t offset,
    const int pattern_length,
    const int text_length) {
  // Locate offset (remember that backtrace is always +1 offset ahead)
  const int v = AFFINE_WAVEFRONT_V(k,offset);
  const int h = AFFINE_WAVEFRONT_H(k,offset);
  return (v > 0 && v <= pattern_length &&
          h > 0 && h <= text_length);
}
void affine_wavefronts_offset_add_trailing_gap(
    edit_cigar_t* const edit_cigar,
    const int k,
    const int alignment_k) {
  // Parameters
  char* const operations = edit_cigar->operations;
  int op_sentinel = edit_cigar->begin_offset;
  // Add trailing gap
  int i;
  if (k < alignment_k) {
    for (i=k;i<alignment_k;++i) operations[op_sentinel--] = 'I';
  } else if (k > alignment_k) {
    for (i=alignment_k;i<k;++i) operations[op_sentinel--] = 'D';
  }
  edit_cigar->begin_offset = op_sentinel;
}
/*
 * Backtrace Paths Offsets
 */
awf_offset_t backtrace_wavefront_trace_deletion_open_offset(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k+1 &&
      k+1 <= mwavefront->hi_base) {
    return mwavefront->offsets[k+1];
  } else {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_deletion_extend_offset(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[score];
  if (dwavefront != NULL &&
      dwavefront->lo_base <= k+1 &&
      k+1 <= dwavefront->hi_base) {
    return dwavefront->offsets[k+1];
  } else {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_insertion_open_offset(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k-1 &&
      k-1 <= mwavefront->hi_base) {
    return mwavefront->offsets[k-1] + 1;
  } else {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_insertion_extend_offset(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[score];
  if (iwavefront != NULL &&
      iwavefront->lo_base <= k-1 &&
      k-1 <= iwavefront->hi_base) {
    return iwavefront->offsets[k-1] + 1;
  } else {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
awf_offset_t backtrace_wavefront_trace_mismatch_offset(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return AFFINE_WAVEFRONT_OFFSET_NULL;
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  if (mwavefront != NULL &&
      mwavefront->lo_base <= k &&
      k <= mwavefront->hi_base) {
    return mwavefront->offsets[k] + 1;
  } else {
    return AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
/*
 * Backtrace Paths Conditions
 */
bool backtrace_wavefront_trace_deletion(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return false;
  affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[score];
  return (dwavefront != NULL &&
          dwavefront->lo_base <= k &&
          k <= dwavefront->hi_base &&
          offset == dwavefront->offsets[k]);
}
bool backtrace_wavefront_trace_deletion_open(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return false;
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  return (mwavefront != NULL &&
          mwavefront->lo_base <= k+1 &&
          k+1 <= mwavefront->hi_base &&
          offset == mwavefront->offsets[k+1]);
}
bool backtrace_wavefront_trace_deletion_extend(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return false;
  affine_wavefront_t* const dwavefront = affine_wavefronts->dwavefronts[score];
  return (dwavefront != NULL &&
          dwavefront->lo_base <= k+1 &&
          k+1 <= dwavefront->hi_base &&
          offset == dwavefront->offsets[k+1]);
}
bool backtrace_wavefront_trace_insertion(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return false;
  affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[score];
  return (iwavefront != NULL &&
          iwavefront->lo_base <= k &&
          k <= iwavefront->hi_base &&
          offset == iwavefront->offsets[k]);
}
bool backtrace_wavefront_trace_insertion_open(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return false;
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  return (mwavefront != NULL &&
          mwavefront->lo_base <= k-1 &&
          k-1 <= mwavefront->hi_base &&
          offset == mwavefront->offsets[k-1]+1);
}
bool backtrace_wavefront_trace_insertion_extend(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return false;
  affine_wavefront_t* const iwavefront = affine_wavefronts->iwavefronts[score];
  return (iwavefront != NULL &&
          iwavefront->lo_base <= k-1 &&
          k-1 <= iwavefront->hi_base &&
          offset == iwavefront->offsets[k-1]+1);
}
bool backtrace_wavefront_trace_mismatch(
    affine_wavefronts_t* const affine_wavefronts,
    const int score,
    const int k,
    const awf_offset_t offset) {
  if (score < 0) return false;
  affine_wavefront_t* const mwavefront = affine_wavefronts->mwavefronts[score];
  return (mwavefront != NULL &&
          mwavefront->lo_base <= k &&
          k <= mwavefront->hi_base &&
          offset == mwavefront->offsets[k]+1);
}
/*
 * Backtrace Operations
 */
int affine_wavefronts_backtrace_compute_max_matches(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const char* const text,
    const int k,
    awf_offset_t offset) {
  // Locate position
  int v = AFFINE_WAVEFRONT_V(k,offset);
  int h = AFFINE_WAVEFRONT_H(k,offset);
  // Check matches
  int num_matches = 0;
  while (v>0 && h>0 && pattern[--v] == text[--h]) ++num_matches;
  // Return max-matches
  return num_matches;
}
void affine_wavefronts_backtrace_matches__check(
    affine_wavefronts_t* const affine_wavefronts,
    const char* const pattern,
    const char* const text,
    const int k,
    awf_offset_t offset,
    const bool valid_location,
    const int num_matches,
    edit_cigar_t* const edit_cigar) {
  int i;
  for (i=0;i<num_matches;++i) {
    // DEBUG
#ifdef AFFINE_WAVEFRONT_DEBUG
    const int v = AFFINE_WAVEFRONT_V(k,offset);
    const int h = AFFINE_WAVEFRONT_H(k,offset);
    if (!valid_location) { // Check inside table
      fprintf(stderr,"Backtrace error: Match outside DP-Table\n");
      exit(1);
    } else if (pattern[v-1] != text[h-1]) { // Check match
      fprintf(stderr,"Backtrace error: Not a match traceback\n");
      exit(1);
    }
#endif
    // Set Match
    edit_cigar->operations[(edit_cigar->begin_offset)--] = 'M';
    // Update state
    --offset;
  }
}
void affine_wavefronts_backtrace_matches(
    edit_cigar_t* const edit_cigar,
    const int num_matches) {
  // Set Matches
  int i;
  for (i=0;i<num_matches;++i) {
    edit_cigar->operations[(edit_cigar->begin_offset)--] = 'M';
  }
}
/*
 * Backtrace (single solution)
 */
void affine_wavefronts_backtrace(
    affine_wavefronts_t* const affine_wavefronts,
    char* const pattern,
    const int pattern_length,
    char* const text,
    const int text_length,
    const int alignment_score) {
  // STATS
  WAVEFRONT_STATS_TIMER_START(affine_wavefronts,wf_time_backtrace);
  // Parameters
  const affine_penalties_t* const wavefront_penalties =
      &(affine_wavefronts->penalties.wavefront_penalties);
  edit_cigar_t* const cigar = &affine_wavefronts->edit_cigar;
  const int alignment_k = AFFINE_WAVEFRONT_DIAGONAL(text_length,pattern_length);
  // Compute starting location
  int score = alignment_score;
  int k = alignment_k;
  awf_offset_t offset = affine_wavefronts->mwavefronts[alignment_score]->offsets[k];
  bool valid_location = affine_wavefronts_valid_location(k,offset,pattern_length,text_length);
  // Trace the alignment back
  backtrace_wavefront_type backtrace_type = backtrace_wavefront_M;
  int v = AFFINE_WAVEFRONT_V(k,offset);
  int h = AFFINE_WAVEFRONT_H(k,offset);
  while (v > 0 && h > 0 && score > 0) {
    // Check location
    if (!valid_location) {
      valid_location = affine_wavefronts_valid_location(k,offset,pattern_length,text_length);
      if (valid_location) {
        affine_wavefronts_offset_add_trailing_gap(cigar,k,alignment_k);
      }
    }
    // Compute scores
    const int gap_open_score = score - wavefront_penalties->gap_opening - wavefront_penalties->gap_extension;
    const int gap_extend_score = score - wavefront_penalties->gap_extension;
    const int mismatch_score = score - wavefront_penalties->mismatch;
    // Compute source offsets
    const awf_offset_t del_ext = (backtrace_type == backtrace_wavefront_I) ? AFFINE_WAVEFRONT_OFFSET_NULL:
        backtrace_wavefront_trace_deletion_extend_offset(affine_wavefronts,gap_extend_score,k,offset);
    const awf_offset_t del_open = (backtrace_type == backtrace_wavefront_I) ? AFFINE_WAVEFRONT_OFFSET_NULL:
        backtrace_wavefront_trace_deletion_open_offset(affine_wavefronts,gap_open_score,k,offset);
    const awf_offset_t ins_ext = (backtrace_type == backtrace_wavefront_D) ? AFFINE_WAVEFRONT_OFFSET_NULL:
        backtrace_wavefront_trace_insertion_extend_offset(affine_wavefronts,gap_extend_score,k,offset);
    const awf_offset_t ins_open = (backtrace_type == backtrace_wavefront_D) ? AFFINE_WAVEFRONT_OFFSET_NULL:
        backtrace_wavefront_trace_insertion_open_offset(affine_wavefronts,gap_open_score,k,offset);
    const awf_offset_t misms = (backtrace_type != backtrace_wavefront_M) ? AFFINE_WAVEFRONT_OFFSET_NULL:
        backtrace_wavefront_trace_mismatch_offset(affine_wavefronts,mismatch_score,k,offset);
    // Compute maximum offset
    const awf_offset_t max_del = MAX(del_ext,del_open);
    const awf_offset_t max_ins = MAX(ins_ext,ins_open);
    const awf_offset_t max_all = MAX(misms,MAX(max_ins,max_del));
    // Traceback Matches
    if (backtrace_type == backtrace_wavefront_M) {
      const int num_matches = offset - max_all;
      affine_wavefronts_backtrace_matches__check(affine_wavefronts,
          pattern,text,k,offset,valid_location,num_matches,cigar);
      offset = max_all;
    }
    // Traceback Operation
    if (max_all == del_ext) {
      // Add Deletion
      if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'D';
      // Update state
      score = gap_extend_score;
      ++k;
      backtrace_type = backtrace_wavefront_D;
    } else if (max_all == del_open) {
      // Add Deletion
      if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'D';
      // Update state
      score = gap_open_score;
      ++k;
      backtrace_type = backtrace_wavefront_M;
    } else if (max_all == ins_ext) {
      // Add Insertion
      if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'I';
      // Update state
      score = gap_extend_score;
      --k;
      --offset;
      backtrace_type = backtrace_wavefront_I;
    } else if (max_all == ins_open) {
      // Add Insertion
      if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'I';
      // Update state
      score = gap_open_score;
      --k;
      --offset;
      backtrace_type = backtrace_wavefront_M;
    } else if (max_all == misms) {
      // Add Mismatch
      if (valid_location) cigar->operations[(cigar->begin_offset)--] = 'X';
      // Update state
      score = mismatch_score;
      --offset;
    } else {
      fprintf(stderr,"Backtrace error: No link found during backtrace\n");
      exit(1);
    }
    // Update coordinates
    v = AFFINE_WAVEFRONT_V(k,offset);
    h = AFFINE_WAVEFRONT_H(k,offset);
  }
  // Account for last operations
  if (score == 0) {
    // Account for last stroke of matches
    affine_wavefronts_backtrace_matches__check(affine_wavefronts,
        pattern,text,k,offset,valid_location,offset,cigar);
  } else {
    // Account for last stroke of insertion/deletion
    while (v > 0) {cigar->operations[(cigar->begin_offset)--] = 'D'; --v;};
    while (h > 0) {cigar->operations[(cigar->begin_offset)--] = 'I'; --h;};
  }
  ++(cigar->begin_offset); // Set CIGAR length
  // STATS
  WAVEFRONT_STATS_TIMER_STOP(affine_wavefronts,wf_time_backtrace);
}
