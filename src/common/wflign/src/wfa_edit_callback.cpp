/*
 *  Wavefront Alignments Algorithms
 *  Copyright (c) 2019 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of Wavefront Alignments Algorithms.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "wfa_edit_callback.hpp"

namespace wflign {

void edit_wavefronts_init(
    edit_wavefronts_t* const wavefronts,
    const int pattern_length,
    const int text_length) {
  // Dimensions
  wavefronts->pattern_length = pattern_length;
  wavefronts->text_length = text_length;
  wavefronts->max_distance = pattern_length + text_length;
  // Allocate wavefronts
  wavefronts->wavefronts = (edit_wavefront_t*)calloc(wavefronts->max_distance,sizeof(edit_wavefront_t));
  wavefronts->wavefronts_allocated = 0;
  // Allocate CIGAR
  wavefronts->edit_cigar = (char*)malloc(wavefronts->max_distance);
}

void edit_wavefronts_clean(
    edit_wavefronts_t* const wavefronts) {
  int i;
  for (i=0;i<wavefronts->wavefronts_allocated;++i) {
    free(wavefronts->wavefronts[i].offsets_mem);
  }
  wavefronts->wavefronts_allocated = 0;
}

edit_wavefront_t* edit_wavefronts_allocate_wavefront(
    edit_wavefronts_t* const edit_wavefronts,
    const int distance,
    const int lo_base,
    const int hi_base) {
  // Compute limits
  const int wavefront_length = hi_base - lo_base + 2; // (+1) for k=0
  // Allocate wavefront
  edit_wavefront_t* const wavefront = edit_wavefronts->wavefronts + distance;
  ++(edit_wavefronts->wavefronts_allocated); // Next
  // Configure offsets
  wavefront->lo = lo_base;
  wavefront->hi = hi_base;
  // Allocate offsets
  ewf_offset_t* const offsets_mem = (ewf_offset_t*)calloc(wavefront_length,sizeof(ewf_offset_t));
  wavefront->offsets_mem = offsets_mem;
  wavefront->offsets = offsets_mem - lo_base; // Center at k=0
  ++(edit_wavefronts->wavefronts_allocated);
  // Return
  return wavefront;
}

/*
 * Edit Wavefront Backtrace
 */
int edit_wavefronts_backtrace(
    edit_wavefronts_t* const wavefronts,
    const int target_k,
    const int target_distance) {
  // Parameters
  int edit_cigar_idx = 0;
  int k = target_k, distance = target_distance;
  ewf_offset_t offset = wavefronts->wavefronts[distance].offsets[k];
  while (distance > 0) {
    // Fetch
    const edit_wavefront_t* const wavefront = &wavefronts->wavefronts[distance-1];
    const ewf_offset_t* const offsets = wavefront->offsets;
    // Traceback operation
    if (wavefront->lo <= k+1 && k+1 <= wavefront->hi && offset == offsets[k+1]) {
      wavefronts->edit_cigar[edit_cigar_idx++] = 'D';
      ++k;
      --distance;
    } else if (wavefront->lo <= k-1 && k-1 <= wavefront->hi && offset == offsets[k-1] + 1) {
      wavefronts->edit_cigar[edit_cigar_idx++] = 'I';
      --k;
      --offset;
      --distance;
    } else if (wavefront->lo <= k && k <= wavefront->hi && offset == offsets[k] + 1) {
      wavefronts->edit_cigar[edit_cigar_idx++] = 'X';
      --distance;
      --offset;
    } else {
      wavefronts->edit_cigar[edit_cigar_idx++] = 'M';
      --offset;
    }
  }
  // Account for last offset of matches
  while (offset > 0) {
    wavefronts->edit_cigar[edit_cigar_idx++] = 'M';
    --offset;
  }
  //std::cerr << std::endl;
  // Return CIGAR length
  return edit_cigar_idx;
}

/*
 * Extend Wavefront
 */
void edit_wavefronts_extend_wavefront(
    edit_wavefronts_t* const wavefronts,
    const std::function<bool(const int&,const int&)>& extend_match,
    const int pattern_length,
    const int text_length,
    const int distance) {
  // Parameters
  edit_wavefront_t* const wavefront = &wavefronts->wavefronts[distance];
  ewf_offset_t* const offsets = wavefront->offsets;
  const int k_min = wavefront->lo;
  const int k_max = wavefront->hi;
  // Extend diagonally each wavefront point
  int k;
  for (k=k_min;k<=k_max;++k) {
    int v = EWAVEFRONT_V(k,offsets[k]);
    int h = EWAVEFRONT_H(k,offsets[k]);
    while (v<pattern_length && h<text_length && extend_match(v++,h++)) {
        ++(offsets[k]);
        //std::cerr << v-1 << "\t" << h-1 << "\t" << distance << "\t" << "aligned" << std::endl;
    }
    //std::cerr << v << "\t" << h << "\t" << distance << "\t" << "unaligned" << std::endl;
  }
}

/*
 * Edit Wavefront Compute
 */
void edit_wavefronts_compute_wavefront(
    edit_wavefronts_t* const wavefronts,
    const int pattern_length,
    const int text_length,
    const int distance) {
  // Fetch wavefronts
  edit_wavefront_t* const wavefront = &wavefronts->wavefronts[distance-1];
  const int hi = wavefront->hi;
  const int lo = wavefront->lo;
  edit_wavefront_t* const next_wavefront = edit_wavefronts_allocate_wavefront(wavefronts,distance,lo-1,hi+1);
  // Fetch offsets
  ewf_offset_t* const offsets = wavefront->offsets;
  ewf_offset_t* const next_offsets = next_wavefront->offsets;
  // Loop peeling (k=lo-1)
  next_offsets[lo-1] = offsets[lo];
  // Loop peeling (k=lo)
  const ewf_offset_t bottom_upper_del = ((lo+1) <= hi) ? offsets[lo+1] : -1;
  next_offsets[lo] = MAX(offsets[lo]+1,bottom_upper_del);
  // Compute next wavefront starting point
  int k;
  #pragma GCC ivdep
  for (k=lo+1;k<=hi-1;++k) {
    /*
     * const int del = offsets[k+1]; // Upper
     * const int sub = offsets[k] + 1; // Mid
     * const int ins = offsets[k-1] + 1; // Lower
     * next_offsets[k] = MAX(sub,ins,del); // MAX
     */
    const ewf_offset_t max_ins_sub = MAX(offsets[k],offsets[k-1]) + 1;
    next_offsets[k] = MAX(max_ins_sub,offsets[k+1]);
  }
  // Loop peeling (k=hi)
  const ewf_offset_t top_lower_ins = (lo <= (hi-1)) ? offsets[hi-1] : -1;
  next_offsets[hi] = MAX(offsets[hi],top_lower_ins) + 1;
  // Loop peeling (k=hi+1)
  next_offsets[hi+1] = offsets[hi] + 1;
}

/*
 * Edit distance alignment using wavefronts
 */
void edit_wavefronts_align(
    edit_wavefronts_t* const wavefronts,
    const std::function<bool(const int&,const int&)>& extend_match,
    const int pattern_length,
    const int text_length) {
  // Parameters
  const int max_distance = pattern_length + text_length;
  const int target_k = EWAVEFRONT_DIAGONAL(text_length,pattern_length);
  const int target_k_abs = ABS(target_k);
  const ewf_offset_t target_offset = EWAVEFRONT_OFFSET(text_length,pattern_length);
  // Init wavefronts
  int distance;
  edit_wavefronts_allocate_wavefront(wavefronts,0,0,0);
  wavefronts->wavefronts[0].offsets[0] = 0;
  // Compute wavefronts for increasing distance
  for (distance=0;distance<max_distance;++distance) {
      //std::cerr << "distance is " << distance << std::endl;
    // Extend diagonally each wavefront point
    edit_wavefronts_extend_wavefront(
        wavefronts,extend_match,pattern_length,text_length,distance);
    // Exit condition
    if (target_k_abs <= distance &&
        wavefronts->wavefronts[distance].offsets[target_k] == target_offset) break;
    // Compute next wavefront starting point
    edit_wavefronts_compute_wavefront(
        wavefronts,pattern_length,text_length,distance+1);
  }
  // Backtrace wavefronts
  edit_wavefronts_backtrace(wavefronts,target_k,distance);
}
/*
int main(int argc,char* argv[]) {
  // Buffers
  char* pattern_mem =
      "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
      "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT"
      "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
  char* text_mem    =
      "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
      "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT"
      "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY";
  // Pattern & Text
  char* pattern = pattern_mem + 64;
  char* text = text_mem + 64;
  const int pattern_length = strlen(pattern_mem)-2*64;
  const int text_length = strlen(text_mem)-2*64;
  const int reps = 10000000;
  // Init Wavefronts
  edit_wavefronts_t wavefronts;
  edit_wavefronts_init(&wavefronts,pattern_length,text_length);
  int i;
  for (i=0;i<reps;++i) {
    edit_wavefronts_clean(&wavefronts);
    edit_wavefronts_align(&wavefronts,pattern,pattern_length,text,text_length);
  }
}
*/

}
