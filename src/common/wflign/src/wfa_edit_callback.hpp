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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <functional>
#include <iostream>

namespace wflign {

/*
 * Translate k and offset to coordinates h,v
 */
#define EWAVEFRONT_V(k,offset) ((offset)-(k))
#define EWAVEFRONT_H(k,offset) (offset)

#define EWAVEFRONT_DIAGONAL(h,v) ((h)-(v))
#define EWAVEFRONT_OFFSET(h,v)   (h)

#define MAX(a,b) (((a)>=(b))?(a):(b))
#define ABS(a) (((a)>=0)?(a):-(a))

/*
 * Wavefront
 */
typedef int ewf_offset_t;  // Edit Wavefront Offset
typedef struct {
  int lo;                      // Effective lowest diagonal (inclusive)
  int hi;                      // Effective highest diagonal (inclusive)
  ewf_offset_t* offsets;       // Offsets
  ewf_offset_t* offsets_mem;   // Offsets memory
} edit_wavefront_t;
/*
 * Edit Wavefronts
 */
typedef struct {
  // Dimensions
  int pattern_length;
  int text_length;
  int max_distance;
  // Waves Offsets
  edit_wavefront_t* wavefronts;
  int wavefronts_allocated;
  // CIGAR
  char* edit_cigar;
  int edit_cigar_length;
} edit_wavefronts_t;

void edit_wavefronts_init(
    edit_wavefronts_t* const wavefronts,
    const int pattern_length,
    const int text_length);

void edit_wavefronts_clean(
    edit_wavefronts_t* const wavefronts);

edit_wavefront_t* edit_wavefronts_allocate_wavefront(
    edit_wavefronts_t* const edit_wavefronts,
    const int distance,
    const int lo_base,
    const int hi_base);

int edit_wavefronts_backtrace(
    edit_wavefronts_t* const wavefronts,
    const int target_k,
    const int target_distance);

void edit_wavefronts_extend_wavefront(
    edit_wavefronts_t* const wavefronts,
    const std::function<bool(const int&,const int&)>& extend_match,
    const int pattern_length,
    const int text_length,
    const int distance);

void edit_wavefronts_compute_wavefront(
    edit_wavefronts_t* const wavefronts,
    const int pattern_length,
    const int text_length,
    const int distance);

void edit_wavefronts_align(
    edit_wavefronts_t* const wavefronts,
    const std::function<bool(const int&,const int&)>& extend_match,
    const int pattern_length,
    const int text_length);

}
