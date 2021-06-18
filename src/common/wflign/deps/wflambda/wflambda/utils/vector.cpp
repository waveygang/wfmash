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
 * DESCRIPTION: Simple linear vector for generic type elements
 */

#include "vector.hpp"

namespace wflambda {

/*
 * Constants
 */
#define WFLAMBDA_VECTOR_EXPAND_FACTOR (3.0/2.0)

/*
 * Setup
 */
wflambda_vector_t* wflambda_vector_new_(const uint64_t num_initial_elements,const uint64_t element_size) {
  wflambda_vector_t* const wflambda_vector_buffer = (wflambda_vector_t*)malloc(sizeof(wflambda_vector_t));
  wflambda_vector_buffer->element_size = element_size;
  wflambda_vector_buffer->elements_allocated = num_initial_elements;
  wflambda_vector_buffer->memory = malloc(num_initial_elements*element_size);
  if (!wflambda_vector_buffer->memory) {
      /*fprintf(stderr, "Could not create new vector (%"PRIu64" bytes requested)",
        num_initial_elements*element_size);*/
    exit(1);
  }
  wflambda_vector_buffer->used = 0;
  return wflambda_vector_buffer;
}
void wflambda_vector_reserve(wflambda_vector_t* const vector,const uint64_t num_elements,const bool zero_mem) {
  if (vector->elements_allocated < num_elements) {
    const uint64_t proposed=(float)vector->elements_allocated*WFLAMBDA_VECTOR_EXPAND_FACTOR;
    vector->elements_allocated = num_elements>proposed?num_elements:proposed;
    vector->memory = realloc(vector->memory,vector->elements_allocated*vector->element_size);
    if (!vector->memory) {
        /*fprintf(stderr, "Could not reserve vector (%"PRIu64" bytes requested)",
          vector->elements_allocated*vector->element_size);*/
      exit(1);
    }
  }
  if (zero_mem) {
      memset((uint8_t*)vector->memory+vector->used*vector->element_size,0,
        (vector->elements_allocated-vector->used)*vector->element_size);
  }
}
void wflambda_vector_resize__clear(wflambda_vector_t* const vector,const uint64_t num_elements) {
  if (vector->elements_allocated < num_elements) {
    const uint64_t proposed = (float)vector->elements_allocated*WFLAMBDA_VECTOR_EXPAND_FACTOR;
    vector->elements_allocated = (num_elements>proposed)?num_elements:proposed;
    // Free previous chunk (no need to pay the cost of reallocating memory)
    free(vector->memory);
    // Allocate new block of memory
    vector->memory = malloc(vector->elements_allocated*vector->element_size);
    if (!vector->memory) {
        /*fprintf(stderr, "Could not reserve vector (%"PRIu64" bytes requested)",
          vector->elements_allocated*vector->element_size);*/
      exit(1);
    }
  }
  vector->used=0;
}
void wflambda_vector_cast__clear_(wflambda_vector_t* const vector,const uint64_t element_size) {
  vector->elements_allocated = (vector->elements_allocated*vector->element_size)/element_size;
  vector->element_size = element_size;
  vector->used = 0;
}
void wflambda_vector_delete(wflambda_vector_t* const vector) {
  free(vector->memory);
  free(vector);
}
/*
 * Accessors
 */
#ifdef WFLAMBDA_VECTOR_DEBUG
void* wflambda_vector_get_mem_element(wflambda_vector_t* const vector,const uint64_t position,const uint64_t element_size) {
  if (position >= (vector)->used) {
      /*fprintf(stderr,"Vector position out-of-range [0,%"PRIu64")",(vector)->used);*/
    exit(1);
  }
  return vector->memory + (position*element_size);
}
#endif
/*
 * Miscellaneous
 */
void wflambda_vector_copy(wflambda_vector_t* const wflambda_vector_to,wflambda_vector_t* const wflambda_vector_from) {
  // Prepare
  wflambda_vector_cast__clear_(wflambda_vector_to,wflambda_vector_from->element_size);
  wflambda_vector_reserve(wflambda_vector_to,wflambda_vector_from->used,false);
  // Copy
  wflambda_vector_set_used(wflambda_vector_to,wflambda_vector_from->used);
  memcpy(wflambda_vector_to->memory,wflambda_vector_from->memory,wflambda_vector_from->used*wflambda_vector_from->element_size);
}
wflambda_vector_t* wflambda_vector_dup(wflambda_vector_t* const wflambda_vector_src) {
  wflambda_vector_t* const wflambda_vector_cpy = wflambda_vector_new_(wflambda_vector_src->used,wflambda_vector_src->element_size);
  // Copy
  wflambda_vector_set_used(wflambda_vector_cpy,wflambda_vector_src->used);
  memcpy(wflambda_vector_cpy->memory,wflambda_vector_src->memory,wflambda_vector_src->used*wflambda_vector_src->element_size);
  return wflambda_vector_cpy;
}

}
