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

#pragma once

#include "commons.hpp"

namespace wflambda {

/*
 * Checkers
 */
//#define WFLAMBDA_VECTOR_DEBUG

/*
 * Data Structures
 */
typedef struct {
  void* memory;
  uint64_t used;
  uint64_t element_size;
  uint64_t elements_allocated;
} wflambda_vector_t;

/*
 * Vector Setup (Initialization & Allocation)
 */
#define wflambda_vector_new(num_initial_elements,type) wflambda_vector_new_(num_initial_elements,sizeof(type))
wflambda_vector_t* wflambda_vector_new_(const uint64_t num_initial_elements,const uint64_t element_size);
void wflambda_vector_reserve(wflambda_vector_t* const vector,const uint64_t num_elements,const bool zero_mem);
void wflambda_vector_resize__clear(wflambda_vector_t* const vector,const uint64_t num_elements);
#define wflambda_vector_cast__clear(vector,type) wflambda_vector_cast__clear_s(vector,sizeof(type))
void wflambda_vector_cast__clear_(wflambda_vector_t* const vector,const uint64_t element_size);
#define wflambda_vector_clear(vector) (vector)->used=0
void wflambda_vector_delete(wflambda_vector_t* const vector);
#define wflambda_vector_is_empty(vector) (wflambda_vector_get_used(vector)==0)
#define wflambda_vector_reserve_additional(vector,additional) wflambda_vector_reserve(vector,wflambda_vector_get_used(vector)+additional,false)
#define wflambda_vector_prepare(vector,num_elements,type) \
  wflambda_vector_cast__clear(vector,sizeof(type)); \
  wflambda_vector_reserve(vector,num_elements,false);

/*
 * Element Getters/Setters
 */
#define wflambda_vector_get_mem(vector,type) ((type*)((vector)->memory))
#define wflambda_vector_get_last_elm(vector,type) (wflambda_vector_get_mem(vector,type)+(vector)->used-1)
#define wflambda_vector_get_free_elm(vector,type) (wflambda_vector_get_mem(vector,type)+(vector)->used)
#define wflambda_vector_set_elm(vector,position,type,elm) *wflambda_vector_get_elm(vector,position,type) = elm
#ifndef WFLAMBDA_VECTOR_DEBUG
  #define wflambda_vector_get_elm(vector,position,type) (wflambda_vector_get_mem(vector,type)+position)
#else
  void* wflambda_vector_get_mem_element(wflambda_vector_t* const vector,const uint64_t position,const uint64_t element_size);
  #define wflambda_vector_get_elm(vector,position,type) ((type*)wflambda_vector_get_mem_element(vector,position,sizeof(type)))
#endif

/*
 * Used elements Getters/Setters
 */
#define wflambda_vector_get_used(vector) ((vector)->used)
#define wflambda_vector_set_used(vector,total_used) (vector)->used=(total_used)
#define wflambda_vector_inc_used(vector) (++((vector)->used))
#define wflambda_vector_dec_used(vector) (--((vector)->used))
#define wflambda_vector_add_used(vector,additional) wflambda_vector_set_used(vector,wflambda_vector_get_used(vector)+additional)
#define wflambda_vector_update_used(vector,pointer_to_next_free_element) \
  (vector)->used = (pointer_to_next_free_element) - ((__typeof__(pointer_to_next_free_element))((vector)->memory))


/*
 * Vector Allocate/Insert (Get a new element or Add an element to the end of the vector)
 */
#define wflambda_vector_alloc_new(vector,type,return_element_pointer) { \
  wflambda_vector_reserve_additional(vector,1); \
  return_element_pointer = wflambda_vector_get_free_elm(vector,type); \
  wflambda_vector_inc_used(vector); \
}
#define wflambda_vector_insert(vector,element,type) { \
  wflambda_vector_reserve_additional(vector,1); \
  *(wflambda_vector_get_free_elm(vector,type))=element; \
  wflambda_vector_inc_used(vector); \
}

/*
 * Macro generic iterator
 *   WFLAMBDA_VECTOR_ITERATE(wflambda_vector_of_ints,elm_iterator,elm_counter,int) {
 *     ..code..
 *   }
 */
#define WFLAMBDA_VECTOR_ITERATE(vector,element,counter,type) \
  const uint64_t wflambda_vector_##element##_used = wflambda_vector_get_used(vector); \
  type* element = wflambda_vector_get_mem(vector,type); \
  uint64_t counter; \
  for (counter=0;counter<wflambda_vector_##element##_used;++element,++counter)
#define WFLAMBDA_VECTOR_ITERATE_OFFSET(vector,element,counter,offset,type) \
  const uint64_t wflambda_vector_##element##_used = wflambda_vector_get_used(vector); \
  type* element = wflambda_vector_get_mem(vector,type)+offset; \
  uint64_t counter; \
  for (counter=offset;counter<wflambda_vector_##element##_used;++counter,++element)
#define WFLAMBDA_VECTOR_ITERATE_CONST(vector,element,counter,type) \
  const uint64_t wflambda_vector_##element##_used = wflambda_vector_get_used(vector); \
  const type* element = wflambda_vector_get_mem(vector,type); \
  uint64_t counter; \
  for (counter=0;counter<wflambda_vector_##element##_used;++element,++counter)
#define WFLAMBDA_VECTOR_ITERATE_ELASTIC(vector,element,counter,type) \
  type* element = wflambda_vector_get_mem(vector,type); \
  uint64_t counter; \
  for (counter=0;counter<wflambda_vector_get_used(vector);++element,++counter)

/*
 * Miscellaneous
 */
void wflambda_vector_copy(wflambda_vector_t* const wflambda_vector_to,wflambda_vector_t* const wflambda_vector_from);
wflambda_vector_t* wflambda_vector_dup(wflambda_vector_t* const wflambda_vector_src);

}
