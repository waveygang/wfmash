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
 * DESCRIPTION: Simple managed-memory allocator that reduces the overhead
 *   of using malloc/calloc/free functions by allocating slabs of memory
 *   and dispatching memory segments in order.
 */

#include "mm_allocator.hpp"

namespace wflambda {

/*
 * Debug
 */
//#define MM_ALLOCATOR_MALLOC

/*
 * Constants
 */
#define MM_ALLOCATOR_SEGMENT_INITIAL_REQUESTS   10000
#define MM_ALLOCATOR_INITIAL_SEGMENTS              10
#define MM_ALLOCATOR_INITIAL_MALLOC_REQUESTS       10
#define MM_ALLOCATOR_INITIAL_STATES                10

/*
 * Allocator Segments Freed Cond
 */
#define MM_ALLOCATOR_FREED_FLAG                 0x80000000ul
#define MM_ALLOCATOR_REQUEST_IS_FREE(request)  ((request)->size & MM_ALLOCATOR_FREED_FLAG)
#define MM_ALLOCATOR_REQUEST_SET_FREE(request) ((request)->size |= MM_ALLOCATOR_FREED_FLAG)
#define MM_ALLOCATOR_REQUEST_SIZE(request)     ((request)->size & ~(MM_ALLOCATOR_FREED_FLAG))

/*
 * Memory Request
 */
typedef struct {
  // Request
  uint32_t offset;
  uint32_t size;
  // Log
#ifdef MM_ALLOCATOR_LOG
  uint64_t timestamp;
  char* func_name;
  uint64_t line_no;
#endif
} mm_allocator_request_t;
/*
 * Memory Segments
 */
typedef struct {
  // Index (ID)
  uint64_t segment_idx;         // Index in the segments vector
  // Memory
  uint64_t segment_size;        // Total memory available
  void* memory;                 // Memory
  uint64_t used;                // Bytes used (offset to memory next free byte)
  // Requests
  vector_t* requests;           // Memory requests (mm_allocator_request_t)
} mm_allocator_segment_t;
/*
 * Reference (Header of every memory allocated)
 */
typedef struct {
  uint32_t segment_idx;
  uint32_t request_idx;
} mm_allocator_reference_t;

/*
 * Segments
 */
mm_allocator_segment_t* mm_allocator_segment_new(
    mm_allocator_t* const mm_allocator) {
  // Allocate new segment
  mm_allocator_segment_t* const segment = (mm_allocator_segment_t*)malloc(sizeof(mm_allocator_segment_t));
  // Index
  const uint64_t segment_idx = vector_get_used(mm_allocator->segments);
  segment->segment_idx = segment_idx;
  // Memory
  segment->segment_size = mm_allocator->segment_size;
  segment->memory = malloc(mm_allocator->segment_size);
  segment->used = 0;
  // Requests
  segment->requests = vector_new(MM_ALLOCATOR_SEGMENT_INITIAL_REQUESTS,mm_allocator_request_t);
  // Add to segments
  vector_insert(mm_allocator->segments,segment,mm_allocator_segment_t*);
  // Return
  return segment;
}
void mm_allocator_segment_clear(
    mm_allocator_segment_t* const segment) {
  segment->used = 0;
  vector_clear(segment->requests);
}
void mm_allocator_segment_delete(
    mm_allocator_segment_t* const segment) {
  vector_delete(segment->requests);
  free(segment->memory);
  free(segment);
}
mm_allocator_request_t* mm_allocator_segment_get_request(
    mm_allocator_segment_t* const segment,
    const uint64_t request_idx) {
  return vector_get_elm(segment->requests,request_idx,mm_allocator_request_t);
}
uint64_t mm_allocator_segment_get_num_requests(
    mm_allocator_segment_t* const segment) {
  return vector_get_used(segment->requests);
}
/*
 * Setup
 */
mm_allocator_t* mm_allocator_new(
    const uint64_t segment_size) {
  // Allocate handler
  mm_allocator_t* const mm_allocator = (mm_allocator_t*)malloc(sizeof(mm_allocator_t));
  mm_allocator->request_ticker = 0;
  // Segments
  mm_allocator->segment_size = segment_size;
  mm_allocator->current_segment_idx = 0;
  mm_allocator->segments = vector_new(MM_ALLOCATOR_INITIAL_SEGMENTS,mm_allocator_segment_t*);
  mm_allocator->segments_free = vector_new(MM_ALLOCATOR_INITIAL_SEGMENTS,mm_allocator_segment_t*);
  // Allocate an initial segment
  mm_allocator_segment_new(mm_allocator);
  // Malloc Memory
  mm_allocator->malloc_requests = vector_new(MM_ALLOCATOR_INITIAL_MALLOC_REQUESTS,void*);
  // Return
  return mm_allocator;
}
void mm_allocator_clear(
    mm_allocator_t* const mm_allocator) {
  // Clear segments
  vector_clear(mm_allocator->segments_free);
  const uint64_t num_segments = vector_get_used(mm_allocator->segments);
  mm_allocator_segment_t** const segments = 
      vector_get_mem(mm_allocator->segments,mm_allocator_segment_t*);
  uint64_t i;
  for (i=1;i<num_segments;++i) {
    mm_allocator_segment_clear(segments[i]); // Clear segment
    vector_insert(mm_allocator->segments_free,segments[i],mm_allocator_segment_t*); // Add to free segments
  }
  mm_allocator->current_segment_idx = 0;
  // Clear malloc memory
  VECTOR_ITERATE(mm_allocator->malloc_requests,malloc_request,m,void*) {
    free(*malloc_request); // Free malloc requests
  }
  vector_clear(mm_allocator->malloc_requests);
}
void mm_allocator_delete(
    mm_allocator_t* const mm_allocator) {
  // Free segments
  VECTOR_ITERATE(mm_allocator->segments,segment_ptr,p,mm_allocator_segment_t*) {
    mm_allocator_segment_delete(*segment_ptr);
  }
  vector_delete(mm_allocator->segments);
  vector_delete(mm_allocator->segments_free);
  // Free malloc memory
  VECTOR_ITERATE(mm_allocator->malloc_requests,malloc_request,m,void*) {
    free(*malloc_request); // Free remaining requests
  }
  vector_delete(mm_allocator->malloc_requests);
  // Free handler
  free(mm_allocator);
}
/*
 * Accessors
 */
mm_allocator_segment_t* mm_allocator_get_segment(
    mm_allocator_t* const mm_allocator,
    const uint64_t segment_idx) {
  return *(vector_get_elm(mm_allocator->segments,segment_idx,mm_allocator_segment_t*));
}
uint64_t mm_allocator_get_num_segments(
    mm_allocator_t* const mm_allocator) {
  return vector_get_used(mm_allocator->segments);
}
/*
 * Allocate
 */
mm_allocator_segment_t* mm_allocator_fetch_segment(
    mm_allocator_t* const mm_allocator,
    const uint64_t num_bytes) {
  // Fetch current segment
  mm_allocator_segment_t* const curr_segment =
      mm_allocator_get_segment(mm_allocator,mm_allocator->current_segment_idx);
  // Check available segment size
  if (curr_segment->used + num_bytes <= curr_segment->segment_size) {
    return curr_segment;
  }
  // Check overall segment size
  if (num_bytes > curr_segment->segment_size) {
    return NULL; // Memory request over segment size
  }
  // Get free segment
  const uint64_t free_segments = vector_get_used(mm_allocator->segments_free);
  if (free_segments > 0) {
    mm_allocator_segment_t* const segment =
        *vector_get_elm(mm_allocator->segments_free,free_segments-1,mm_allocator_segment_t*);
    vector_dec_used(mm_allocator->segments_free);
    mm_allocator->current_segment_idx = segment->segment_idx;
    return segment;
  }
  // Allocate new segment
  mm_allocator_segment_t* const segment = mm_allocator_segment_new(mm_allocator);
  mm_allocator->current_segment_idx = segment->segment_idx;
  return segment;
}
void* mm_allocator_allocate(
    mm_allocator_t* const mm_allocator,
    uint64_t num_bytes,
    const bool zero_mem
#ifdef MM_ALLOCATOR_LOG
    ,const char* func_name,
    uint64_t line_no
#endif
    ) {
  // Zero check
  if (num_bytes == 0) {
    fprintf(stderr,"MM-Allocator error. Zero bytes request\n");
    exit(1);
  }
  // Add mm_reference payload
  num_bytes += sizeof(mm_allocator_reference_t);
  if (num_bytes%16 != 0) { // Ensure alignment to 16 Bytes boundaries (128bits)
    num_bytes += 16 - (num_bytes%16); // Serve more bytes
  }
  // Fetch segment
#ifdef MM_ALLOCATOR_MALLOC
  mm_allocator_segment_t* const segment = NULL; // Force malloc memory
#else
  mm_allocator_segment_t* const segment = mm_allocator_fetch_segment(mm_allocator,num_bytes);
#endif
  if (segment != NULL) {
    // Allocate memory
    void* const memory = segment->memory + segment->used;
    if (zero_mem) memset(memory,0,num_bytes); // Set zero
    // Set mm_reference
    mm_allocator_reference_t* const mm_reference = (mm_allocator_reference_t* const)memory;
    mm_reference->segment_idx = segment->segment_idx;
    mm_reference->request_idx = mm_allocator_segment_get_num_requests(segment);
    // Add request
    mm_allocator_request_t* request;
    vector_alloc_new(segment->requests,mm_allocator_request_t,request);
    request->offset = segment->used;
    request->size = num_bytes;
#ifdef MM_ALLOCATOR_LOG
    request->timestamp = (mm_allocator->request_ticker)++;
    request->func_name = (char*)func_name;
    request->line_no = line_no;
#endif
    // Update segment
    segment->used += num_bytes;
    // Return memory
    return memory + sizeof(mm_allocator_reference_t);
  } else {
    // Allocate memory
    void* const memory = malloc(num_bytes);
    if (zero_mem) memset(memory,0,num_bytes); // Set zero
    vector_insert(mm_allocator->malloc_requests,memory,void*);
    // Set reference
    mm_allocator_reference_t* const mm_reference = (mm_allocator_reference_t* const)memory;
    mm_reference->segment_idx = UINT32_MAX;
    // Return memory
    return memory + sizeof(mm_allocator_reference_t);
  }
}
/*
 * Allocator Free
 */
void mm_allocator_free_malloc_request(
    mm_allocator_t* const mm_allocator,
    void* memory) {
  // Search for address in malloc requests
  const uint64_t num_malloc_requests = vector_get_used(mm_allocator->malloc_requests);
  void** const malloc_requests = vector_get_mem(mm_allocator->malloc_requests,void*);
  uint64_t i;
  for (i=0;i<num_malloc_requests;++i) { // TODO Needed just for report
    if (malloc_requests[i] == memory) {
      // Free memory
      free(memory);
      // Shift requests
      for (;i<num_malloc_requests-1;++i) {
        malloc_requests[i] = malloc_requests[i+1];
      }
      vector_dec_used(mm_allocator->malloc_requests);
      // Ok
      return;
    }
  }
  // Address not found. Raise error
  fprintf(stderr,"MM-Allocator error. Invalid address freed (request not found)\n");
  exit(1);
}
void mm_allocator_free_allocator_request(
    mm_allocator_t* const mm_allocator,
    mm_allocator_segment_t* const segment,
    const uint32_t request_idx,
    mm_allocator_request_t* const request) {
  if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
    fprintf(stderr,"MM-Allocator error: double free\n");
    exit(1);
  }
  // Free request
  MM_ALLOCATOR_REQUEST_SET_FREE(request);
  // Free contiguous request(s) at the end of the segment
  uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
  if (request_idx == num_requests-1) { // Is the last request?
    --num_requests;
    mm_allocator_request_t* request =
        vector_get_mem(segment->requests,mm_allocator_request_t) + (num_requests-1);
    while (num_requests>0 && MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
      --num_requests; // Free request
      --request;
    }
    // Update segment used
    if (num_requests > 0) {
      segment->used = request->offset + request->size;
      vector_set_used(segment->requests,num_requests);
    } else {
      // Segment fully freed
      mm_allocator_segment_clear(segment); // Clear
      // Add to free segments (if it is not the current segment)
      if (segment->segment_idx != mm_allocator->current_segment_idx) {
        vector_insert(mm_allocator->segments_free,segment,mm_allocator_segment_t*);
      }
    }
  }
}
void mm_allocator_free(
    mm_allocator_t* const mm_allocator,
    void* const memory) {
  // Get mm_reference
  void* const effective_memory = memory - sizeof(mm_allocator_reference_t);
  mm_allocator_reference_t* const mm_reference = (mm_allocator_reference_t* const)effective_memory;
  if (mm_reference->segment_idx == UINT32_MAX) {
    // Malloc memory
    mm_allocator_free_malloc_request(mm_allocator,effective_memory);
  } else {
    // Allocator Memory
    mm_allocator_segment_t* const segment =
        mm_allocator_get_segment(mm_allocator,mm_reference->segment_idx);
    mm_allocator_request_t* const request =
        mm_allocator_segment_get_request(segment,mm_reference->request_idx);
    mm_allocator_free_allocator_request(
        mm_allocator,segment,mm_reference->request_idx,request);
  }
}
/*
 * Utils
 */
void mm_allocator_get_occupation(
    mm_allocator_t* const mm_allocator,
    uint64_t* const bytes_used,
    uint64_t* const bytes_free_available,
    uint64_t* const bytes_free_fragmented) {
  // Init
  *bytes_used = 0;
  *bytes_free_available = 0;
  *bytes_free_fragmented = 0;
  // Traverse all segments
  const int64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  int64_t segment_idx, request_idx;
  for (segment_idx=0;segment_idx<num_segments;++segment_idx) {
    mm_allocator_segment_t* const segment = mm_allocator_get_segment(mm_allocator,segment_idx);
    const uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
    bool free_memory = true;
    for (request_idx=num_requests-1;request_idx>=0;--request_idx) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,request_idx);
      const uint64_t size = MM_ALLOCATOR_REQUEST_SIZE(request);
      if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
        if (free_memory) {
          *bytes_free_available += size;
        } else {
          *bytes_free_fragmented += size;
        }
      } else {
        free_memory = false;
        *bytes_used += size;
      }
    }
    // Account for free space at the end of the segment
    if (num_requests > 0) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,num_requests-1);
      *bytes_free_available += segment->used - (request->offset+request->size);
    }
  }
}
/*
 * Display
 */
void mm_allocator_print_request(
    FILE* const stream,
    mm_allocator_request_t* const request,
    const uint64_t segment_idx,
    const uint64_t request_idx) {
      fprintf(stream,"    [#%03" PRIu64 "/%05" PRIu64 "\t%s\t@%08u\t(%" PRIu64 " Bytes)"
#ifdef MM_ALLOCATOR_LOG
          "\t%s:%" PRIu64 "\t{ts=%" PRIu64 "}"
#endif
          "\n",
          segment_idx,
          request_idx,
          MM_ALLOCATOR_REQUEST_IS_FREE(request) ? "Free]     " : "Allocated]",
          request->offset,
          (uint64_t)MM_ALLOCATOR_REQUEST_SIZE(request)
#ifdef MM_ALLOCATOR_LOG
          ,request->func_name,
          request->line_no,
          request->timestamp
#endif
      );
}
void mm_allocator_print_requests(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool compact_free) {
  uint64_t segment_idx, request_idx;
  uint64_t free_block = 0;
  fprintf(stream,"  => Memory.requests\n");
  // Traverse all segments
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  for (segment_idx=0;segment_idx<num_segments;++segment_idx) {
    mm_allocator_segment_t* const segment = mm_allocator_get_segment(mm_allocator,segment_idx);
    const uint64_t num_requests = mm_allocator_segment_get_num_requests(segment);
    for (request_idx=0;request_idx<num_requests;++request_idx) {
      mm_allocator_request_t* const request = mm_allocator_segment_get_request(segment,request_idx);
      if (compact_free) {
        if (MM_ALLOCATOR_REQUEST_IS_FREE(request)) {
          free_block += MM_ALLOCATOR_REQUEST_SIZE(request);
        } else {
          if (free_block > 0) {
            fprintf(stream,"    [n/a\tFree]      \t(%" PRIu64 " Bytes)\n",free_block);
            free_block = 0;
          }
          mm_allocator_print_request(stream,request,segment_idx,request_idx);
        }
      } else {
        mm_allocator_print_request(stream,request,segment_idx,request_idx);
      }
    }
    if (request_idx==0) {
      fprintf(stream,"    -- No requests --\n");
    }
  }
}
void mm_allocator_print(
    FILE* const stream,
    mm_allocator_t* const mm_allocator,
    const bool display_requests) {
  // Print header
  fprintf(stream,"MM-Allocator.report\n");
  // Print segment information
  const uint64_t num_segments = mm_allocator_get_num_segments(mm_allocator);
  const uint64_t segment_size = mm_allocator_get_segment(mm_allocator,0)->segment_size;
  fprintf(stream,"  => Segments.allocated %" PRIu64 "\n",num_segments);
  fprintf(stream,"    => Segments.size %" PRIu64 " MB\n",segment_size/(1024*1024));
  fprintf(stream,"    => Memory.available %" PRIu64 " MB\n",num_segments*(segment_size/(1024*1024)));
  // Print memory information
  uint64_t bytes_used, bytes_free_available, bytes_free_fragmented;
  mm_allocator_get_occupation(mm_allocator,&bytes_used,&bytes_free_available,&bytes_free_fragmented);
  fprintf(stream,"  => Memory.used %" PRIu64 "\n",bytes_used);
  fprintf(stream,"  => Memory.free %" PRIu64 "\n",bytes_free_available+bytes_free_fragmented);
  fprintf(stream,"    => Memory.free.available  %" PRIu64 "\n",bytes_free_available);
  fprintf(stream,"    => Memory.free.fragmented %" PRIu64 "\n",bytes_free_fragmented);
  // Print memory requests
  if (display_requests) {
    mm_allocator_print_requests(stream,mm_allocator,false);
  }
}


}
