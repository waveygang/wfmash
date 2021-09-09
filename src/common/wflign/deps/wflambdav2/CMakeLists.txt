cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wflambda VERSION 1.0.0)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcx16 -march=native -g -Wno-pointer-arith -D__STDC_FORMAT_MACROS")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcx16 -march=native -g -Wno-pointer-arith -D__STDC_FORMAT_MACROS")

# Build wflambda as static library by default
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all libraries as shared")

# Define output directories for created binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(SOURCES
  wflambda/alignment/cigar.c
  wflambda/alignment/score_matrix.c
  wflambda/edit/edit_dp.c
  wflambda/gap_affine2p/affine2p_dp.c
  wflambda/gap_affine2p/affine2p_matrix.c
  wflambda/gap_affine2p/affine2p_penalties.c
  wflambda/gap_affine/affine_matrix.c
  wflambda/gap_affine/affine_penalties.c
  wflambda/gap_affine/swg.c
  wflambda/gap_lineal/nw.c
  wflambda/system/mm_allocator.c
  wflambda/system/mm_stack.c
  wflambda/system/profiler_counter.c
  wflambda/system/profiler_timer.c
  wflambda/utils/commons.c
  wflambda/utils/string_padded.c
  wflambda/utils/vector.c
  wflambda/wavefront/wavefront.c
  wflambda/wavefront/wavefront_align.c
  wflambda/wavefront/wavefront_aligner.c
  wflambda/wavefront/wavefront_backtrace.c
  wflambda/wavefront/wavefront_backtrace_buffer.c
  wflambda/wavefront/wavefront_components.c
  wflambda/wavefront/wavefront_compute.c
  wflambda/wavefront/wavefront_compute_affine.c
  wflambda/wavefront/wavefront_compute_affine2p.c
  wflambda/wavefront/wavefront_display.c
  wflambda/wavefront/wavefront_extend.c
  wflambda/wavefront/wavefront_pcigar.c
  wflambda/wavefront/wavefront_penalties.c
  wflambda/wavefront/wavefront_reduction.c
  wflambda/wavefront/wavefront_slab.c)

set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CXX)

# Create libraries.
add_library(wflambda ${SOURCES})
set_target_properties(wflambda
    PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION})

target_include_directories(wflambda PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

