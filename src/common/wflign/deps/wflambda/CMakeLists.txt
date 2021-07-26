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
  wflambda/utils/commons.cpp
  wflambda/utils/vector.cpp
  wflambda/utils/string_padded.cpp
  wflambda/gap_affine/affine_wavefront_backtrace.cpp
  wflambda/gap_affine/affine_wavefront_reduction.cpp
  wflambda/gap_affine/affine_wavefront_penalties.cpp
  wflambda/gap_affine/affine_penalties.cpp
  wflambda/gap_affine/affine_wavefront_utils.cpp
  wflambda/gap_affine/affine_table.cpp
  wflambda/gap_affine/wavefront_stats.cpp
  wflambda/gap_affine/swg.cpp
  wflambda/gap_affine/affine_wavefront_display.cpp
  wflambda/gap_affine/affine_wavefront_align.cpp
  wflambda/gap_affine/affine_wavefront.cpp
  wflambda/gap_affine/affine_wavefront_extend.cpp
  wflambda/gap_lineal/nw.cpp
  wflambda/benchmark/benchmark_edit.cpp
  wflambda/benchmark/benchmark_utils.cpp
  wflambda/benchmark/benchmark_gap_lineal.cpp
  wflambda/benchmark/benchmark_gap_affine.cpp
  wflambda/system/profiler_timer.cpp
  wflambda/system/mm_allocator.cpp
  wflambda/system/profiler_counter.cpp
  wflambda/edit/edit_dp.cpp
  wflambda/edit/edit_cigar.cpp
  wflambda/edit/edit_table.cpp)

# Create libraries.
add_library(wflambda ${SOURCES})
set_target_properties(wflambda
    PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION})

target_include_directories(wflambda PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

