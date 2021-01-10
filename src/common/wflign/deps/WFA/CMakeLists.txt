cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wfa VERSION 1.0.0)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcx16 -march=native -g -Wno-pointer-arith")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcx16 -march=native -g -Wno-pointer-arith")

# Build wfa as static library by default
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all libraries as shared")

# Define output directories for created binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(SOURCES
  WFA/utils/commons.cpp
  WFA/utils/vector.cpp
  WFA/utils/string_padded.cpp
  WFA/gap_affine/affine_wavefront_backtrace.cpp
  WFA/gap_affine/affine_wavefront_reduction.cpp
  WFA/gap_affine/affine_wavefront_penalties.cpp
  WFA/gap_affine/affine_penalties.cpp
  WFA/gap_affine/affine_wavefront_utils.cpp
  WFA/gap_affine/affine_table.cpp
  WFA/gap_affine/wavefront_stats.cpp
  WFA/gap_affine/swg.cpp
  WFA/gap_affine/affine_wavefront_display.cpp
  WFA/gap_affine/affine_wavefront_align.cpp
  WFA/gap_affine/affine_wavefront.cpp
  WFA/gap_affine/affine_wavefront_extend.cpp
  WFA/gap_lineal/nw.cpp
  #WFA/benchmark/benchmark_edit.cpp
  #WFA/benchmark/benchmark_utils.cpp
  #WFA/benchmark/benchmark_gap_lineal.cpp
  #WFA/benchmark/benchmark_gap_affine.cpp
  WFA/system/profiler_timer.cpp
  WFA/system/mm_allocator.cpp
  WFA/system/profiler_counter.cpp
  WFA/edit/edit_dp.cpp
  WFA/edit/edit_cigar.cpp
  WFA/edit/edit_table.cpp)

# Create libraries.
add_library(wfa ${SOURCES})
set_target_properties(wfa
    PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION})

target_include_directories(wfa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

