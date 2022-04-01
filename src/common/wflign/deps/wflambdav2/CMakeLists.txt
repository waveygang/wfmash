cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wflambda VERSION 1.0.0)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.

#set(CMAKE_BUILD_TYPE Release)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Release Debug Generic." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if (${CMAKE_BUILD_TYPE} MATCHES Release)
    set(EXTRA_FLAGS "-Ofast -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    # Debug use the defaults - so commenting out:
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -g -fsanitize=address -Wno-pointer-arith -D__STDC_FORMAT_MACROS")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -g -fsanitize=address -Wno-pointer-arith -D__STDC_FORMAT_MACROS")
else()
    # Use all standard-compliant optimizations - always add these:
    set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS} -Wno-pointer-arith -D__STDC_FORMAT_MACROS")
    set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS} -Wno-pointer-arith -D__STDC_FORMAT_MACROS")
endif ()

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

