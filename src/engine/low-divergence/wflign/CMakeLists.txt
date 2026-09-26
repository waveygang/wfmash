# lde_wflign: the wfmash 0.14-snapshot (7bf8988) wflign aligner, vendored
# for the low-divergence engine.  All first-party symbols carry the lde_
# prefix (see the file headers), and WFA2-lib is vendored separately as
# ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib with every linked C symbol renamed to lde_* and the
# C++ binding namespace renamed wfa -> lde_wfa, so the engine never
# touches the mainline WFA2 build.  Robin-hood hashing and PicoSHA2
# remain shared with the mainline build (header-only, identical trees).
add_library(lde_wfa2 STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_align.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_aligner.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_attributes.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_backtrace.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_backtrace_buffer.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_backtrace_offload.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_bialign.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_bialigner.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_components.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_compute.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_compute_affine.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_compute_affine2p.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_compute_edit.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_compute_linear.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_debug.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_display.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_extend.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_extend_kernels.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_extend_kernels_avx.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_heuristic.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_pcigar.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_penalties.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_plot.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_sequences.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_slab.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_termination.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront/wavefront_unialign.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/alignment/affine2p_penalties.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/alignment/affine_penalties.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/alignment/cigar.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/alignment/cigar_utils.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/alignment/score_matrix.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/system/mm_allocator.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/system/mm_stack.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/system/profiler_counter.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/system/profiler_timer.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/utils/bitmap.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/utils/commons.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/utils/dna_text.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/utils/heatmap.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/utils/sequence_buffer.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/utils/vector.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/bindings/cpp/WFAligner.cpp)
target_include_directories(lde_wfa2 PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/wavefront
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/alignment
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/system
  ${CMAKE_CURRENT_SOURCE_DIR}/../deps/lde_WFA2-lib/utils)
target_compile_options(lde_wfa2 PRIVATE -Wno-everything)

add_library(lde_wflign STATIC
  src/lde_wflign.cpp
  src/lde_wflign_alignment.cpp
  src/lde_wflign_patch.cpp
  src/lde_rkmh.cpp
  src/lde_murmur3.cpp
  deps/lde_atomic_image.cpp
  deps/lde_lodepng.cpp)
target_include_directories(lde_wflign PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/deps)
target_include_directories(lde_wflign PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/engine/low-divergence/deps)
target_include_directories(lde_wflign PRIVATE
  ${CMAKE_SOURCE_DIR}/src/common/wflign/deps
  ${CMAKE_SOURCE_DIR}/src/common/wflign/deps/robin-hood-hashing)
target_compile_options(lde_wflign PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(lde_wflign PUBLIC lde_wfa2 z)

# Standalone engine test binary: only the low-divergence engine, no
# mainline sources linked at all. Used to isolate mainline-vs-engine
# link-context effects.  It assumes system htslib/gsl/z libraries and
# default include paths, so it is only configured for non-vendored,
# non-static builds; under VENDOR_EVERYTHING/BUILD_STATIC the vendored
# headers and archives live outside the default paths and the target
# would not resolve (this broke the macOS CI build).
if (NOT BUILD_STATIC AND NOT VENDOR_EVERYTHING)
add_executable(lde_standalone ${CMAKE_SOURCE_DIR}/src/engine/low-divergence/lde_main.cpp)
target_include_directories(lde_standalone PRIVATE
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/common
  ${CMAKE_SOURCE_DIR}/src/common/wflign/deps
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/src/engine/low-divergence/deps)
target_link_libraries(lde_standalone lde_wflign hts z pthread gsl gslcblas)

# tiny main() for the standalone engine test binary
file(WRITE ${CMAKE_BINARY_DIR}/lde_standalone_main.cpp
"int lde_main(int argc, char** argv);\nint main(int argc, char** argv) { return lde_main(argc, argv); }\n")
target_sources(lde_standalone PRIVATE ${CMAKE_BINARY_DIR}/lde_standalone_main.cpp)
endif()
