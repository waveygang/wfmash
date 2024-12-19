# Usage:
#
#   mkdir build ; cd build
#   cmake -DCMAKE_BUILD_TYPE=Debug -DASAN=ON ..
#   make -j 8 VERBOSE=1
#
#   Build types are Release (default), Debug, Generic, RelWithDebInfo
#
#   Debug builds have options:
#
#      cmake -DCMAKE_BUILD_TYPE=Debug -DASAN=ON     - check for memory leaks
#      cmake -DCMAKE_BUILD_TYPE=Debug -DPROFILER=ON - optionally use profiler
#
#   See also the header of ./guix.scm

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(wfmash)

enable_testing()

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD 20) # because of sequenceIds.hpp:101:23: warning: captured structured bindings are a C++20 extension [-Wc++20-extensions]
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard is not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


option(BUILD_STATIC "Build static binary" OFF)
option(BUILD_DEPS "Build external dependencies (not recommended)" OFF)
option(BUILD_OPTIMIZED "Build optimized binary" OFF)
option(PROFILER "Enable profiling" OFF)
option(ASAN "Use address sanitiser (Debug build only)" OFF)
option(DISABLE_LTO "Disable IPO/LTO" OFF)
option(STOP_ON_ERROR "Stop compiling on first error" OFF)

if (NOT DISABLE_LTO)
  include(CheckIPOSupported) # adds lto
  check_ipo_supported(RESULT ipo_supported OUTPUT output) # lto
endif()

if (STOP_ON_ERROR)
  add_definitions( -Wfatal-errors )
endif()

if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_POSITION_INDEPENDENT_CODE OFF)
else ()
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

find_package(PkgConfig REQUIRED)
include(FindCURL) # for htslib?
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Release|Debug|Generic|RelWithDebInfo (default for distros)." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  if (BUILD_OPTIMIZED)
    # set(EXTRA_FLAGS "-Ofast -march=native -pipe -msse4.2 -funroll-all-loops") #  -fprofile-generate=../pgo")
    # set(EXTRA_FLAGS "-Ofast -march=x86-64-v3 -funroll-all-loops")
    set(EXTRA_FLAGS "-O3 -march=native -funroll-all-loops")
    # set(EXTRA_FLAGS "-O3 -march=native -funroll-all-loops -fprofile-generate=${CMAKE_BINARY_DIR}/../pgo")
    # set(EXTRA_FLAGS "-O3 -march=native -funroll-all-loops -fprofile-use=${CMAKE_BINARY_DIR}/../pgo")
  else()
    set(EXTRA_FLAGS "-Ofast -march=native")
  endif()
  set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG ${EXTRA_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG ${EXTRA_FLAGS}")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Generic)
  set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
  set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
  if (ASAN)
    # Enable ASan
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
  endif ()
else ()
    set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif ()

if( NOT DISABLE_LTO AND ipo_supported )
  message(STATUS "IPO / LTO enabled")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${ipo_supported})
endif()

if(PROFILER)
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -g") # add debug switch for resolution
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

if(GPROF) # for the GNU profiler
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
endif(GPROF)

if (WFA_PNG_AND_TSV)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWFA_PNG_TSV_TIMING")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWFA_PNG_TSV_TIMING")
endif ()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src/common/wflign EXCLUDE_FROM_ALL)

include(ExternalProject)

if (BUILD_DEPS)
  ExternalProject_Add(htslib
    URL https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/htslib
    CONFIGURE_COMMAND autoreconf -i && ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/htslib --disable-libcurl --disable-s3
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install
    BUILD_IN_SOURCE 1
  )

  ExternalProject_Add(gsl
    URL https://mirror.ibcp.fr/pub/gnu/gsl/gsl-2.8.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gsl
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/gsl
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install
    BUILD_IN_SOURCE 1
  )

  ExternalProject_Add(libdeflate
    URL https://github.com/ebiggers/libdeflate/releases/download/v1.20/libdeflate-1.20.tar.gz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/libdeflate
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/libdeflate
      -DCMAKE_BUILD_TYPE=Release
    BUILD_COMMAND cmake --build . --config Release
    INSTALL_COMMAND cmake --install . --config Release
    BUILD_IN_SOURCE 1
  )

  add_dependencies(wfmash htslib gsl libdeflate)
endif()

add_executable(wfmash
  src/common/utils.cpp
  src/interface/main.cpp)

if (BUILD_DEPS)
  target_include_directories(wfmash PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/htslib/include
    ${CMAKE_CURRENT_BINARY_DIR}/gsl/include
    ${CMAKE_CURRENT_BINARY_DIR}/libdeflate/include
  )

  target_link_libraries(wfmash
    ${CMAKE_CURRENT_BINARY_DIR}/gsl/lib/libgsl.a
    ${CMAKE_CURRENT_BINARY_DIR}/gsl/lib/libgslcblas.a
    ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
    ${CMAKE_CURRENT_BINARY_DIR}/libdeflate/lib/libdeflate.a
  )
else()
  #find_package(HTSLIB REQUIRED)
  #find_package(GSL REQUIRED)
  #find_package(LibDeflate REQUIRED)

  #target_include_directories(wfmash PRIVATE
  #  ${HTSLIB_INCLUDE_DIR}
  #  ${GSL_INCLUDE_DIR}
  #  ${LIBDEFLATE_INCLUDE_DIR}
  #)

  target_link_libraries(wfmash
    gsl
    gslcblas
    hts
    deflate
  )
endif()



target_include_directories(wfmash PRIVATE
  src
  src/common
  src/common/wflign/deps
  src/common/wflign/deps/WFA2-lib
)

target_link_libraries(wfmash
  m
  pthread
  rt
  lzma
  bz2
  z
  Threads::Threads
)

if (PROFILER)
  target_link_libraries(wfmash profiler)
endif()

if (BUILD_STATIC)
  target_link_libraries(wfmash
    wflign_static
    wfa2cpp_static
  )
else()
  target_link_libraries(wfmash
    wflign
    wfa2cpp
  )
endif()

# This is to disable tests defined in CTestCustom.cmake:
configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR})

if (ASAN)
  set(INVOKE "ASAN_OPTIONS=detect_leaks=1:symbolize=1:verify_asan_link_order=0 LSAN_OPTIONS=verbosity=0:log_threads=1:suppressions=${CMAKE_SOURCE_DIR}/test/libasan-suppress.txt ${CMAKE_BINARY_DIR}/bin/wfmash")
else ()
  set(INVOKE "${CMAKE_BINARY_DIR}/bin/wfmash")
endif()

if (PROFILER)
  set(INVOKE "CPUPROFILE=wfmash.prof ${INVOKE}")
endif()

set(REGRESSION_TEST_DIR "test/data/regression")

add_test(
  NAME wfmash-time-LPA
  COMMAND bash -c "${INVOKE} data/LPA.subset.fa.gz -p 80 -n 5 -t 8 > wfmash-time-LPA.output"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-subset-LPA-to-SAM
  COMMAND bash -c "${INVOKE} data/LPA.subset.fa.gz -N -a -L > LPA.subset.sam && samtools view LPA.subset.sam -bS | samtools sort > LPA.subset.bam && samtools index LPA.subset.bam && samtools view LPA.subset.bam | head | cut -f 1-9"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF
  COMMAND bash -c "${INVOKE} data/scerevisiae8.fa.gz -p 95 -n 7 -m -L -Y \\# > scerevisiae8.paf && ./scripts/test.sh data/scerevisiae8.fa.gz.fai scerevisiae8.paf 0.92 && head -3000 scerevisiae8.paf > scerevisiae8.paf.output && ${CMAKE_COMMAND} -E compare_files ${REGRESSION_TEST_DIR}/scerevisiae8.paf.output scerevisiae8.paf.output"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-short-reads-500bps-to-SAM
  COMMAND bash -c "${INVOKE} data/reference.fa.gz data/reads.500bps.fa.gz -s 0.5k -N -a > reads.500bps.sam && samtools view reads.500bps.sam -bS | samtools sort > reads.500bps.bam && samtools index reads.500bps.bam && samtools view reads.500bps.bam | head > wfmash-short-reads-500bps-to-SAM.output && ${CMAKE_COMMAND} -E compare_files ${REGRESSION_TEST_DIR}/wfmash-short-reads-500bps-to-SAM.output wfmash-short-reads-500bps-to-SAM.output"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-short-reads-255bps-to-PAF
  COMMAND bash -c "${INVOKE} data/reads.255bps.fa.gz -w 16 -s 100 -L > reads.255bps.paf && head reads.255bps.paf && ${CMAKE_COMMAND} -E compare_files ${REGRESSION_TEST_DIR}/reads.255bps.paf reads.255bps.paf"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-input-mapping
  COMMAND bash -c "${INVOKE} data/scerevisiae8.fa.gz -p 95 -T S288C -Q SK1 -m > mappings.paf && \
${INVOKE} data/scerevisiae8.fa.gz -i mappings.paf|sort > aligned.paf.output && \
${CMAKE_COMMAND} -E compare_files ${REGRESSION_TEST_DIR}/aligned.paf.output aligned.paf.output"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-all2all
  COMMAND bash -c "${INVOKE} -t 8 data/scerevisiae8.fa.gz > all2all.paf && \
  sort all2all.paf | head -300 > all2all-300.paf.output && \
  ${CMAKE_COMMAND} -E compare_files ${REGRESSION_TEST_DIR}/all2all-300.paf.output all2all-300.paf.output"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS wfmash DESTINATION bin)

if (BUILD_STATIC)

  install(TARGETS wfa2cpp_static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(TARGETS wflign_static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(TARGETS wfa2_static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
else()

  install(TARGETS wfa2cpp
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(TARGETS wflign
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(TARGETS wfa2
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/generate_git_version.sh ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/common/wflign/src)
