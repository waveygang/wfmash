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

include(GNUInstallDirs) # For CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_BINDIR etc.
include(CheckCXXCompilerFlag)

# --- Start: RPATH settings for macOS ---
if(APPLE)
  set(CMAKE_MACOSX_RPATH ON) # Enable RPATH generation on macOS.
                           # This uses @rpath for linking against dylibs in the build tree
                           # and sets install_name to use @rpath.

  # Set the RPATH for installed targets.
  # For an executable in ${CMAKE_INSTALL_PREFIX}/bin and libraries in ${CMAKE_INSTALL_PREFIX}/lib,
  # @executable_path/../${CMAKE_INSTALL_LIBDIR} allows the executable to find the libraries.
  # ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR} is an absolute path.
  set(CMAKE_INSTALL_RPATH "@executable_path/../${CMAKE_INSTALL_LIBDIR};${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  # It's also good practice to use the install RPATHs during the build phase.
  # This helps ensure that targets run correctly from the build directory
  # and can help catch RPATH issues earlier.
  set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
endif()
# --- End: RPATH settings for macOS ---

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(BUILD_STATIC "Build static binary" OFF)
option(BUILD_DEPS "Build external dependencies (not recommended)" OFF)
option(BUILD_OPTIMIZED "Build optimized binary" OFF)
option(PROFILER "Enable profiling" OFF)
option(ASAN "Use address sanitiser (Debug build only)" OFF)
option(DISABLE_LTO "Disable IPO/LTO" OFF)
option(STOP_ON_ERROR "Stop compiling on first error" OFF)
option(WFA_PNG_AND_TSV "Enable PNG and TSV output for WFA (for debugging/visualization)" OFF)

if (NOT DISABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT output)
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
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# --- Start: Logic for using Homebrew HTSlib, GSL, and LibDeflate on Apple Silicon ---
set(WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON FALSE)
set(WFMASH_USE_HOMEBREW_GSL_ON_APPLE_SILICON FALSE)
set(WFMASH_USE_HOMEBREW_LIBDEFLATE_ON_APPLE_SILICON FALSE)

if (APPLE AND (CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64"))
    # Check for HTSlib
    pkg_check_modules(HOMEBREW_HTSLIB QUIET htslib)
    if(HOMEBREW_HTSLIB_FOUND)
        set(WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON TRUE)
        message(STATUS "Apple Silicon: Found Homebrew HTSlib. Will use it.")
        message(STATUS "  Homebrew HTSlib Includes: ${HOMEBREW_HTSLIB_INCLUDE_DIRS}")
        message(STATUS "  Homebrew HTSlib Link Libraries: ${HOMEBREW_HTSLIB_LINK_LIBRARIES}")
    else()
        message(WARNING "Apple Silicon: Homebrew HTSlib not found via pkg-config. Will try other HTSlib options.")
    endif()

    # Check for GSL
    pkg_check_modules(HOMEBREW_GSL QUIET gsl)
    if(HOMEBREW_GSL_FOUND)
        set(WFMASH_USE_HOMEBREW_GSL_ON_APPLE_SILICON TRUE)
        message(STATUS "Apple Silicon: Found Homebrew GSL. Will use it.")
        message(STATUS "  Homebrew GSL Includes: ${HOMEBREW_GSL_INCLUDE_DIRS}")
        message(STATUS "  Homebrew GSL Link Libraries: ${HOMEBREW_GSL_LINK_LIBRARIES}")
    else()
        message(WARNING "Apple Silicon: Homebrew GSL not found via pkg-config. Will try other GSL options.")
    endif()

    # Check for LibDeflate
    pkg_check_modules(HOMEBREW_LIBDEFLATE QUIET libdeflate)
    if(HOMEBREW_LIBDEFLATE_FOUND)
        set(WFMASH_USE_HOMEBREW_LIBDEFLATE_ON_APPLE_SILICON TRUE)
        message(STATUS "Apple Silicon: Found Homebrew LibDeflate. Will use it.")
        message(STATUS "  Homebrew LibDeflate Includes: ${HOMEBREW_LIBDEFLATE_INCLUDE_DIRS}")
        message(STATUS "  Homebrew LibDeflate Link Libraries: ${HOMEBREW_LIBDEFLATE_LINK_LIBRARIES}")
    else()
        message(WARNING "Apple Silicon: Homebrew LibDeflate not found via pkg-config. Will try other LibDeflate options.")
    endif()
endif()
# --- End: Logic for using Homebrew HTSlib, GSL, and LibDeflate on Apple Silicon ---

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Release|Debug|Generic|RelWithDebInfo (default for distros)." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wno-unused-parameter")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic -Wno-unused-parameter")

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  SET(CMAKE_CXX_FLAGS_RELEASE  "")
  SET(CMAKE_C_FLAGS_RELEASE  "")
  set(EXTRA_FLAGS "-Ofast -march=native -funroll-all-loops")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG ${EXTRA_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG ${EXTRA_FLAGS}")
elseif (${CMAKE_BUILD_TYPE} MATCHES Generic)
  set(EXTRA_FLAGS "-O3 -funroll-all-loops")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG ${EXTRA_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG ${EXTRA_FLAGS}")
elseif (${CMAKE_BUILD_TYPE} MATCHES RelWithDebInfo)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
elseif (${CMAKE_BUILD_TYPE} MATCHES Debug)
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
  if (ASAN)
    set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-common")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -fno-common")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address")
  endif ()
else ()
    message(FATAL_ERROR "Unknown CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif ()

if(OpenMP_FOUND AND NOT (${CMAKE_BUILD_TYPE} MATCHES Debug))
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if( NOT DISABLE_LTO AND ipo_supported )
  message(STATUS "IPO / LTO enabled")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${ipo_supported})
endif()

if(PROFILER)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

if(GPROF)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -pg")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
endif(GPROF)

if (WFA_PNG_AND_TSV)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWFA_PNG_TSV_TIMING")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWFA_PNG_TSV_TIMING")
endif ()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Conditionally add HTSLIB submodule
if(NOT WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON AND NOT BUILD_DEPS)
    message(STATUS "HTSLIB: Using submodule from deps/htslib.")
    add_subdirectory(deps/htslib)
elseif(WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON)
    message(STATUS "Apple Silicon: Using Homebrew HTSlib, skipping build of deps/htslib submodule.")
else() 
    message(STATUS "HTSLIB: BUILD_DEPS is ON, htslib will be managed by ExternalProject_Add if defined below.")
endif()

add_subdirectory(deps/WFA2-lib)
set_target_properties(align_benchmark PROPERTIES EXCLUDE_FROM_ALL 1)
if (BUILD_STATIC)
  set_target_properties(wfa2cpp PROPERTIES EXCLUDE_FROM_ALL 1)
else()
  set_target_properties(wfa2cpp_static PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

add_subdirectory(src/common/wflign EXCLUDE_FROM_ALL)

if (BUILD_DEPS)
  include(ExternalProject)
  if(NOT WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON)
    ExternalProject_Add(htslib_ext
      URL https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/htslib_ext_build
      CONFIGURE_COMMAND autoreconf -i && ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/htslib_ext_install --disable-libcurl --disable-s3
      BUILD_COMMAND $(MAKE)
      INSTALL_COMMAND $(MAKE) install
      BUILD_IN_SOURCE 1
    )
    add_dependencies(wfmash htslib_ext)
  endif()

  if(NOT WFMASH_USE_HOMEBREW_GSL_ON_APPLE_SILICON)
    ExternalProject_Add(gsl_ext
      URL https://mirror.ibcp.fr/pub/gnu/gsl/gsl-2.8.tar.gz
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gsl_ext_build
      CONFIGURE_COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/gsl_ext_install
      BUILD_COMMAND $(MAKE)
      INSTALL_COMMAND $(MAKE) install
      BUILD_IN_SOURCE 1
    )
    add_dependencies(wfmash gsl_ext)
  endif()

  if(NOT WFMASH_USE_HOMEBREW_LIBDEFLATE_ON_APPLE_SILICON)
    ExternalProject_Add(libdeflate_ext
      URL https://github.com/ebiggers/libdeflate/releases/download/v1.20/libdeflate-1.20.tar.gz
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/libdeflate_ext_build
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/libdeflate_ext_install -DCMAKE_BUILD_TYPE=Release
    )
    add_dependencies(wfmash libdeflate_ext)
  endif()
endif()

add_executable(wfmash
  src/common/utils.cpp
  src/interface/main.cpp)

target_include_directories(wfmash PRIVATE
  src
  src/common
  src/common/wflign/deps
  deps/WFA2-lib
  deps/WFA2-lib/bindings/cpp
)

# Add Homebrew include directories if applicable
if(WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON)
    target_include_directories(wfmash PRIVATE ${HOMEBREW_HTSLIB_INCLUDE_DIRS})
endif()
if(WFMASH_USE_HOMEBREW_GSL_ON_APPLE_SILICON)
    target_include_directories(wfmash PRIVATE ${HOMEBREW_GSL_INCLUDE_DIRS})
endif()
if(WFMASH_USE_HOMEBREW_LIBDEFLATE_ON_APPLE_SILICON)
    target_include_directories(wfmash PRIVATE ${HOMEBREW_LIBDEFLATE_INCLUDE_DIRS})
endif()

# Linking logic for dependencies
if (BUILD_DEPS)
  # Includes for ExternalProjects (if not already covered by Homebrew logic for Apple Silicon)
  if(NOT WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON)
    target_include_directories(wfmash PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/htslib_ext_install/include)
  endif()
  if(NOT WFMASH_USE_HOMEBREW_GSL_ON_APPLE_SILICON)
     target_include_directories(wfmash PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/gsl_ext_install/include)
  endif()
  if(NOT WFMASH_USE_HOMEBREW_LIBDEFLATE_ON_APPLE_SILICON)
    target_include_directories(wfmash PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/libdeflate_ext_install/include)
  endif()

  # Linking for ExternalProjects
  if(NOT WFMASH_USE_HOMEBREW_GSL_ON_APPLE_SILICON)
    target_link_libraries(wfmash PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/gsl_ext_install/lib/libgsl.a
        ${CMAKE_CURRENT_BINARY_DIR}/gsl_ext_install/lib/libgslcblas.a
    )
  endif()
  if(NOT WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON)
    target_link_libraries(wfmash PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/htslib_ext_install/lib/libhts.a
    )
  endif()
  if(NOT WFMASH_USE_HOMEBREW_LIBDEFLATE_ON_APPLE_SILICON)
    target_link_libraries(wfmash PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}/libdeflate_ext_install/lib/libdeflate.a
    )
  endif()
else() # BUILD_DEPS is OFF
  # HTSlib
  if(WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON)
    message(STATUS "Apple Silicon (BUILD_DEPS OFF): Linking wfmash against Homebrew HTSlib using LINK_LIBRARIES.")
    target_link_libraries(wfmash PRIVATE ${HOMEBREW_HTSLIB_LINK_LIBRARIES})
  elseif(NOT WFMASH_USE_HOMEBREW_HTSLIB_ON_APPLE_SILICON AND TARGET htslib)
    message(STATUS "(BUILD_DEPS OFF): Linking wfmash against 'htslib' target (from submodule).")
    target_link_libraries(wfmash PRIVATE htslib)
  else()
    message(STATUS "(BUILD_DEPS OFF): Linking wfmash against system 'hts'.")
    target_link_libraries(wfmash PRIVATE hts) 
  endif()

  # GSL
  if(WFMASH_USE_HOMEBREW_GSL_ON_APPLE_SILICON)
      message(STATUS "Apple Silicon (BUILD_DEPS OFF): Linking wfmash against Homebrew GSL using LINK_LIBRARIES.")
      target_link_libraries(wfmash PRIVATE ${HOMEBREW_GSL_LINK_LIBRARIES})
  else()
      find_package(GSL REQUIRED) 
      if(GSL_FOUND)
          message(STATUS "(BUILD_DEPS OFF, Not Apple Silicon or Homebrew GSL not found): Linking wfmash against system GSL.")
          if(TARGET GSL::gsl AND TARGET GSL::gslcblas)
            target_link_libraries(wfmash PRIVATE GSL::gsl GSL::gslcblas)
          else()
            target_link_libraries(wfmash PRIVATE ${GSL_LIBRARIES})
          endif()
      else()
          message(FATAL_ERROR "System GSL not found, and not using Homebrew GSL. Please install GSL or enable BUILD_DEPS.")
      endif()
  endif()

  # LibDeflate
  if(WFMASH_USE_HOMEBREW_LIBDEFLATE_ON_APPLE_SILICON)
      message(STATUS "Apple Silicon (BUILD_DEPS OFF): Linking wfmash against Homebrew LibDeflate using LINK_LIBRARIES.")
      target_link_libraries(wfmash PRIVATE ${HOMEBREW_LIBDEFLATE_LINK_LIBRARIES})
  else()
      message(STATUS "(BUILD_DEPS OFF, Not Apple Silicon or Homebrew LibDeflate not found): Linking wfmash against system 'deflate'.")
      target_link_libraries(wfmash PRIVATE deflate)
  endif()
endif()

# Common system libraries for wfmash
target_link_libraries(wfmash PRIVATE
  m
  pthread
  # rt # rt is often Linux-specific, ensure it's needed/available on macOS
  ${LIBLZMA_LIBRARIES}
  ${BZIP2_LIBRARIES}
  ${ZLIB_LIBRARIES}
  Threads::Threads
)

# Link wflign and wfa2cpp (static or shared based on BUILD_STATIC)
if (BUILD_STATIC)
  target_link_libraries(wfmash PRIVATE wflign_static wfa2cpp_static)
else()
  target_link_libraries(wfmash PRIVATE wflign wfa2cpp)
endif()

# Link OpenMP using modern imported target
if(OpenMP_FOUND)
    target_link_libraries(wfmash PRIVATE OpenMP::OpenMP_CXX)
endif()

set(INVOKE ${CMAKE_BINARY_DIR}/bin/wfmash)

if(EXISTS ${CMAKE_SOURCE_DIR}/CTestCustom.cmake)
    configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR})
endif()

if (ASAN)
  set(INVOKE "ASAN_OPTIONS=detect_leaks=1:symbolize=1:verify_asan_link_order=0 LSAN_OPTIONS=verbosity=0:log_threads=1:suppressions=${CMAKE_SOURCE_DIR}/test/libasan-suppress.txt ${CMAKE_BINARY_DIR}/bin/wfmash")
else ()
  set(INVOKE "${CMAKE_BINARY_DIR}/bin/wfmash")
endif()

if (PROFILER)
  set(INVOKE "CPUPROFILE=${CMAKE_BINARY_DIR}/wfmash.prof ${INVOKE}")
endif()

set(REGRESSION_TEST_DIR "test/data/regression")

find_program(WGATOOLS "wgatools")
if (NOT WGATOOLS)
   set(WGATOOLS "env LD_LIBRARY_PATH=$GUIX_PROFILE/lib /wgatools/target/release/wgatools")
endif()
MESSAGE(STATUS "WGATOOLS found at: ${WGATOOLS}")

add_test(
  NAME wfmash-time-LPA
  COMMAND bash -c "${INVOKE} data/LPA.subset.fa.gz -p 80 -n 5 -t 4 > wfmash-time-LPA.output"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-subset-LPA-to-SAM
  COMMAND bash -c "${INVOKE} data/LPA.subset.fa.gz -N -a -L > LPA.subset.sam && samtools view LPA.subset.sam -bS | samtools sort > LPA.subset.bam && samtools index LPA.subset.bam && samtools view LPA.subset.bam | head | cut -f 1-9"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-mapping-coverage-with-8-yeast-genomes-to-PAF
  COMMAND bash -c "${INVOKE} data/scerevisiae8.fa.gz -p 95 -n 7 -m -L -Y \\# > scerevisiae8.paf && ./scripts/test.sh data/scerevisiae8.fa.gz.fai scerevisiae8.paf 0.89"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-pafcheck-yeast
  COMMAND bash -c "${INVOKE} data/scerevisiae8.fa.gz -t 4 -T S288C -Q Y12 > x.paf && pafcheck --query-fasta data/scerevisiae8.fa.gz --target-fasta data/scerevisiae8.fa.gz --paf x.paf"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-maf-validity
  COMMAND bash -c "${INVOKE} data/scerevisiae8.fa.gz -t 4 -T S288C -Q Y12 > x.paf && $(${WGATOOLS} paf2maf --target data/scerevisiae8.fa.gz --query data/scerevisiae8.fa.gz x.paf > x.maf) && test $(${WGATOOLS} stat x.maf | cut -f 7 | awk '{ s+=$1 } END { print s }') -gt 11000000"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
  NAME wfmash-multi-subset-index
  COMMAND bash -c "${INVOKE} data/scerevisiae8.fa.gz -t 4 -b 1m -T S288C -W index.idx && ${INVOKE} data/scerevisiae8.fa.gz -t 4 -T S288C -I index.idx -Q Y12 > index.paf && ./scripts/test.sh data/scerevisiae8.fa.gz.fai index.paf 0.9 'Y12\\|S288C'"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS wfmash DESTINATION bin)

if (BUILD_STATIC)
  install(TARGETS wfa2cpp_static LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(TARGETS wflign_static LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(TARGETS wfa2_static LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
else()
  install(TARGETS wfa2cpp LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(TARGETS wflign LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(TARGETS wfa2 LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/generate_git_version.sh ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/common/wflign/src)
