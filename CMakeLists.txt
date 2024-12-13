# Usage:
#
#   mkdir build ; cd build
#   cmake -DCMAKE_BUILD_TYPE=DEBUG ..
#   make -j 8 VERBOSE=1
#

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(wfmash)

enable_testing()

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard is not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(BUILD_STATIC "Build static binary" OFF)
option(BUILD_DEPS "Build external dependencies" OFF)
option(BUILD_OPTIMIZED "Build optimized binary" OFF)
option(GPROF "Enable gprof profiling" OFF)
option(DISABLE_LTO "Disable IPO/LTO" OFF)

if (NOT DISABLE_LTO)
  include(CheckIPOSupported) # adds lto
  check_ipo_supported(RESULT ipo_supported OUTPUT output) # lto
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
            "Choose the type of build, options are: Release Debug Generic." FORCE)
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
    # Enable debug symbols and ASan with no optimizations
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -fsanitize=address -fno-omit-frame-pointer")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fsanitize=address -fno-omit-frame-pointer")
    # Ensure that ASan is linked explicitly
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
else ()
    set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif ()

if( NOT DISABLE_LTO AND ipo_supported )
  message(STATUS "IPO / LTO enabled")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${ipo_supported})
endif()

if(GPROF)
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

add_test(
  NAME wfmash-test
  COMMAND wfmash data/LPA.subset.fa.gz -p 80 -n 5 -t 8
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
