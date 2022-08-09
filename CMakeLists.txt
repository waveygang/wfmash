cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wfmash VERSION 0.0.1)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE) # does not work

#set(CMAKE_BUILD_TYPE Release)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Release Debug Generic." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if (${CMAKE_BUILD_TYPE} MATCHES Release)
    #set(EXTRA_FLAGS "-Ofast -march=native -flto -fno-fat-lto-objects")
    set(EXTRA_FLAGS "-Ofast -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    # Debug use the defaults - so commenting out:
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -g -fsanitize=address")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -g -fsanitize=address")
else()
    # Use all standard-compliant optimizations - always add these:
    set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
    set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
endif ()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#add_subdirectory(src/common/wflign/deps/WFA2-lib EXCLUDE_FROM_ALL)
#add_subdirectory(src/common/wflign/deps/wflambda EXCLUDE_FROM_ALL)
add_subdirectory(src/common/wflign EXCLUDE_FROM_ALL)
#add_subdirectory(src/common/abseil-cpp)

add_executable(wfmash
  src/common/utils.cpp
  src/yeet/yeet_main.cpp
  src/version.cpp) #$(CXX) $(CXXFLAGS) $(CPPFLAGS) -Lsrc/common/wflign/build/lib -Isrc/common/wflign/deps -Isrc/common $(SOURCE_1) -o wfmash @mathlib@ -lstdc++ -lz -lm -lpthread -lwflambda -lwflign

target_include_directories(wfmash PRIVATE
  src
  src/common
  src/common/wflign/deps
  src/common/wflign/deps/WFA2-lib
  #src/common/wflign/abseil-cpp
  )
target_link_directories(wfmash PRIVATE
 src/common/wflign/deps/WFA2-lib/lib)
target_link_libraries(wfmash
  gsl
  gslcblas
  m
  pthread
  libwflign_static
  hts
  rt
  wfacpp
  jemalloc
  z
  #absl::flat_hash_map
)

install(TARGETS wfmash DESTINATION bin)

# version stuff
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/generate_git_version.sh ${CMAKE_SOURCE_DIR}/src)
