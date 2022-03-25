cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wfmash VERSION 0.0.1)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcx16 -march=native -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcx16 -march=native -g")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -mcx16  -march=native -g -fsanitize=address")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -mcx16  -march=native -g -fsanitize=address")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#add_subdirectory(src/common/wflign/deps/WFA EXCLUDE_FROM_ALL)
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
  src/common/wflign/deps/wflambdav2
  src/common/wflign/deps/WFAv2
  #src/common/wflign/abseil-cpp
  )

target_link_libraries(wfmash
  gsl
  gslcblas
  m
  pthread
  libwflign_static
  hts
  wflambda
  wfa
  jemalloc
  z
  #absl::flat_hash_map
)

install(TARGETS wfmash DESTINATION bin)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# version stuff
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/generate_git_version.sh ${CMAKE_SOURCE_DIR}/src)
