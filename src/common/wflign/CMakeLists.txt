#
# CMAKE
# 
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wflign VERSION 0.0.1)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

#
# Config
# 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcx16 -march=native -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcx16 -march=native -g")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -mcx16 -march=native -g -fsanitize=address")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -mcx16 -march=native -g -fsanitize=address")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#
# WFA Build
# 
add_custom_target(
   wfa_lib
   COMMAND make clean all
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/deps/WFA
)

#
# WFling
# 
add_executable(
  wflign
  src/main.cpp
  src/wflign.cpp
  src/wflign_alignment.cpp
  src/wflign_patch.cpp
  src/rkmh.cpp
  src/murmur3.cpp)
target_include_directories(
  wflign PRIVATE deps/
  deps/ 
  deps/robin-hood-hashing
  deps/WFA)

#
# WFling library
# 
target_link_libraries(wflign z)
add_library(
  libwflign_static STATIC
  src/wflign.cpp
  src/wflign_alignment.cpp
  src/wflign_patch.cpp
  src/rkmh.cpp
  src/murmur3.cpp)
set_target_properties(libwflign_static PROPERTIES OUTPUT_NAME "wflign")
target_include_directories(
  libwflign_static PRIVATE
  deps/
  deps/robin-hood-hashing
  deps/WFA)
add_dependencies(libwflign_static wfa_lib)

