cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wflign VERSION 0.0.1)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mcx16 -march=native -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mcx16 -march=native -g")
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -mcx16 -march=native -g -fsanitize=address")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -mcx16 -march=native -g -fsanitize=address")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(deps/wflambdav2 EXCLUDE_FROM_ALL)
add_subdirectory(deps/WFAv2 EXCLUDE_FROM_ALL)

add_executable(wflign
        deps/atomic_image.cpp
        deps/lodepng/lodepng.cpp
        deps/atomic_image.cpp
  src/wflign_wfa.cpp
  src/wfa_edit_callback.cpp
  src/rkmh.cpp
  src/murmur3.cpp)

target_include_directories(wflign PRIVATE
  deps/)

target_link_libraries(wflign
        lodepng
        picosha256
  wfa
  wflambda
  z
atomic_image
  )

add_library(libwflign_static STATIC
        deps/atomic_image.cpp
  deps/lodepng/lodepng.cpp
  deps/atomic_image.cpp
        src/wflign_wfa.cpp
        src/wfa_edit_callback.cpp
        src/rkmh.cpp
        src/murmur3.cpp
  )
set_target_properties(libwflign_static PROPERTIES OUTPUT_NAME "wflign")
target_include_directories(libwflign_static PRIVATE
  deps/
  deps/lodepng
        deps/PicoSHA2
  deps/robin-hood-hashing
  deps/wflambdav2
  deps/WFAv2
  )
# not clear if necessary
#set_target_properties(libodgi_static PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

#set(CMAKE_BUILD_TYPE Release)
