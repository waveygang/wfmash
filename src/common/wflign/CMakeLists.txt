cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(wflign VERSION 0.0.1)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#set(CMAKE_BUILD_TYPE Release)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Release Debug." FORCE)
endif()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build.")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
    message("Release build.")
ELSE()
    message("Some other build type.")
ENDIF()

# set(EXTRA_FLAGS "-Ofast -mcx16 -march=native ")

if (${CMAKE_BUILD_TYPE} MATCHES "Release")
    if (EXTRA_FLAGS)
        set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
    endif (EXTRA_FLAGS)
    # Use all standard-compliant optimizations - always add these:
    set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
    set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    # Debug use the defaults - so commenting out:
    # set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -g -fsanitize=address")
    # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -g -fsanitize=address")
endif ()

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
  wfa
  wflambda
  z
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
