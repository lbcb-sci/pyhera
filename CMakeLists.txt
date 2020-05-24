cmake_minimum_required (VERSION 3.2)
project (ScaRa)

# The version number.
set (ScaRa_VERSION_MAJOR 1)
set (ScaRa_VERSION_MINOR 0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

# Add libraries
# add_library(sbridger src/SBridger.cpp)
# add_library(loader src/Loader.cpp)
# add_library(overlap src/Overlap.cpp)
# add_library(sequence src/Sequence.cpp)

set(SOURCE_FILES src/scara.cpp src/Sequence.cpp src/Overlap.cpp src/Loader.cpp src/SBridger.cpp 
src/GraphBasic.cpp src/GraphAlgorithms.cpp)

# Add the executable
add_executable(scara ${SOURCE_FILES})

# Adding bioparser
add_subdirectory(ezra/vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool)
# target_link_libraries(scara bioparser sbridger loader overlap sequence)
target_link_libraries(scara bioparser thread_pool)
