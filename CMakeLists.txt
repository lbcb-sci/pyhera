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

set(SOURCE_FILES_SCARA src/scara.cpp src/Sequence.cpp src/Overlap.cpp src/Loader.cpp src/SBridger.cpp 
src/GraphBasic.cpp src/GraphAlgorithms.cpp)
set(SOURCE_FILES_LOAD1 src/Load_HiC.cpp src/Sequence.cpp src/Overlap.cpp src/Loader.cpp)
# Adding bioparser
add_subdirectory(ezra/vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool)

# Add the scara executable
add_executable(scara ${SOURCE_FILES_SCARA})
# target_link_libraries(scara bioparser sbridger loader overlap sequence)
target_link_libraries(scara bioparser thread_pool)

# Add the load_HiC executable
add_executable(load_HiC ${SOURCE_FILES_LOAD1})
target_link_libraries(load_HiC bioparser)
