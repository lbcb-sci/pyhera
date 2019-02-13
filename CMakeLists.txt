cmake_minimum_required (VERSION 3.2)
project (ScaRa)

# The version number.
set (BSckktest_VERSION_MAJOR 1)
set (BSckktest_VERSION_MINOR 0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
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

# Add the executable
add_executable(scara src/ScaRa.cpp)

# Adding bioparser
add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
target_link_libraries(scara bioparser)

# Adding Ezra
add_subdirectory(Ezra EXCLUDE_FROM_ALL)
target_link_libraries(scara Ezra)
