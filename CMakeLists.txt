# CMake project
cmake_minimum_required(VERSION 2.8)
project("vmath" C CXX)

# Options
option(ENABLE_SAMPLES "Enable compilation of sample programs" OFF)
option(ENABLE_TESTS "Enable compilation of the Google tests" ON)



# Compiler-version check
#if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
#    message(FATAL_ERROR "GCC version must be at least 4.7 (for C++11)")
#  endif()
#endif()

# C++11 compiler settings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wall -mavx2")

# Package scripts location
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# The includes
include_directories(
  ${vmath_SOURCE_DIR}/include
  )

if (ENABLE_SAMPLES)
endif()

if (ENABLE_TESTS)
  add_subdirectory(external/gtest)
  enable_testing()
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

  add_executable(unit_tests 
    unit_tests/simdmath.cpp
    unit_tests/simdmem.cpp
    )

  target_link_libraries(unit_tests gtest gtest_main)

  # Adds the tests
  add_test(
    name unit_tests
    command unit_tests
    )




endif()

