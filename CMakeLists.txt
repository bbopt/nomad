# ----------------------------------------------------------------------------
#  Root CMake file for Nomad
# ----------------------------------------------------------------------------

# cmake version : the latest one
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# name of the project
project (Nomad VERSION 4.0.0)

# use standard compilers parameters for c++14
SET(CMAKE_CXX_STANDARD 14 )
SET(CMAKE_CXX_STANDARD_REQUIRED ON )

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed.
    You should create a separate directory for build files.
    ")
endif()

#check compiler version
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
        message(FATAL_ERROR "GCC version < 8 has not been tested for Nomad!")
    endif()
elseif (("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") 
    OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
    # require at least clang 5
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "Clang version has not been tested for Nomad!")
    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()

#options

# choose to build with unit tests
# option(BUILD_TESTS "Build unit tests programs" OFF)

# choose to build with benchmarks
# option(BUILD_BENCHMARKS "Build benchmark tests" OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# special flag required
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# add sgtelib
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext)

#add nomad directory
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

#enable testing
#if (BUILD_TESTS)
#    enable_testing()
#	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_tests)
#endif(BUILD_TESTS)
