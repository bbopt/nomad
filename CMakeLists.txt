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


# choose to build with benchmarks
# option(BUILD_BENCHMARKS "Build benchmark tests" OFF)

# Load some modules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# special flag required
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#
# Message for starting configuration
#
message(CHECK_START " Configuring custom options")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

#
# Modify the install prefix if not specified on the command
#
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE PATH "..." FORCE) 
message(STATUS "  Installation prefix set to ${CMAKE_INSTALL_PREFIX}")

#
# Use sgtelib
#
set(USE_SGTELIB ON)
if(USE_SGTELIB MATCHES ON)
  message(CHECK_START "  Configuring for use of Sgtelib library")

  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ext/sgtelib) 
     message(CHECK_FAIL "  failed. Directory ${CMAKE_CURRENT_SOURCE_DIR}/ext/sgtelib not found. Let's configure without Sgtelib library")
     set(USE_SGTELIB OFF)
  else()
     add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext)
     add_compile_definitions(USE_SGTELIB)
     list(APPEND EXTRA_LIBS sgtelib)
     message(CHECK_PASS "   Sgtelib found. done")
  endif()
else()
   message(STATUS "  Sgtelib library will NOT be used\n")
endif()

#
# Choose to build with time stats enabled
#
option(TIME_STATS "Option to build with time stats enabled" OFF)
if(TIME_STATS MATCHES ON)
   message(CHECK_START "  Enabling time stats for build")
   add_compile_definitions(TIME_STATS)
endif()

#
# Choose to build with unit tests or not
#
option(BUILD_TESTS "Option to build the unit tests" OFF)
if(BUILD_TESTS MATCHES ON)
   message(CHECK_START "  Configuring unit tests build")

   if(NOT DEFINED ENV{GTEST_DIR} OR NOT EXISTS $ENV{GTEST_DIR})
        message(FATAL_ERROR "The NOMAD 4 unit tests need the Google Test framework. Please set environment variable GTEST_DIR to your local googletest directory. The Google Test framework can be downloaded from https://github.com/google/googletest")
   endif()

   enable_testing()

   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_tests)
   message(CHECK_PASS " done")
else()
   message(STATUS "  Unit tests NOT built")
endif()

#
# Choose to build the examples/tests using library mode or not
#
option(BUILD_LIBMODE_EXAMPLES "Option to build library mode examples" ON)
if(BUILD_LIBMODE_EXAMPLES MATCHES ON)
   message(CHECK_START "  Configuring build for library mode examples and tests")

   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/integration_tests)
   message(CHECK_PASS " done")
else()
   message(STATUS "  Examples/tests NOT built")
endif()

#
# Choose to build the interface
#
option(BUILD_INTERFACES "Option to build interfaces to Nomad" OFF)
if(BUILD_INTERFACES MATCHES ON)
   message(CHECK_START "  Configuring build for interfaces")

   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/interfaces)
   message(CHECK_PASS " done")
else()
   message(STATUS "  Interfaces to Nomad NOT built")
endif()

#
# Test openMP package
#
find_package(OpenMP QUIET)
if(OpenMP_FOUND)
   message(STATUS "  Test OpenMP for parallel functionalities  -- found ")
else()
   message(STATUS "  OpenMP not found. Parallel functionalities NOT available")
endif()

#
# Custom options final message
#
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS " done")

#
# Add nomad app directory for building 
#
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)

