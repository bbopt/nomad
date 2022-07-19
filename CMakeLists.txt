# ----------------------------------------------------------------------------
#  Root CMake file for Nomad
# ----------------------------------------------------------------------------

# cmake version : the latest one
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# name of the project
project (Nomad VERSION 4.2.0)

# Nomad src template (copy of template defined in src/nomad_version.hpp, used by jNomad) 
# Need to update when version changes
set(NOMAD_SRC_TEMPLATE NOMAD_4_2)

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
        message(FATAL_ERROR "GCC version < 8 has not been tested for Nomad")
    endif()
elseif (("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") 
    OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
    # require at least clang 5
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "Clang version has not been tested for Nomad")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # require at least 15.0 (MSVC 2017) for C++14 support
    if (MSVC_TOOLSET_VERSION VERSION_LESS 141)
        message(FATAL_ERROR "MSVC version ${CMAKE_CXX_COMPILER_VERSION} has not been tested for Nomad")
    endif()
else()
    message(WARNING "You are using an unsupported compiler. Compilation has only been tested with Clang, GCC, and MSVC.")
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
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE PATH "..." FORCE) 
endif()
message(STATUS "  Installation prefix set to ${CMAKE_INSTALL_PREFIX}")

#
# Use sgtelib
#
option(USE_SGTELIB "Option to use the Sgtelib library" ON)
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
# Choose to disable display degree tests for faster code
#
option(REDUCED_DISPLAY_TEST "Option to disable display tests for faster run" OFF)
if(REDUCED_DISPLAY_TEST MATCHES ON)
   message(CHECK_START "  Disable display test for build")
   add_compile_definitions(REDUCED_DISPLAY_TEST)
endif()

#
# Test openMP package
#
option(TEST_OPENMP "Option to compile Nomad with OpenMP" ON)
if(TEST_OPENMP MATCHES ON)
   find_package(OpenMP QUIET)
   if(OpenMP_FOUND)
      message(STATUS "  Test OpenMP for parallel functionalities  -- found ")
   else()
      message(STATUS "  OpenMP not found. Parallel functionalities NOT available")
   endif()
else()
    message(STATUS "  Compiling Nomad without OpenMP")
endif()

#
# Choose to build the examples/tests using batch and library mode 
#
option(BUILD_EXAMPLES "Option to build library mode examples" ON)
if(BUILD_EXAMPLES MATCHES ON)
   message(CHECK_START "  Configuring build for library and batch mode examples")

   enable_testing()

   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
   message(CHECK_PASS " done")
else()
   message(STATUS "  Examples NOT built")
endif()


#
# Choose to build the C interface
#
option(BUILD_INTERFACE_C "Option to build C interfaces to Nomad" OFF)
if(BUILD_INTERFACE_C MATCHES ON)
   message(CHECK_START "  Configuring build for C interface")
   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/CInterface)
   message(CHECK_PASS " done")
else()
   message(STATUS "  C interface to Nomad NOT built")
endif()

#
# Choose to build the Python interface
#
option(BUILD_INTERFACE_PYTHON "Option to build Python interface to Nomad" OFF)
if(BUILD_INTERFACE_PYTHON MATCHES ON)
   message(CHECK_START "  Configuring build for Python interface")
   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/PyNomad)
   message(CHECK_PASS " done")
else()
   message(STATUS "  Python interface to Nomad NOT built")
endif()

#
# Choose to build the Matlab interface
#
option(BUILD_INTERFACE_MATLAB "Option to build Matlab interface to Nomad" OFF)
if(BUILD_INTERFACE_MATLAB MATCHES ON)
   if(OpenMP_FOUND)
       message(STATUS "  Warning: Cannot build Matlab interface with OpenMP enabled")
   else()
       message(CHECK_START "  Configuring build for Matlab interface")
       find_package(Matlab QUIET REQUIRED COMPONENTS MEX_COMPILER)
       if(MATLAB_FOUND)
          add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/Matlab_MEX)
          message(CHECK_PASS " done")
       else()
          message(STATUS "  Matlab no found... interface not built")
       endif()
   endif()
else()
   message(STATUS "  Matlab interface to Nomad NOT built")
endif()

#
# Choose to build the Java interface
#
option(BUILD_INTERFACE_JAVA "Option to build Java Swig interface to Nomad" OFF)
if(BUILD_INTERFACE_JAVA MATCHES ON)
   find_package(Java REQUIRED COMPONENTS Runtime Development)
   find_package(SWIG 4.0 REQUIRED COMPONENTS java)
   find_package(JNI REQUIRED)
   if(SWIG_FOUND)
       if(NOT SWIG_java_FOUND)
           message(WARNING "   SWIG java bindings cannot be generated")
       else()
           add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/interfaces/jNomad)
       endif()
   else()
       message(STATUS "  Swig not found... java interface not built")
   endif()

else()
   message(STATUS "  Java Swig interface to Nomad NOT built")
endif()


#
# Custom options final message
#
list(REMOVE_ITEM CMAKE_MESSAGE_INDENT "  ")
message(CHECK_PASS " done")

#
# Add nomad app directory for building 
#
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
