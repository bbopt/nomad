# ----------------------------------------------------------------------------
#  Root CMake file for sgtelib library
#
#    From the off-tree build directory, invoke:
#      $ cmake <PATH_TO_OPENCV_ROOT>
#
# ----------------------------------------------------------------------------
# cmake version : the latest one
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set(SGTELIB_VERSION_MAJOR 2)
set(SGTELIB_VERSION_MINOR 0)
set(SGTELIB_VERSION_PATCH 4)
set(SGTELIB_VERSION ${SGTELIB_VERSION_MAJOR}.${SGTELIB_VERSION_MINOR}.${SGTELIB_VERSION_PATCH})

# name of the project
project (sgtelib LANGUAGES CXX VERSION ${SGTELIB_VERSION})

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed.
    You should create a separate directory for build files.")
endif()

# check compiler version
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 8
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
        message(FATAL_ERROR "GCC version < 8 has not been tested for sgtelib")
    endif()
elseif (("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") 
    OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
    # require at least clang 5
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
        message(FATAL_ERROR "Clang version has not been tested for sgtelib")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # require at least 15.8 (MSVC 2017) for C++17 support
    if (MSVC_TOOLSET_VERSION VERSION_LESS 141)
        message(FATAL_ERROR "MSVC version ${CMAKE_CXX_COMPILER_VERSION} has not been tested for sgtelib")
    endif()
else()
    message(WARNING "You are using an unsupported compiler. Compilation has only been tested with Clang, GCC, and MSVC.")
endif()

# require c++17
set(CXX_STANDARD_REQUIRED 17)


# header files
set(SGTELIB_HEADERS
    src/Defines.hpp
    src/Exception.hpp
    src/Kernel.hpp
    src/Matrix.hpp
    src/Metrics.hpp
    src/Surrogate.hpp
    src/Surrogate_CN.hpp
    src/Surrogate_Ensemble.hpp
    src/Surrogate_Ensemble_Stat.hpp
    src/Surrogate_Factory.hpp
    src/Surrogate_KS.hpp
    src/Surrogate_Kriging.hpp
    src/Surrogate_LOWESS.hpp
    src/Surrogate_PRS.hpp
    src/Surrogate_PRS_CAT.hpp
    src/Surrogate_PRS_EDGE.hpp
    src/Surrogate_Parameters.hpp
    src/Surrogate_RBF.hpp
    src/Surrogate_Utils.hpp
    src/Tests.hpp
    src/TrainingSet.hpp
    src/sgtelib_help.hpp
    )

# source files
set(SGTELIB_SOURCES
    src/Kernel.cpp
    src/Matrix.cpp
    src/Metrics.cpp
    src/Surrogate.cpp
    src/Surrogate_CN.cpp
    src/Surrogate_Ensemble.cpp
    src/Surrogate_Ensemble_Stat.cpp
    src/Surrogate_Factory.cpp
    src/Surrogate_KS.cpp
    src/Surrogate_Kriging.cpp
    src/Surrogate_LOWESS.cpp
    src/Surrogate_PRS.cpp
    src/Surrogate_PRS_CAT.cpp
    src/Surrogate_PRS_EDGE.cpp
    src/Surrogate_Parameters.cpp
    src/Surrogate_RBF.cpp
    src/Surrogate_Utils.cpp
    src/Tests.cpp
    src/TrainingSet.cpp
    src/sgtelib_help.cpp
    )

set(SGTELIB_MAIN_SOURCE src/sgtelib.cpp)

# load classic directories for unix systems
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Build shared libraries (sgtelib)
include(./CMakeListsLibs.txt)

# Build static libraries (sgtelibStatic)
if(BUILD_INTERFACE_PYTHON MATCHES ON)
  message(STATUS "  Warning: build sgtelib static libraries for Python interface.")
  include(./CMakeListsLibsStatic.txt)
endif()



# install header files
install(
  FILES ${SGTELIB_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sgtelib
)




