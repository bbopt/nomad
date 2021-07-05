# ----------------------------------------------------------------------------
#  Root CMake file for sgtelib library
#
#    From the off-tree build directory, invoke:
#      $ cmake <PATH_TO_OPENCV_ROOT>
#
# ----------------------------------------------------------------------------
# cmake version : the latest one
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# name of the project
project (sgtelib LANGUAGES CXX VERSION 2.0.2)

# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed.
    You should create a separate directory for build files.")
endif()

# check compiler version
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4)
        message(FATAL_ERROR "GCC version < 4 has not been tested for Nomad")
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

# require c++14
set(CXX_STANDARD_REQUIRED 14)


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
    src/sgtelib.hpp
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
    src/sgtelib.cpp
    src/sgtelib_help.cpp
    )

# load classic directories for unix systems
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# build shared library sgtelib
add_library (sgtelib SHARED ${SGTELIB_SOURCES})
target_include_directories(sgtelib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/sgtelib>
    )

# export symbols on Windows
if (WIN32)
    target_compile_definitions(sgtelib PRIVATE -DDLL_EXPORTS)
endif()

## build executable sgtelib_main
#add_executable(sgtelib_main ${SGTELIB_SOURCES})
#target_include_directories(sgtelib_main PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
#    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/sgtelib>
#    )

## installing executables and libraries
#install(TARGETS sgtelib_main sgtelib 
#    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sgtelib
#    )

# installing executables and libraries
install(TARGETS sgtelib 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/sgtelib
    )


# install header files
install(FILES ${SGTELIB_HEADERS}
    DESTINATION include/sgtelib
    )




