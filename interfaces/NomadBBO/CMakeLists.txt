set(NOMADBBO_SOURCE
${CMAKE_CURRENT_SOURCE_DIR}/nomadCySimpleInterface.cpp
${CMAKE_CURRENT_SOURCE_DIR}/NomadBBO_src.pyx
${CMAKE_CURRENT_SOURCE_DIR}/setup.py
)


# Detect sysroot for macOS, otherwise leave empty
if(APPLE)
    set(NOMADBBO_SYSROOT ${CMAKE_OSX_SYSROOT})
else()
    set(NOMADBBO_SYSROOT "")
endif()

# Detect architecture
if(CMAKE_OSX_ARCHITECTURES)
    set(PY_ARCH ${CMAKE_OSX_ARCHITECTURES})
elseif(CMAKE_SYSTEM_PROCESSOR)
    set(PY_ARCH ${CMAKE_SYSTEM_PROCESSOR})
else()
    set(PY_ARCH "unknown")
endif()


add_custom_target(
    NomadBBO ALL
    COMMAND
    ${CMAKE_COMMAND} -E env
      "CC=${CMAKE_C_COMPILER}"
      "CXX=${CMAKE_CXX_COMPILER}"
      "SYSROOT=${NOMADBBO_SYSROOT}"
      "NOMAD_SRC=${CMAKE_SOURCE_DIR}/src"
      "NOMAD_BUILD_DIR=${CMAKE_BINARY_DIR}"
      "NOMAD_MSVC_FLAG=${MSVC}"
      "NOMAD_MSVC_CONF=$<CONFIG>"
      "BUILD_OPENMP=${OpenMP_FOUND}"
      ARCHFLAGS=-arch\ ${PY_ARCH}
      ${Python_EXECUTABLE} -m pip wheel -w dist -- .
    DEPENDS
      ${NOMADBBO_SOURCE}
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
)

# NOMAD must be built before building the NomadBBO wheel
add_dependencies(NomadBBO nomadStatic)
