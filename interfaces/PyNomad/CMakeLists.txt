set(PYNOMAD_SOURCE
${CMAKE_CURRENT_SOURCE_DIR}/nomadCySimpleInterface.cpp
${CMAKE_CURRENT_SOURCE_DIR}/PyNomad.pyx
${CMAKE_CURRENT_SOURCE_DIR}/setup.py
)


# Modern-ish approach that does not (openly) complain about being deprecated.
# The build process supports both GCC and MSVC.

add_custom_target(
    PyNomad ALL
    COMMAND
    ${CMAKE_COMMAND} -E env "NOMAD_SRC=${CMAKE_SOURCE_DIR}/src" "NOMAD_BUILD_DIR=${CMAKE_BINARY_DIR}" "NOMAD_MSVC_FLAG=${MSVC}" "NOMAD_MSVC_CONF=$<CONFIG>" ${Python_EXECUTABLE} -m pip wheel -w dist -- .
    DEPENDS 
      nomadStatic
      ${PYNOMAD_SOURCE}
    WORKING_DIRECTORY 
      ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
)
