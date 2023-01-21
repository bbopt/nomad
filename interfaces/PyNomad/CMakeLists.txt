set(PYNOMAD_SOURCE
  ${CMAKE_CURRENT_SOURCE_DIR}/nomadCySimpleInterface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PyNomad.pyx
  ${CMAKE_CURRENT_SOURCE_DIR}/setup_PyNomad.py
)

# Modern-ish approach that does not (openly) complain about being deprecated.
# The build process supports both GCC and MSVC.
#
# NOTES (@FIXME)
# - Uses 'python -m pip wheel' which requires 'setup.py' to exist. 
# - Since we use 'setup_PyNomadStatic.py' instead, we must make a copy first.

set(PyNomad_MSVC_VIA "")

if (MSVC)
  set(PyNomad_MSVC_VIA "YES")
endif()

add_custom_target(
  PyNomad
  COMMAND
    ${CMAKE_COMMAND} -E copy ./setup_PyNomadStatic.py setup.py
  COMMAND
  ${CMAKE_COMMAND} -E env "NOMAD_INSTALL=${CMAKE_BINARY_DIR}" "NOMAD_VERSION=${NOMAD_VERSION}" "NOMAD_VIAMSVC=${PyNomad_MSVC_VIA}" python -m pip wheel -w dist -- .
  DEPENDS 
    nomadStatic
    ${PYNOMAD_SOURCE}
  WORKING_DIRECTORY 
    ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)

