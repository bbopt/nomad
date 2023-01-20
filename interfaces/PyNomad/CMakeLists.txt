set(PYNOMAD_SOURCE
${CMAKE_CURRENT_SOURCE_DIR}/nomadCySimpleInterface.cpp
${CMAKE_CURRENT_SOURCE_DIR}/PyNomad.pyx
${CMAKE_CURRENT_SOURCE_DIR}/setup_PyNomad.py
)

if (WITH_STATIC_BUILD MATCHES ON)

  # Modern-ish approach that does not (openly) complain about being deprecated.
  # The build process supports both GCC and MSVC.
  #
  # NOTES
  # - Uses 'python -m pip wheel' which requires 'setup.py' to exist. 
  # - Since we use 'setup_PyNomadStatic.py' instead, we must make a copy first.

  set(PyNomadStatic_MSVC_VIA "")

  if (MSVC)
    set(PyNomadStatic_MSVC_VIA "YES")
  endif()

  add_custom_target(
    PyNomadStatic ALL
    COMMAND
      ${CMAKE_COMMAND} -E copy ./setup_PyNomadStatic.py setup.py
    COMMAND
      ${CMAKE_COMMAND} -E env "NOMAD_SRC=${CMAKE_SOURCE_DIR}/src" "NOMAD_BUILD_DIR=${CMAKE_BINARY_DIR}" "NOMAD_VERSION=${NOMAD_VERSION}" "NOMAD_VIAMSVC=${PyNomadStatic_MSVC_VIA}" python -m pip wheel -w dist -- .
    DEPENDS 
      nomadStatic
      ${PYNOMAD_SOURCE}
    WORKING_DIRECTORY 
      ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
)

else()

# Command to create PyNomad interface
add_custom_command(
    OUTPUT PyNomad.${NOMAD_VERSION}.so # Test without it (not plateform indep.)
    COMMAND python setup_PyNomad.py ${CMAKE_BINARY_DIR} ${NOMAD_VERSION} build_ext --inplace
    COMMAND python setup_PyNomad.py ${CMAKE_BINARY_DIR} ${NOMAD_VERSION} install --user
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS ${PYNOMAD_SOURCE}
)

# target pynomad is always built
add_custom_target(pynomad ALL
    COMMAND echo "IMPORTANT: need to run cmake --install before using PyNomad module"
    DEPENDS PyNomad.${NOMAD_VERSION}.so nomadUtils nomadEval nomadAlgos
    VERBATIM
)
endif()

