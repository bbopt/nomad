set(PYNOMAD_SOURCE
${CMAKE_CURRENT_SOURCE_DIR}/nomadCySimpleInterface.cpp
${CMAKE_CURRENT_SOURCE_DIR}/PyNomad.pyx
${CMAKE_CURRENT_SOURCE_DIR}/setup_PyNomad.py
)


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

