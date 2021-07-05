set(PYNOMAD_SOURCE
${CMAKE_CURRENT_SOURCE_DIR}/nomadCySimpleInterface.cpp
${CMAKE_CURRENT_SOURCE_DIR}/PyNomad.pyx
${CMAKE_CURRENT_SOURCE_DIR}/setup_PyNomad.py
)


# Command to create PyNomad interface
if(OpenMP_CXX_FOUND)
    set(USE_OPENMP 1)
else()
    set(USE_OPENMP 0)
endif()

add_custom_command(
    OUTPUT PyNomad.so
    COMMAND python setup_PyNomad.py ${USE_OPENMP} ${CMAKE_BINARY_DIR} build_ext --inplace
    COMMAND python setup_PyNomad.py ${USE_OPENMP} ${CMAKE_BINARY_DIR} install --user
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS ${PYNOMAD_SOURCE}
)

# target pynomad is always built
add_custom_target(pynomad ALL
    # COMMAND echo "This is ALL target 'pynomad', and it depends on ${PYNOMAD_SOURCE}"
    COMMAND echo "IMPORTANT: need to run cmake --install before using PyNomad module"
    DEPENDS PyNomad.so nomadUtils nomadEval nomadAlgos
    VERBATIM
)

