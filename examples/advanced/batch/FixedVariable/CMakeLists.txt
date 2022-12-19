set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(uu.exe uu.cpp )
set_target_properties(uu.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS uu.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    if (WIN32)
        add_test(NAME ExampleAdvancedBatchFixedVariables
        COMMAND bash.exe runFixedWin.sh ${CMAKE_INSTALL_PREFIX}/bin/nomad.exe
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    else()
        message(STATUS "    Add example batch fixed variables")

        add_test(NAME ExampleAdvancedBatchFixedVariables
        COMMAND ./runFixed.sh ${CMAKE_INSTALL_PREFIX}/bin/nomad
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    endif()
endif()
