set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(uu.exe uu.cpp )
set_target_properties(uu.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS uu.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)

    message(STATUS "    Add example batch fixed variables") 
    # Test run in working directory AFTER install of bb executable

    if (WIN32)
        add_test(NAME ExampleAdvancedBatchFixedVariables
        COMMAND bash.exe runFixedWin.sh ${CMAKE_INSTALL_PREFIX}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    # Temp comment
    #else()
    #    add_test(NAME ExampleAdvancedBatchFixedVariables
    #    COMMAND ./runFixed.sh ${CMAKE_INSTALL_PREFIX}
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    endif()
endif()
