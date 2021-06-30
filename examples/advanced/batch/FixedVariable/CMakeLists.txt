set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(uu uu.cpp )

# installing executables and libraries
install(TARGETS uu
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)

    message(STATUS "    Add example batch test fixed variables") 
    # Test run in working directory AFTER install of bb executable

    if (WIN32)
        add_test(NAME ExampleAdvancedBatchFixedVariables
        COMMAND bash.exe runFixedWin.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    else()
        add_test(NAME ExampleAdvancedBatchFixedVariables
        COMMAND ./runFixed.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    endif()
endif()
