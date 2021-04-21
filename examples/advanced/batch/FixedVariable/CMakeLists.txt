add_executable(uu.exe uu.cpp )

# installing executables and libraries
install(TARGETS uu.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)

    message(STATUS "    Add example batch test fixed variables") 

    # Test run in working directory AFTER install of bb executable
    add_test(NAME ExampleAdvancedBatchFixedVariables
	COMMAND ./runFixed.sh 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
