set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(u u.cpp )

# installing executables and libraries
install(TARGETS u
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example batch test LHOnly") 

    # Test run in working directory AFTER install of bb executable
    add_test(NAME ExampleAdvancedBatchLHOnly
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    )
endif()
