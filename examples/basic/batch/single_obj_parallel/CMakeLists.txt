add_executable(bb3.exe bb.cpp )

# installing executables and libraries
install(TARGETS bb3.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example batch test single_obj_parallel") 

    # Test run in working directory AFTER install of bb executable
    add_test(NAME ExampleBasicBatchSingleObjParallel
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    )
endif()
