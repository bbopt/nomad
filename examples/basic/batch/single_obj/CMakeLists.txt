add_executable(bb.exe bb.cpp )

# installing executables and libraries
install(TARGETS bb.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example batch test single_obj") 

    # Test run in working directory AFTER install of bb executable
    add_test(NAME ExampleSingleObjBasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    )
endif()
