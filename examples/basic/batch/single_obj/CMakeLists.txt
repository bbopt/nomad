set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb bb.cpp )

# installing executables and libraries
install(TARGETS bb
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch test single_obj") 

# Test run in working directory AFTER install of bb executable
add_test(NAME ExampleSingleObjBasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
