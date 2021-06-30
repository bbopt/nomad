set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb3 bb.cpp )

target_link_libraries(bb3 PUBLIC OpenMP::OpenMP_CXX)

# installing executables and libraries
install(TARGETS bb3
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch test single_obj_parallel") 

# Test run in working directory AFTER install of bb executable
add_test(NAME ExampleBasicBatchSingleObjParallel
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
