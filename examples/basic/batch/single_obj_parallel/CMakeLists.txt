set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb_s.exe bb.cpp )
set_target_properties(bb_s.exe PROPERTIES SUFFIX "")

target_link_libraries(bb_s.exe PUBLIC OpenMP::OpenMP_CXX)

# installing executables and libraries
install(TARGETS bb_s.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch single_obj_parallel") 

# Test run in working directory AFTER install of bb executable
add_test(NAME ExampleBasicBatchSingleObjParallel
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
