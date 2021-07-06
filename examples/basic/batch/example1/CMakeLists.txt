set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb1.exe bb.cpp )
set_target_properties(bb1.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb1.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch test 1") 

# Test run in working directory AFTER install of bb1.exe executable
add_test(NAME Example1BasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
