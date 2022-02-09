set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb2.exe bb.cpp )
set_target_properties(bb2.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb2.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch #2") 

# Test run in working directory AFTER install of bb executable
add_test(NAME Example2BasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
