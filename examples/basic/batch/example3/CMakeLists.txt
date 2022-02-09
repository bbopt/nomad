set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb3.exe bb3.cpp )
set_target_properties(bb3.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb3.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch #3") 

# Test run in working directory AFTER install of bb1.exe executable
add_test(NAME Example3BasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
