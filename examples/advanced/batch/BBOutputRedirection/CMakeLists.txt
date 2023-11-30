set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb5.exe bb5.cpp )
set_target_properties(bb5.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb5.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch #5") 

# Test run in working directory AFTER install of bb5.exe executable
add_test(NAME Example5BasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
