set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb_redirect.exe bb_redirect.cpp )
set_target_properties(bb_redirect.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb_redirect.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example advanced batch BB redirection") 

# Test run in working directory AFTER install of bb_redirect.exe executable
add_test(NAME ExampleAdvancedBatchRedirect
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
