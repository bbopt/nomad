add_executable(ss_bb.exe bb.cpp )
add_executable(ss_sgte.exe sgte.cpp )

# installing executables and libraries
install(TARGETS ss_bb.exe ss_sgte.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch test surrogate sort") 

# Test run in working directory AFTER install of bb executable
add_test(NAME ExampleSurrogateSortBasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
