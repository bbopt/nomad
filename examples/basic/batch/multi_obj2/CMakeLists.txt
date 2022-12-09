set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bbMO2.exe bb.cpp )
set_target_properties(bbMO2.exe PROPERTIES SUFFIX "")


# installing executables and libraries
install(TARGETS bbMO2.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch multi_obj 2") 

# Test run in working directory AFTER install of bb.exe executable
add_test(NAME ExampleMultiObj2BasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
