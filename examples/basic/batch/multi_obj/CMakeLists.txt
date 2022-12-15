set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bbMO.exe bb.cpp )
set_target_properties(bbMO.exe PROPERTIES SUFFIX "")


# installing executables and libraries
install(TARGETS bbMO.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch multi_obj") 

# Test run in working directory AFTER install of bb.exe executable
add_test(NAME ExampleMultiObjBasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
