set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb_disco.exe bb_disco.cpp)
set_target_properties(bb_disco.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb_disco.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example batch DiscoMads") 
    # Test run in working directory AFTER install of bb executable
    add_test(NAME ExampleAdvancedBatchDiscoMadsDiscont
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad paramEscapeDiscont.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
