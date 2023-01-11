add_executable(restart.exe restart.cpp )

target_include_directories(restart.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(restart.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(restart.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared OpenMP::OpenMP_CXX)
else()
    target_link_libraries(restart.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared)
endif()

# installing executables and libraries
install(TARGETS restart.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example test for restart") 
    # Can run this test after install
    if (WIN32)
# Todo Find a way to detect error
#        add_test(NAME ExampleAdvancedRestart
#	    COMMAND ./restart.exe # Let ctest detect error. Not sure it can!
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
    else()
        add_test(NAME ExampleAdvancedRestart
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./restart.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    endif()
endif()
