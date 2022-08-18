add_executable(example2_c_api.exe example2_c_api.c )

target_include_directories(example2_c_api.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/interfaces/CInterface)

set_target_properties(example2_c_api.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(example2_c_api.exe PUBLIC nomadAlgos nomadUtils nomadEval nomadCInterface OpenMP::OpenMP_CXX)
else()
    target_link_libraries(example2_c_api.exe PUBLIC nomadAlgos nomadUtils nomadEval nomadCInterface)
endif()

# installing executables and libraries
install(TARGETS example2_c_api.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example test for c-api 2") 
    # Can run this test after install

    if (WIN32)
        # Todo Find a way to detect error
        #add_test(NAME ExampleAdvancedC-API2
	#    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example2_c_api.exe
	#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        #)
    else()
        add_test(NAME ExampleAdvancedC-API2
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example2_c_api.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    endif()
endif()
