add_executable(rosenbrock.exe rosenbrock.cpp )

target_include_directories(rosenbrock.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(rosenbrock.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(rosenbrock.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(rosenbrock.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS rosenbrock.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example test for PSDMads") 
    # Can run this test after install

    if (WIN32)
        add_test(NAME ExampleAdvancedPSDMads
	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./rosenbrock.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    else()
        add_test(NAME ExampleAdvancedPSDMads
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./rosenbrock.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    endif()
endif()
