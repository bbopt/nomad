add_executable(example.exe example.cpp )

target_include_directories(example.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(example.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")


if(OpenMP_CXX_FOUND)
    target_link_libraries(example.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(example.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS example.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example library for custom h computation") 

    # Can run this test after install
    add_test(NAME ExampleAdvancedCustomHCompute
	COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example.exe
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
    )
endif()
