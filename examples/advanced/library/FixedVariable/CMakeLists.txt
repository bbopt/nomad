add_executable(fixedVariable.exe fixedVariable.cpp )
add_executable(ufl.exe u.cpp )

target_include_directories(fixedVariable.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(fixedVariable.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")
set_target_properties(ufl.exe PROPERTIES SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(fixedVariable.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(fixedVariable.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS fixedVariable.exe ufl.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example test for fixed variables") 
    # Can run this test after install

    if (WIN32)
        add_test(NAME ExampleAdvancedFixedVariables
	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./fixedVariable.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    else()
        add_test(NAME ExampleAdvancedFixedVariables
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./fixedVariable.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    endif()
endif()
