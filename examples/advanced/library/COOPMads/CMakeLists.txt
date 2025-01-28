# Example only for OPENMP build

add_executable(rosenbrockCoop.exe rosenbrockCoop.cpp )

target_include_directories(rosenbrockCoop.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(rosenbrockCoop.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")

target_link_libraries(rosenbrockCoop.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)

# installing executables and libraries
install(TARGETS rosenbrockCoop.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add a test for this example
message(STATUS "    Add example test for COOPMads") 

if (WIN32)
    add_test(NAME ExampleAdvancedCOOPMads
	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./rosenbrockCoop.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
else()
    add_test(NAME ExampleAdvancedCOOPMads
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./rosenbrockCoop.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
endif()

