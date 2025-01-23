
add_executable(nestedOptim.exe nestedOptim.cpp)

target_include_directories(nestedOptim.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(nestedOptim.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")

target_link_libraries(nestedOptim.exe PUBLIC nomadAlgos nomadUtils nomadEval)

# installing executables and libraries
install(TARGETS nestedOptim.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add a test for this example
message(STATUS "    Add example for nested optimization") 

if (WIN32)
     add_test(NAME ExampleAdvancedNestedMads
	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./nestedOptim.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
else()
     add_test(NAME ExampleAdvancedNestedMads
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./nestedOptim.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
endif()
endif()
