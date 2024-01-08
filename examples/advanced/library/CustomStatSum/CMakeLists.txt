add_executable(customStatSum.exe customStatSum.cpp )

target_include_directories(customStatSum.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(customStatSum.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(customStatSum.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(customStatSum.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS customStatSum.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add test for this example
#if(BUILD_TESTS MATCHES ON)
#    message(STATUS "    Add example test for stop on stat sum criterion") 
#    # Can run this test after install
#    if (WIN32)
#        add_test(NAME ExampleAdvancedcustomStatSum
#	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./customStatSum.exe
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
#    else()
#        add_test(NAME ExampleAdvancedcustomStatSum
#	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./customStatSum.exe
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
#    endif()
#endif()
