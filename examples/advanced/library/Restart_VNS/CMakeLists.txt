add_executable(restart_vns.exe restart_vns.cpp )

target_include_directories(restart_vns.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(restart_vns.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(restart_vns.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared OpenMP::OpenMP_CXX)
else()
    target_link_libraries(restart_vns.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared)
endif()

# installing executables and libraries
install(TARGETS restart_vns.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add a test for this example
#if(BUILD_TESTS MATCHES ON)
#    message(STATUS "    Add example test for restart VNS") 
#    # Can run this test after install
#    if (WIN32)
#        add_test(NAME ExampleAdvancedRestart
#	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./restart_vns.exe
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
#    else()
#        add_test(NAME ExampleAdvancedRestart
#	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./restart_vns.exe
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
#    endif()
#endif()
