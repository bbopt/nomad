add_executable(NMOpt.exe NMOpt.cpp )

target_include_directories(NMOpt.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(NMOpt.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(NMOpt.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared OpenMP::OpenMP_CXX)
else()
    target_link_libraries(NMOpt.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared)
endif()

# installing executables and libraries
install(TARGETS NMOpt.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
if(BUILD_TESTS MATCHES ON)
    message(STATUS "    Add example test for NM only") 
    # Can run this test after install

    if (WIN32)
# Todo Find a way to detect error
#        add_test(NAME ExampleAdvancedNMOnly
#	    COMMAND ./NMOpt.exe # Let ctest detect error. Not sure it can!
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
    else()
        add_test(NAME ExampleAdvancedNMOnly
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./NMOpt.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
    endif()
endif()
