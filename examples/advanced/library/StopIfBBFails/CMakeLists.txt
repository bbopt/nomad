add_executable(stopIfBBFails.exe stopIfBBFails.cpp )

target_include_directories(stopIfBBFails.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(stopIfBBFails.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(stopIfBBFails.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(stopIfBBFails.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS stopIfBBFails.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# Add test for this example
#if(BUILD_TESTS MATCHES ON)
#    message(STATUS "    Add example test for handling bb fail") 
#    # Can run this test after install
#    if (WIN32)
#        add_test(NAME ExampleAdvancedstopIfBBFails
#	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./stopIfBBFails.exe
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
#    else()
#        add_test(NAME ExampleAdvancedstopIfBBFails
#	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./stopIfBBFails.exe
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
#    endif()
#endif()
