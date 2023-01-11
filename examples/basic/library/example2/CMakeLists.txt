add_executable(example2_lib.exe example2_lib.cpp )

target_include_directories(example2_lib.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(example2_lib.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(example2_lib.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared OpenMP::OpenMP_CXX)
else()
    target_link_libraries(example2_lib.exe PUBLIC nomadAlgosShared nomadUtilsShared nomadEvalShared)
endif()

# installing executables and libraries
install(TARGETS example2_lib.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example library #2") 

# Can run this test after install
if (WIN32)
# Todo Find a way to detect error
#    add_test(NAME Example2BasicLib
#	    COMMAND ./example2_lib.exe # Let ctest detect error. Not sure it can!
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
else()
    add_test(NAME Example2BasicLib
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example2_lib.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
