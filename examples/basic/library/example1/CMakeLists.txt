add_executable(example1_lib.exe example1_lib.cpp )

target_include_directories(example1_lib.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(example1_lib.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(example1_lib.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(example1_lib.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS example1_lib.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example library #1") 

# Can run this test after install
if (WIN32)
    add_test(NAME Example1BasicLib
	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example1_lib.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
else()
    add_test(NAME Example1BasicLib
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example1_lib.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
