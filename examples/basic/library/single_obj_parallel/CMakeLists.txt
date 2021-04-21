add_executable(single_obj_parallel.exe basic_lib.cpp )

target_include_directories(single_obj_parallel.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(single_obj_parallel.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")


if(OpenMP_CXX_FOUND)
    target_link_libraries(single_obj_parallel.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(single_obj_parallel.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS single_obj_parallel.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example library test single_obj_parallel") 

# Can run this test after install
add_test(NAME SingleObjParallelBasicLib
	COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./single_obj_parallel.exe
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
