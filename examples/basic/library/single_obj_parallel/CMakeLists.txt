add_executable(single_obj_parallel.exe basic_lib.cpp )

target_include_directories(single_obj_parallel.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(single_obj_parallel.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(single_obj_parallel.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(single_obj_parallel.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS single_obj_parallel.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example library single_obj_parallel") 

# Can run this test after install
if (WIN32)
    add_test(NAME SingleObjParallelBasicLib
	    COMMAND ./single_obj_parallel.exe # Let ctest detect error. Not sure it can!
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
else()
    add_test(NAME SingleObjParallelBasicLib
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./single_obj_parallel.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
