add_executable(multi_obj_lib.exe multi_obj_lib.cpp)

target_include_directories(multi_obj_lib.exe PRIVATE ${CMAKE_SOURCE_DIR}/src)

set_target_properties(multi_obj_lib.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")

if(OpenMP_CXX_FOUND)
    target_link_libraries(multi_obj_lib.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(multi_obj_lib.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS multi_obj_lib.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add multi obj library example")

# Can run this test after install
if (WIN32)
    add_test(NAME MultiObjBasicLib
	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./multi_obj_lib.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
else()
    add_test(NAME MultiObjBasicLib
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./multi_obj_lib.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
