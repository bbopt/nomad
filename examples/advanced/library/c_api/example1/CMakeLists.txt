add_executable(example1_c_api.exe example1_c_api.c )

target_include_directories(
  example1_c_api.exe 
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/interfaces/CInterface
)

set_target_properties(
  example1_c_api.exe 
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" 
    SUFFIX ""
)


if(OpenMP_CXX_FOUND)
   target_link_libraries(
     example1_c_api.exe 
     PUBLIC 
       nomadCInterface 
       OpenMP::OpenMP_CXX
   )
else()
   target_link_libraries(
     example1_c_api.exe
     PUBLIC nomadCInterface
   )
endif()

if (NOT WIN32)
  target_link_libraries(
    example1_c_api.exe 
    PUBLIC m
  )
endif()

# installing executables and libraries
install(TARGETS example1_c_api.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example test for c-api 1") 
    
# Can run this test after install

if (WIN32)
#  MAYBE ENABLE the test when script is adapted to windows
#        add_test(NAME ExampleAdvancedC-API1
#	    COMMAND bash.exe ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example1_c_api.exe
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
#        )
else()
    add_test(NAME ExampleAdvancedC-API1
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example1_c_api.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        )
endif()
