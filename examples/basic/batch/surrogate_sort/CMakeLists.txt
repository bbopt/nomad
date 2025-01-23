add_executable(ss_bb.exe bb.cpp )
add_executable(ss_sgte.exe sgte.cpp )
set_target_properties(ss_bb.exe PROPERTIES SUFFIX "")
set_target_properties(ss_sgte.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS ss_bb.exe ss_sgte.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# No test for this example
