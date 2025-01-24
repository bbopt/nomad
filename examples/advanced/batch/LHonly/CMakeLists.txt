set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(u.exe u.cpp )
set_target_properties(u.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS u.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# No test for this example
