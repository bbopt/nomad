set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(uu.exe uu.cpp )
set_target_properties(uu.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS uu.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# No test for this example