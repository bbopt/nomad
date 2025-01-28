set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb_disco.exe bb_disco.cpp)
set_target_properties(bb_disco.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb_disco.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# No test for this example
