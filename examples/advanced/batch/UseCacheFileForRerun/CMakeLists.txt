set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb_cache.exe bb.cpp )
set_target_properties(bb_cache.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb_cache.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# No test for this example
