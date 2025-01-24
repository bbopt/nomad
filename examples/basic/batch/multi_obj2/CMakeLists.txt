set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bbMO2.exe bb.cpp )
set_target_properties(bbMO2.exe PROPERTIES SUFFIX "")


# installing executables and libraries
install(TARGETS bbMO2.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# No test for this example