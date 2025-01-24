add_executable(customOpport.exe customOpport.cpp )

target_include_directories(customOpport.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(customOpport.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(customOpport.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(customOpport.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS customOpport.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# No test for this example