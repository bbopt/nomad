add_executable(customSurrogateOrdering.exe customSurrogateOrdering.cpp )

target_include_directories(customSurrogateOrdering.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(customSurrogateOrdering.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(customSurrogateOrdering.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(customSurrogateOrdering.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS customSurrogateOrdering.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# No test for this example