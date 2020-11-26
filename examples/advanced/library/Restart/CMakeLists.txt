add_executable(restart.exe restart.cpp )

target_include_directories(restart.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(restart.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")


if(OpenMP_CXX_FOUND)
    target_link_libraries(restart.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(restart.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS restart.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )
