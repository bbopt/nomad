add_executable(fixedVariable.exe fixedVariable.cpp )

target_include_directories(fixedVariable.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)

set_target_properties(fixedVariable.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")


if(OpenMP_CXX_FOUND)
    target_link_libraries(fixedVariable.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(fixedVariable.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS fixedVariable.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )
