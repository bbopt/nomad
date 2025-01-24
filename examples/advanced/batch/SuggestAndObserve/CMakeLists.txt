set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(suggest.exe suggest.cpp )
add_executable(observe.exe observe.cpp )
add_executable(bbr.exe bbr.cpp )

add_library (utils OBJECT utils.cpp)

set_target_properties(suggest.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")
set_target_properties(observe.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")
set_target_properties(bbr.exe PROPERTIES SUFFIX "")

if(OpenMP_CXX_FOUND)
    target_link_libraries(utils PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
    target_link_libraries(suggest.exe PUBLIC utils nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
    target_link_libraries(observe.exe PUBLIC utils nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(utils PUBLIC nomadAlgos nomadUtils nomadEval)
    target_link_libraries(suggest.exe PUBLIC utils nomadAlgos nomadUtils nomadEval)
    target_link_libraries(observe.exe PUBLIC utils nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS suggest.exe observe.exe bbr.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# No test for this example
