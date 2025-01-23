add_executable(example3_c_api.exe example3_c_api.c)

target_include_directories(
  example3_c_api.exe
    PRIVATE
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_SOURCE_DIR}/interfaces/CInterface
)

set_target_properties(
  example3_c_api.exe
  PROPERTIES
    INSTALL_RPATH
      "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    SUFFIX ""
)


if(OpenMP_CXX_FOUND)
  target_link_libraries(
    example3_c_api.exe
    PUBLIC
      nomadCInterface
      OpenMP::OpenMP_CXX
  )
else()
  target_link_libraries(
    example3_c_api.exe
    PUBLIC nomadCInterface
  )
endif()

if (NOT WIN32)
  target_link_libraries(
    example3_c_api.exe
    PUBLIC m
  )
endif()

# installing executables and libraries
install(TARGETS example3_c_api.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )


# No test for this example
