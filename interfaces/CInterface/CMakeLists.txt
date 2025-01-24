#
# Source and header
#
set(
  INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/NomadStdCInterface.cpp
)

set(
  INTERFACE_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/NomadStdCInterface.h
)


#
# Build
#
add_library(
  nomadCInterface
  SHARED ${INTERFACE_SOURCES}
) 

target_include_directories(
  nomadCInterface
  PUBLIC
    $<BUILD_INTERFACE: ${CMAKE_BINARY_DIR}/interfaces/CInterface>
    $<INSTALL_INTERFACE: ${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
  nomadCInterface
  PUBLIC
    nomadAlgos
    nomadUtils
    nomadEval
)

if (WIN32)
  target_compile_definitions(
    nomadCInterface 
    PRIVATE -DNOMAD_INTERFACE_C_DLL
  )
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(
    nomadCInterface 
    PUBLIC OpenMP::OpenMP_CXX
  )
endif()

if (NOT WIN32)
  target_link_libraries(
    nomadCInterface
    PRIVATE m
  )
endif()

set_target_properties(
  nomadCInterface
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    VERSION "${NOMAD_VERSION}"
    SOVERSION "${NOMAD_VERSION_MAJOR}.${NOMAD_VERSION_MINOR}"
    NO_SONAME ${NO_SONAME}
)


#
# Install binary and header
#
install(
  TARGETS nomadCInterface
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/
)

install(
  FILES ${INTERFACE_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/interfaces/CInterface
)

