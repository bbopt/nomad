#
# Source and header
#
set(INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/NomadStdCInterface.cpp)

set(INTERFACE_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/NomadStdCInterface.h)


# 
# Build 
#
add_library (nomadCInterface SHARED ${INTERFACE_SOURCES}) 
target_include_directories(nomadCInterface PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/interfaces>
    $<INSTALL_INTERFACE:${CMAKE_BINARY_DIR}/interfaces>
)
target_link_libraries(nomadCInterface PUBLIC nomadAlgos nomadUtils nomadEval)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nomadCInterface PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties(nomadCInterface PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# 
# Install binary and header
#
install(TARGETS nomadCInterface
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
    )

install(FILES ${INTERFACE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/interfaces
    )
