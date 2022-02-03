# 
# Build 
#


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   matlab_add_mex(NAME nomadOpt SHARED SRC Source/nomadmex.cpp LINK_TO nomadAlgos nomadUtils nomadEval ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft/libut.lib)
else()
   matlab_add_mex(NAME nomadOpt SHARED SRC Source/nomadmex.cpp LINK_TO nomadAlgos nomadUtils nomadEval)
endif()


set_target_properties(nomadOpt PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR};${Matlab_ROOT_DIR}/sys/os/glnxa64/;${Matlab_ROOT_DIR}/sys/os/maci64/" )

# 
# Install binary and header
#
install(TARGETS nomadOpt
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/
    )

