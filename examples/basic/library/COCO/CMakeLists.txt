if((NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/coco.c) OR (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/coco.h))
    return()
endif()

find_library(MATH_LIBRARY m)

include(CheckCompilerFlag)
check_compiler_flag(C "-pedantic -Wall -Wextra -Wstrict-prototypes -Wshadow -Wno-sign-compare -Wconversion" CC_HAS_WALL_ETC)

add_library(coco STATIC coco.c coco.h)
target_include_directories(coco PUBLIC .)
if(MATH_LIBRARY)
    target_link_libraries(coco PUBLIC ${MATH_LIBRARY})
endif()
# Add warning flags
if (MSVC)
  target_compile_options(coco PRIVATE "/W3")
elseif (CC_HAS_WALL_ETC)
  target_compile_options(coco PRIVATE -pedantic -Wall -Wextra -Wstrict-prototypes -Wshadow -Wno-sign-compare -Wconversion)
endif()

add_executable(coco_Nomad_experiment.exe coco_Nomad_experiment.cpp)

set_target_properties(coco_Nomad_experiment.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")

if(OpenMP_CXX_FOUND)
    target_link_libraries(coco_Nomad_experiment.exe PUBLIC coco nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(coco_Nomad_experiment.exe PUBLIC coco nomadAlgos nomadUtils nomadEval)    
endif()

if(MATH_LIBRARY)
    target_link_libraries(coco_Nomad_experiment.exe PUBLIC ${MATH_LIBRARY})
endif()




# 
# Install binary and header
#
install(TARGETS coco_Nomad_experiment.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
    )
