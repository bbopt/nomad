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


set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(single_bbob-constrained.exe single_bbob-constrained.c )

set_target_properties(single_bbob-constrained.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")
target_link_libraries(single_bbob-constrained.exe PUBLIC coco)

if(MATH_LIBRARY)
    target_link_libraries(single_bbob-constrained.exe PUBLIC ${MATH_LIBRARY})
endif()

# installing executables and libraries
install(TARGETS single_bbob-constrained.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example batch single obj bbob-constrained") 

# Test run in working directory AFTER install of single_bbob-constrained.exe executable
add_test(NAME ExampleSingleBbobConsBasicBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
