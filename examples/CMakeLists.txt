# Basic cpp examples in library mode
add_executable(basic_lib_example1 basic/library/example1/example1_lib.cpp)
target_link_libraries(basic_lib_example1 PUBLIC nomadAlgos nomadUtils nomadEval)
target_include_directories(basic_lib_example1 PUBLIC #${CMAKE_SOURCE_DIR}/src)
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_executable(basic_lib_example2 basic/library/example2/example2_lib.cpp)
target_link_libraries(basic_lib_example2 PUBLIC nomadAlgos nomadUtils nomadEval)
target_include_directories(basic_lib_example2 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_executable(basic_lib_example3 basic/library/example3/example3_lib.cpp)
target_link_libraries(basic_lib_example3 PUBLIC nomadAlgos nomadUtils nomadEval)
target_include_directories(basic_lib_example3 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Advanced cpp examples in library mode
add_executable(exampleForCustomHComputation advanced/library/exampleForCustomHComputation/example.cpp)
target_link_libraries(exampleForCustomHComputation PUBLIC nomadAlgos nomadUtils nomadEval)
target_include_directories(exampleForCustomHComputation PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
