add_executable(suggest.exe suggest.cpp )
add_executable(observe.exe observe.cpp )
add_executable(bbr.exe bbr.cpp )

add_library (utils OBJECT utils.cpp)

set_target_properties(suggest.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set_target_properties(observe.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

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

# Add a test for this example
if(BUILD_TESTS MATCHES ON)

    message(STATUS "    Add example batch test suggest and observe") 

    # Test run in working directory AFTER install of bb executable
    add_test(NAME ExampleAdvancedBatchSuggestAndObserve
	COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./loopSuggestAndObserve.sh
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
