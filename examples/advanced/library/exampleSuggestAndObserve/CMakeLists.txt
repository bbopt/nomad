add_executable(suggestAndObserve.exe suggestAndObserve.cpp )
add_executable(loopSuggestAndObserve.exe loopSuggestAndObserve.cpp )

target_include_directories(suggestAndObserve.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)
target_include_directories(loopSuggestAndObserve.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src)


set_target_properties(suggestAndObserve.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")
set_target_properties(loopSuggestAndObserve.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" SUFFIX "")


if(OpenMP_CXX_FOUND)
    target_link_libraries(suggestAndObserve.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
    target_link_libraries(loopSuggestAndObserve.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(suggestAndObserve.exe PUBLIC nomadAlgos nomadUtils nomadEval)
    target_link_libraries(loopSuggestAndObserve.exe PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS suggestAndObserve.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )
install(TARGETS loopSuggestAndObserve.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example library tests using Suggest And Observe") 

# Can run this test after install
if (WIN32)
# Todo Find a way to detect error
#    add_test(NAME ExampleAdvancedSuggestAndObserve
#	    COMMAND ./suggestAndObserve.exe # Let ctest detect error. Not sure it can!
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
#    add_test(NAME ExampleAdvancedLoopSuggestAndObserve
#	    COMMAND ./loopSuggestAndObserve.exe # Let ctest detect error. Not sure it can!
#	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
else()
    add_test(NAME ExampleAdvancedSuggestAndObserve
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./suggestAndObserve.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_test(NAME ExampleAdvancedLoopSuggestAndObserve
	    COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./loopSuggestAndObserve.exe
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()

