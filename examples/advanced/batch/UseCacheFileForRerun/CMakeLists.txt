set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb_cache.exe bb.cpp )
set_target_properties(bb_cache.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb_cache.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
# message(STATUS "    Add example batch for using cache file for rerun") 

# Test run in working directory AFTER install of bb_cache.exe executable
# add_test(NAME ExampleAdvancedBatchCacheFileForRerun
#	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
#	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
