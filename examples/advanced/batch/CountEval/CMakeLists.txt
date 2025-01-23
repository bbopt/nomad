set(CMAKE_EXECUTABLE_SUFFIX .exe)
add_executable(bb_cnt_eval.exe bb_cnt_eval.cpp )
set_target_properties(bb_cnt_eval.exe PROPERTIES SUFFIX "")

# installing executables and libraries
install(TARGETS bb_cnt_eval.exe
    RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

# Add a test for this example
message(STATUS "    Add example advanced batch count eval") 

# Test run in working directory AFTER install of bb_cnt_eval.exe executable
add_test(NAME CountEvalAdvancedBatch
	COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
