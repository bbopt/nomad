find_package(MPI QUIET)

if(MPI_FOUND)

   include_directories(${MPI_INCLUDE_PATH})

   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")
   set(CMAKE_EXECUTABLE_SUFFIX .exe)

   add_executable(bb4.exe bb.cpp )
   set_target_properties(bb4.exe PROPERTIES SUFFIX "")

   target_link_libraries(bb4.exe PUBLIC MPI::MPI_CXX)

   # installing executables and libraries
   install(TARGETS bb4.exe
       RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

   message(STATUS "    Add example batch single_obj_MPIparallel (MPI found).") 

   # Test run in working directory AFTER install of bb executable
   add_test(NAME ExampleBasicBatchSingleObjMPIParallel
            COMMAND ${CMAKE_INSTALL_PREFIX}/bin/nomad param.txt 
	    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
endif()
