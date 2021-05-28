add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/library/example1)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/library/example2)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/library/example3)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/library/single_obj_parallel)


# Issue #491
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/exampleForCustomHComputation)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/FixedVariable)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/NMonly)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/Restart)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/exampleSuggestAndObserve)
if(OpenMP_CXX_FOUND)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/PSDMads)
endif()

if (BUILD_INTERFACES MATCHES ON)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/c_api/example1)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/library/c_api/example2)
endif()

# For running the tests on the batch examples we need to build the blackbox
# AND install
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/batch/example1)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/batch/example2)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/batch/single_obj)
if(OpenMP_CXX_FOUND)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/batch/single_obj_parallel)
endif()
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/batch/LHonly)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/batch/FixedVariable)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/advanced/batch/SuggestAndObserve)

# The script for running library examples is created in a temp directory
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/tmp/runExampleTest.sh
    "#!/bin/bash

    executable=$1
    param_file=\"empty\" 
    if [ \$# == 2 ]
    then
        param_file =$2  
    fi
    echo \"Running Example test $executable\"
    $executable $param_file | tee runTestLog.txt
    if [ $? != 0 ]
    then
        exit 1
    else
        test_passed=`grep termination runTestLog.txt`
        test_failed=`grep -e exception -e Exception -e Error -e abort -e terminate runTestLog.txt`
        if [ \"$test_passed\" != \"\" ] \&\& [ \"$test_failed\" == \"\" ]
        then
            exit 0
        else
            exit 1
        fi
    fi
    \n")

# This manipulation is for changing permission of the script
file (COPY ${CMAKE_CURRENT_BINARY_DIR}/tmp/runExampleTest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)

# Remove the temp directory
file (REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/tmp)

