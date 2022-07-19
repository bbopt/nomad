include(UseSWIG)
include(UseJava)

set(CMAKE_SWIG_FLAGS -package jNomad)
set_property(SOURCE jNomad.i PROPERTY CPLUSPLUS ON)
swig_add_library(jNomad TYPE SHARED LANGUAGE java SOURCES jNomad.i)
target_include_directories(jNomad PRIVATE ${JNI_INCLUDE_DIRS})

#
# Build
#
swig_link_libraries(jNomad PUBLIC nomadAlgos nomadUtils nomadEval)
set_target_properties(jNomad PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

#
# Jar
#
set(CMAKE_JNI_TARGET TRUE)
set(SWIG_JAVA_FILES
${CMAKE_CURRENT_BINARY_DIR}/AllParameters.java
${CMAKE_CURRENT_BINARY_DIR}/ArrayOfDouble.java
${CMAKE_CURRENT_BINARY_DIR}/Double.java
${CMAKE_CURRENT_BINARY_DIR}/EvalPoint.java
${CMAKE_CURRENT_BINARY_DIR}/EvalParameters.java
${CMAKE_CURRENT_BINARY_DIR}/EvalType.java
${CMAKE_CURRENT_BINARY_DIR}/Evaluator.java
${CMAKE_CURRENT_BINARY_DIR}/EvalXDefined.java
${CMAKE_CURRENT_BINARY_DIR}/jNomad.java
${CMAKE_CURRENT_BINARY_DIR}/jNomadConstants.java
${CMAKE_CURRENT_BINARY_DIR}/jNomadJNI.java
${CMAKE_CURRENT_BINARY_DIR}/MainStep.java
${CMAKE_CURRENT_BINARY_DIR}/Point.java
${CMAKE_CURRENT_BINARY_DIR}/Step.java
${CMAKE_CURRENT_BINARY_DIR}/SWIGTYPE_p_${NOMAD_SRC_TEMPLATE}__EvalParameters.java
${CMAKE_CURRENT_BINARY_DIR}/SWIGTYPE_p_std__vectorT_bool_t.java
${CMAKE_CURRENT_BINARY_DIR}/SWIGTYPE_p_std__vectorT_double_t.java
${CMAKE_CURRENT_BINARY_DIR}/SWIGTYPE_p_std__vectorT_std__shared_ptrT_${NOMAD_SRC_TEMPLATE}__EvalPoint_t_t.java
)
add_jar(jNomad_jar SOURCES ${SWIG_JAVA_FILES} OUTPUT_NAME "jNomad")

set_target_properties(jNomad_jar PROPERTIES CLASSDIR "${CMAKE_CURRENT_BINARY_DIR}/")

add_dependencies(jNomad_jar jNomad)

#
# Install binary, java, class and jar
#
install(TARGETS jNomad
       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
       PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/)

install_jar(jNomad_jar DESTINATION ${CMAKE_INSTALL_LIBDIR})

# This one is to run the example in the current source directory
install_jar(jNomad_jar DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS jNomad DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
