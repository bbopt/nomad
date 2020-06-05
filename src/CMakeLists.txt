# Main modules

# Algos
set(ALGOS_HEADERS
    Algos/AlgoStopReasons.hpp
    Algos/Algorithm.hpp
    Algos/AllStopReasons.hpp
    Algos/CacheInterface.hpp
    Algos/EvcInterface.hpp
    Algos/Initialization.hpp
    Algos/Iteration.hpp
    Algos/IterationUtils.hpp
    Algos/MainStep.hpp
    Algos/MegaIteration.hpp
    Algos/MeshBase.hpp
    Algos/Projection.hpp
    Algos/Step.hpp
    Algos/Subproblem.hpp
    Algos/Termination.hpp)

set(ALGOS_SOURCES
    Algos/AlgoStopReasons.cpp
    Algos/Algorithm.cpp
    Algos/AllStopReasons.cpp
    Algos/CacheInterface.cpp
    Algos/EvcInterface.cpp
    Algos/Initialization.cpp
    Algos/Iteration.cpp
    Algos/IterationUtils.cpp
    Algos/MainStep.cpp
    Algos/MegaIteration.cpp
    Algos/MeshBase.cpp
    Algos/Projection.cpp
    Algos/Step.cpp
    Algos/Subproblem.cpp
    Algos/Termination.cpp)

# Attribute
set(ATTRIBUTE_HEADERS
    Attribute/cacheAttributesDefinition.hpp
    Attribute/displayAttributesDefinition.hpp
    Attribute/evalAttributesDefinition.hpp
    Attribute/evaluatorControlAttributesDefinition.hpp
    Attribute/pbAttributesDefinition.hpp
    Attribute/runAttributesDefinition.hpp)

set(ATTRIBUTE_SOURCES
    Attribute/WriteAttributeDefinitionFile.cpp)

# Cache
set(CACHE_HEADERS
    Cache/CacheBase.hpp
    Cache/CacheSet.hpp)

set(CACHE_SOURCES
    Cache/CacheBase.cpp
    Cache/CacheSet.cpp)

# Eval
set(EVAL_HEADERS
    Eval/BBInput.hpp
    Eval/BBOutput.hpp
    Eval/Barrier.hpp
    Eval/Eval.hpp
    Eval/EvalPoint.hpp
    Eval/EvalQueuePoint.hpp
    Eval/Evaluator.hpp
    Eval/EvaluatorControl.hpp)

set(EVAL_SOURCES
    Eval/BBInput.cpp
    Eval/BBOutput.cpp
    Eval/Barrier.cpp
    Eval/Eval.cpp
    Eval/EvalPoint.cpp
    Eval/EvalQueuePoint.cpp
    Eval/Evaluator.cpp
    Eval/EvaluatorControl.cpp)

# Interfaces
set(INTERFACE_SOURCES
    Interfaces/NomadStdCInterface.cpp)

set(INTERFACE_HEADERS
    Interfaces/NomadStdCInterface.h)

# Math
set(MATH_HEADERS
    Math/ArrayOfDouble.hpp
    Math/ArrayOfPoint.hpp
    Math/Direction.hpp
    Math/Double.hpp
    Math/LHS.hpp
    Math/MatrixUtils.hpp
    Math/Point.hpp
    Math/RNG.hpp)

set(MATH_SOURCES
    Math/ArrayOfDouble.cpp
    Math/ArrayOfPoint.cpp
    Math/Direction.cpp
    Math/Double.cpp
    Math/LHS.cpp
    Math/MatrixUtils.cpp
    Math/Point.cpp
    Math/RNG.cpp)

# Nomad
set(NOMAD_HEADERS
    Nomad/nomad.hpp)

set(NOMAD_SOURCES
    Nomad/nomad.cpp)

# Output
set(OUTPUT_HEADERS
    Output/OutputInfo.hpp
    Output/OutputQueue.hpp
    Output/StatsInfo.hpp)

set(OUTPUT_SOURCES
    Output/OutputInfo.cpp
    Output/OutputQueue.cpp
    Output/StatsInfo.cpp)

# Param
set(PARAM_HEADERS
    Param/AllParameters.hpp
    Param/Attribute.hpp
    Param/AttributeFactory.hpp
    Param/CacheParameters.hpp
    Param/DisplayParameters.hpp
    Param/EvalParameters.hpp
    Param/EvaluatorControlParameters.hpp
    Param/ParameterEntries.hpp
    Param/ParameterEntry.hpp
    Param/Parameters.hpp
    Param/ParametersNomad3.hpp
    Param/PbParameters.hpp
    Param/RunParameters.hpp
    Param/TypeAttribute.hpp)

set(PARAM_SOURCES
    Param/AllParameters.cpp
    Param/Attribute.cpp
    Param/AttributeFactory.cpp
    Param/CacheParameters.cpp
    Param/DisplayParameters.cpp
    Param/EvalParameters.cpp
    Param/EvaluatorControlParameters.cpp
    Param/ParameterEntries.cpp
    Param/ParameterEntry.cpp
    Param/Parameters.cpp
    Param/ParametersNomad3.cpp
    Param/PbParameters.cpp
    Param/RunParameters.cpp
    Param/TypeAttribute.cpp)

# Type
set(TYPE_HEADERS
    Type/BBInputType.hpp
    Type/BBOutputType.hpp
    Type/CallbackType.hpp
    Type/EvalType.hpp
    Type/LHSearchType.hpp
    Type/SgtelibModelFeasibilityType.hpp
    Type/SgtelibModelFormulationType.hpp)

set(TYPE_SOURCES
    Type/BBInputType.cpp
    Type/BBOutputType.cpp
    Type/CallbackType.cpp
    Type/EvalType.cpp
    Type/LHSearchType.cpp
    Type/SgtelibModelFeasibilityType.cpp
    Type/SgtelibModelFormulationType.cpp)

# Util
set(UTIL_HEADERS
    Util/ArrayOfString.hpp
    Util/Clock.hpp
    Util/Copyright.hpp
    Util/Exception.hpp
    Util/StopReason.hpp
    Util/Uncopyable.hpp
    Util/defines.hpp
    Util/fileutils.hpp
    Util/utils.hpp)

set(UTIL_SOURCES
    Util/ArrayOfString.cpp
    Util/Clock.cpp
    Util/Copyright.cpp
    Util/Exception.cpp
    Util/StopReason.cpp
    Util/Uncopyable.cpp
    Util/defines.cpp
    Util/fileutils.cpp
    Util/utils.cpp)

# Submodules

# LH
set(LH_HEADERS
    Algos/LatinHypercubeSampling/LH.hpp)

set(LH_SOURCES
    Algos/LatinHypercubeSampling/LH.cpp)

# Mads
set(MADS_HEADERS
    Algos/Mads/GMesh.hpp
    Algos/Mads/LHSearchMethod.hpp
    Algos/Mads/Mads.hpp
    Algos/Mads/MadsInitialization.hpp
    Algos/Mads/MadsIteration.hpp
    Algos/Mads/MadsMegaIteration.hpp
    Algos/Mads/MadsUpdate.hpp
    Algos/Mads/MegaSearchPoll.hpp
    Algos/Mads/NMSearchMethod.hpp
    Algos/Mads/Poll.hpp
    Algos/Mads/QuadSearchMethod.hpp
    Algos/Mads/Search.hpp
    Algos/Mads/SearchMethodAlgo.hpp
    Algos/Mads/SearchMethodBase.hpp
    Algos/Mads/SearchMethodSimple.hpp
    Algos/Mads/SgtelibSearchMethod.hpp
    Algos/Mads/SpeculativeSearchMethod.hpp
    Algos/Mads/UserSearchMethod.hpp)

set(MADS_SOURCES
    Algos/Mads/GMesh.cpp
    Algos/Mads/LHSearchMethod.cpp
    Algos/Mads/Mads.cpp
    Algos/Mads/MadsInitialization.cpp
    Algos/Mads/MadsIteration.cpp
    Algos/Mads/MadsMegaIteration.cpp
    Algos/Mads/MadsUpdate.cpp
    Algos/Mads/MegaSearchPoll.cpp
    Algos/Mads/NMSearchMethod.cpp
    Algos/Mads/Poll.cpp
    Algos/Mads/QuadSearchMethod.cpp
    Algos/Mads/Search.cpp
    Algos/Mads/SearchMethodAlgo.cpp
    Algos/Mads/SearchMethodBase.cpp
    Algos/Mads/SearchMethodSimple.cpp
    Algos/Mads/SgtelibSearchMethod.cpp
    Algos/Mads/SpeculativeSearchMethod.cpp
    Algos/Mads/UserSearchMethod.cpp)

# NM
set(NM_HEADERS
    Algos/NelderMead/NM.hpp
    Algos/NelderMead/NMAllReflective.hpp
    Algos/NelderMead/NMInitialization.hpp
    Algos/NelderMead/NMInitializeSimplex.hpp
    Algos/NelderMead/NMIteration.hpp
    Algos/NelderMead/NMIterationUtils.hpp
    Algos/NelderMead/NMMegaIteration.hpp
    Algos/NelderMead/NMReflective.hpp
    Algos/NelderMead/NMShrink.hpp
    Algos/NelderMead/NMSimplexEvalPoint.hpp
    Algos/NelderMead/NMUpdate.hpp)

set(NM_SOURCES
    Algos/NelderMead/NM.cpp
    Algos/NelderMead/NMAllReflective.cpp
    Algos/NelderMead/NMInitialization.cpp
    Algos/NelderMead/NMInitializeSimplex.cpp
    Algos/NelderMead/NMIteration.cpp
    Algos/NelderMead/NMIterationUtils.cpp
    Algos/NelderMead/NMMegaIteration.cpp
    Algos/NelderMead/NMReflective.cpp
    Algos/NelderMead/NMShrink.cpp
    Algos/NelderMead/NMSimplexEvalPoint.cpp
    Algos/NelderMead/NMUpdate.cpp)

# PhaseOne
set(PHASEONE_HEADERS
    Algos/PhaseOne/PhaseOne.hpp)

set(PHASEONE_SOURCES
    Algos/PhaseOne/PhaseOne.cpp)

# QuadModel
set(QUADMODEL_HEADERS
    Algos/QuadModel/QuadModelAlgo.hpp
    Algos/QuadModel/QuadModelEvaluator.hpp
    Algos/QuadModel/QuadModelInitialization.hpp
    Algos/QuadModel/QuadModelIteration.hpp
    Algos/QuadModel/QuadModelIterationUtils.hpp
    Algos/QuadModel/QuadModelMegaIteration.hpp
    Algos/QuadModel/QuadModelOptimize.hpp
    Algos/QuadModel/QuadModelSinglePass.hpp
    Algos/QuadModel/QuadModelUpdate.hpp)

set(QUADMODEL_SOURCES
    Algos/QuadModel/QuadModelAlgo.cpp
    Algos/QuadModel/QuadModelEvaluator.cpp
    Algos/QuadModel/QuadModelInitialization.cpp
    Algos/QuadModel/QuadModelIteration.cpp
    Algos/QuadModel/QuadModelIterationUtils.cpp
    Algos/QuadModel/QuadModelMegaIteration.cpp
    Algos/QuadModel/QuadModelOptimize.cpp
    Algos/QuadModel/QuadModelSinglePass.cpp
    Algos/QuadModel/QuadModelUpdate.cpp)

# SgtelibModel
set(SGTELIBMODEL_HEADERS
    Algos/SgtelibModel/SgtelibModel.hpp
    Algos/SgtelibModel/SgtelibModelEvaluator.hpp
    Algos/SgtelibModel/SgtelibModelFilterCache.hpp
    Algos/SgtelibModel/SgtelibModelInitialization.hpp
    Algos/SgtelibModel/SgtelibModelIteration.hpp
    Algos/SgtelibModel/SgtelibModelMegaIteration.hpp
    Algos/SgtelibModel/SgtelibModelOptimize.hpp
    Algos/SgtelibModel/SgtelibModelUpdate.hpp)

set(SGTELIBMODEL_SOURCES
    Algos/SgtelibModel/SgtelibModel.cpp
    Algos/SgtelibModel/SgtelibModelEvaluator.cpp
    Algos/SgtelibModel/SgtelibModelFilterCache.cpp
    Algos/SgtelibModel/SgtelibModelInitialization.cpp
    Algos/SgtelibModel/SgtelibModelIteration.cpp
    Algos/SgtelibModel/SgtelibModelMegaIteration.cpp
    Algos/SgtelibModel/SgtelibModelOptimize.cpp
    Algos/SgtelibModel/SgtelibModelUpdate.cpp)

# require openMP package
if(NOMAD_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# build libraries

# Some headers are dynamically generated and have to be taken into account before
add_custom_command(
    PRE_BUILD
    OUTPUT cacheAttributesDefinition.hpp displayAttributesDefinition.hpp evalAttributesDefinition.hpp evaluatorControlAttributesDefinition.hpp pbAttributesDefinition.hpp runAttributesDefinition.hpp
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/cacheAttributesDefinition.txt ${CMAKE_CURRENT_BINARY_DIR}/cacheAttributesDefinition.txt
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/displayAttributesDefinition.txt ${CMAKE_CURRENT_BINARY_DIR}/displayAttributesDefinition.txt
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/evalAttributesDefinition.txt ${CMAKE_CURRENT_BINARY_DIR}/evalAttributesDefinition.txt
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/evaluatorControlAttributesDefinition.txt ${CMAKE_CURRENT_BINARY_DIR}/evaluatorControlAttributesDefinition.txt
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/pbAttributesDefinition.txt ${CMAKE_CURRENT_BINARY_DIR}/pbAttributesDefinition.txt
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinition.txt ${CMAKE_CURRENT_BINARY_DIR}/runAttributesDefinition.txt
    COMMAND WriteAttributeDefinitionFile
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/cacheAttributesDefinition.hpp ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/cacheAttributesDefinition.hpp
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/displayAttributesDefinition.hpp ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/displayAttributesDefinition.hpp
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/evalAttributesDefinition.hpp ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/evalAttributesDefinition.hpp
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/evaluatorControlAttributesDefinition.hpp ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/evaluatorControlAttributesDefinition.hpp
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/pbAttributesDefinition.hpp ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/pbAttributesDefinition.hpp
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/runAttributesDefinition.hpp ${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinition.hpp
    DEPENDS WriteAttributeDefinitionFile)

add_custom_target(write_attributes
    DEPENDS cacheAttributesDefinition.hpp displayAttributesDefinition.hpp evalAttributesDefinition.hpp evaluatorControlAttributesDefinition.hpp pbAttributesDefinition.hpp runAttributesDefinition.hpp)

# Utils
add_library (nomadUtils SHARED ${MATH_SOURCES} ${OUTPUT_SOURCES} ${PARAM_SOURCES} ${TYPE_SOURCES} ${UTIL_SOURCES})
target_include_directories(nomadUtils PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Math
    ${CMAKE_CURRENT_SOURCE_DIR}/Output
    ${CMAKE_CURRENT_SOURCE_DIR}/Param
    ${CMAKE_CURRENT_SOURCE_DIR}/Type
    ${CMAKE_CURRENT_SOURCE_DIR}/Util>
    $<INSTALL_INTERFACE:include/nomad>)
add_dependencies(nomadUtils write_attributes)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nomadUtils PUBLIC OpenMP::OpenMP_CXX)
endif()

# Evals
add_library (nomadEval SHARED ${CACHE_SOURCES} ${EVAL_SOURCES}) 
target_include_directories(nomadEval PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Cache
    ${CMAKE_CURRENT_SOURCE_DIR}/Eval>
    $<INSTALL_INTERFACE:include/nomad>)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nomadEval PUBLIC nomadUtils OpenMP::OpenMP_CXX)
else()
    target_link_libraries(nomadEval nomadUtils)
endif()

# Algos
add_library (nomadAlgos SHARED ${ALGOS_SOURCES} ${MADS_SOURCES} ${LH_SOURCES} ${NM_SOURCES} ${PHASEONE_SOURCES} ${QUADMODEL_SOURCES} ${SGTELIBMODEL_SOURCES}) 
target_include_directories(nomadAlgos PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Algos 
    ${CMAKE_CURRENT_SOURCE_DIR}/Algos/LatinHypercubeSampling
    ${CMAKE_CURRENT_SOURCE_DIR}/Algos/Mads
    ${CMAKE_CURRENT_SOURCE_DIR}/Algos/NelderMead
    ${CMAKE_CURRENT_SOURCE_DIR}/Algos/PhaseOne
    ${CMAKE_CURRENT_SOURCE_DIR}/Algos/QuadModel
    ${CMAKE_CURRENT_SOURCE_DIR}/Algos/SgtelibModel>
    $<INSTALL_INTERFACE:include/nomad>)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nomadAlgos PUBLIC nomadUtils nomadEval sgtelib PRIVATE OpenMP::OpenMP_CXX)#PRIVATE
else()
    target_link_libraries(nomadAlgos nomadUtils nomadEval sgtelib)
endif()

# build executables

# build attribute executable
add_executable(WriteAttributeDefinitionFile ${ATTRIBUTE_SOURCES})
target_include_directories(WriteAttributeDefinitionFile PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Attribute>
    $<INSTALL_INTERFACE:include/nomad>)

add_executable (nomad ${NOMAD_SOURCES}) 
target_include_directories(nomad PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ 
    ${CMAKE_CURRENT_SOURCE_DIR}/Nomad>
    $<INSTALL_INTERFACE:include/nomad>)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nomad PUBLIC nomadAlgos nomadUtils nomadEval INTERFACE OpenMP::OpenMP_CXX)#PRIVATE
else()
    target_link_libraries(nomad PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# external C interface
add_library (nomadInterface SHARED ${INTERFACE_SOURCES})
target_include_directories(nomadInterface PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Interfaces>
    $<INSTALL_INTERFACE:include/nomad>
    )
if(OpenMP_CXX_FOUND)
    target_link_libraries(nomadInterface PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
else()
    target_link_libraries(nomadInterface PUBLIC nomadAlgos nomadUtils nomadEval)
endif()

# installing executables and libraries
install(TARGETS nomad WriteAttributeDefinitionFile nomadUtils nomadEval nomadAlgos nomadInterface
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad)

# # installing headers

# General
install(FILES ${ALGOS_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos)

install(FILES ${ATTRIBUTE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Attribute)

install(FILES ${CACHE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Cache)

install(FILES ${EVAL_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Eval)

install(FILES ${INTERFACE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Interfaces)

install(FILES ${MATH_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Math)

install(FILES ${NOMAD_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Nomad)

install(FILES ${OUTPUT_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Output)

install(FILES ${PARAM_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Param)

install(FILES ${TYPE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Type)

install(FILES ${UTIL_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Util)

# Submodules
install(FILES ${LH_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/LatinHypercubeSampling)

install(FILES ${MADS_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/Mads)

install(FILES ${NM_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/NelderMead)

install(FILES ${PHASEONE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/PhaseOne)

install(FILES ${QUADMODEL_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/QuadModel)

install(FILES ${SGTELIBMODEL_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/SgtelibModel)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/nomad_nsbegin.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nomad_nsend.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nomad_version.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/)
