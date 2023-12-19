# Application modules

#
# Algos
#
set(ALGOS_HEADERS
Algos/Algorithm.hpp
Algos/AlgoStopReasons.hpp
Algos/CacheInterface.hpp
Algos/EvcInterface.hpp
Algos/Initialization.hpp
Algos/Iteration.hpp
Algos/IterationUtils.hpp
Algos/MainStep.hpp
Algos/MegaIteration.hpp
Algos/Step.hpp
Algos/Subproblem.hpp
Algos/SubproblemManager.hpp
Algos/SurrogateEvaluation.hpp
Algos/SurrogateEvaluator.hpp
Algos/Termination.hpp
Algos/TrialPointStats.hpp
)

set(ALGOS_SOURCES
Algos/Algorithm.cpp
Algos/AlgoStopReasons.cpp
Algos/CacheInterface.cpp
Algos/EvcInterface.cpp
Algos/Initialization.cpp
Algos/Iteration.cpp
Algos/IterationUtils.cpp
Algos/MainStep.cpp
Algos/MegaIteration.cpp
Algos/Step.cpp
Algos/Subproblem.cpp
Algos/SubproblemManager.cpp
Algos/SurrogateEvaluation.cpp
Algos/Termination.cpp
Algos/TrialPointStats.cpp)

#
# Attribute
#
set(ATTRIBUTE_TEXT
Attribute/cacheAttributesDefinition.txt
Attribute/deprecatedAttributesDefinition.txt
Attribute/displayAttributesDefinition.txt
Attribute/evalAttributesDefinition.txt
Attribute/evaluatorControlAttributesDefinition.txt
Attribute/evaluatorControlGlobalAttributesDefinition.txt
Attribute/pbAttributesDefinition.txt
Attribute/runAttributesDefinition.txt
Attribute/runAttributesDefinitionIBEX.txt
Attribute/runAttributesDefinitionNM.txt
Attribute/runAttributesDefinitionCS.txt
Attribute/runAttributesDefinitionPSDSSD.txt
Attribute/runAttributesDefinitionQPSolver.txt
Attribute/runAttributesDefinitionQuadModel.txt
Attribute/runAttributesDefinitionSgtelibModel.txt
Attribute/runAttributesDefinitionVNS.txt
Attribute/runAttributesDefinitionDisco.txt
)

set(ATTRIBUTE_HEADERS
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/cacheAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/deprecatedAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/displayAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/evalAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/evaluatorControlAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/evaluatorControlGlobalAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/pbAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinition.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinitionIBEX.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinitionNM.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinitionCS.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinitionPSDSSD.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinitionQuadModel.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinitionSgtelibModel.hpp
${CMAKE_CURRENT_SOURCE_DIR}/Attribute/runAttributesDefinitionVNS.hpp)

set(ATTRIBUTE_HEADERS_GENERATOR
Attribute/WriteAttributeDefinitionFile.cpp
)

#
# Cache
#
set(CACHE_HEADERS
Cache/CacheBase.hpp
Cache/CacheSet.hpp
)

set(CACHE_SOURCES
Cache/CacheBase.cpp
Cache/CacheSet.cpp
)

#
# Eval
#
set(EVAL_HEADERS
Eval/BarrierBase.hpp
Eval/BBInput.hpp
Eval/BBOutput.hpp
Eval/ComparePriority.hpp
Eval/ComputeSuccessType.hpp
Eval/Eval.hpp
Eval/EvalPoint.hpp
Eval/EvalQueuePoint.hpp
Eval/Evaluator.hpp
Eval/EvaluatorControl.hpp
Eval/EvcMainThreadInfo.hpp
Eval/MeshBase.hpp
Eval/ProgressiveBarrier.hpp
Eval/SuccessStats.hpp)

set(EVAL_SOURCES
Eval/BarrierBase.cpp
Eval/BBInput.cpp
Eval/BBOutput.cpp
Eval/ComparePriority.cpp
Eval/ComputeSuccessType.cpp
Eval/Eval.cpp
Eval/EvalPoint.cpp
Eval/EvalQueuePoint.cpp
Eval/Evaluator.cpp
Eval/EvaluatorControl.cpp
Eval/EvcMainThreadInfo.cpp
Eval/MeshBase.cpp
Eval/ProgressiveBarrier.cpp
Eval/SuccessStats.cpp
)

#
# Math
#
set(MATH_HEADERS
Math/ArrayOfDouble.hpp
Math/ArrayOfPoint.hpp
Math/Direction.hpp
Math/Double.hpp
Math/LHS.hpp
Math/MatrixUtils.hpp
Math/Point.hpp
Math/RandomPickup.hpp
Math/RNG.hpp)

set(MATH_SOURCES
Math/ArrayOfDouble.cpp
Math/ArrayOfPoint.cpp
Math/Direction.cpp
Math/Double.cpp
Math/LHS.cpp
Math/MatrixUtils.cpp
Math/Point.cpp
Math/RandomPickup.cpp
Math/RNG.cpp
)

#
# Nomad
#
set(NOMAD_HEADERS
Nomad/nomad.hpp
)

set(NOMAD_SOURCES
Nomad/nomad.cpp
)

#
# Output
#
set(OUTPUT_HEADERS
Output/OutputDirectToFile.hpp
Output/OutputInfo.hpp
Output/OutputQueue.hpp
Output/StatsInfo.hpp
)

set(OUTPUT_SOURCES
Output/OutputDirectToFile.cpp
Output/OutputInfo.cpp
Output/OutputQueue.cpp
Output/StatsInfo.cpp
)

#
# Param
#
set(PARAM_HEADERS
Param/AllParameters.hpp
Param/Attribute.hpp
Param/AttributeFactory.hpp
Param/CacheParameters.hpp
Param/DeprecatedParameters.hpp
Param/DisplayParameters.hpp
Param/EvalParameters.hpp
Param/EvaluatorControlParameters.hpp
Param/EvaluatorControlGlobalParameters.hpp
Param/ParameterEntries.hpp
Param/ParameterEntry.hpp
Param/Parameters.hpp
Param/ParametersNomad3.hpp
Param/PbParameters.hpp
Param/RunParameters.hpp
Param/TypeAttribute.hpp
)

set(PARAM_SOURCES
Param/AllParameters.cpp
Param/Attribute.cpp
Param/AttributeFactory.cpp
Param/CacheParameters.cpp
Param/DeprecatedParameters.cpp
Param/DisplayParameters.cpp
Param/EvalParameters.cpp
Param/EvaluatorControlParameters.cpp
Param/EvaluatorControlGlobalParameters.cpp
Param/ParameterEntries.cpp
Param/ParameterEntry.cpp
Param/Parameters.cpp
Param/ParametersNomad3.cpp
Param/PbParameters.cpp
Param/RunParameters.cpp
Param/TypeAttribute.cpp
)

#
# Type
#
set(TYPE_HEADERS
Type/BBInputType.hpp
Type/BBOutputType.hpp
Type/CallbackType.hpp
Type/CompareType.hpp
Type/ComputeType.hpp
Type/DirectionType.hpp
Type/EvalSortType.hpp
Type/EvalType.hpp
Type/LHSearchType.hpp
Type/ListOfVariableGroup.hpp
Type/SgtelibModelFeasibilityType.hpp
Type/SgtelibModelFormulationType.hpp
Type/StepType.hpp
)

set(TYPE_SOURCES
Type/BBInputType.cpp
Type/BBOutputType.cpp
Type/CallbackType.cpp
Type/CompareType.cpp
Type/ComputeType.cpp
Type/DirectionType.cpp
Type/EvalSortType.cpp
Type/EvalType.cpp
Type/LHSearchType.cpp
Type/ListOfVariableGroup.cpp
Type/SgtelibModelFeasibilityType.cpp
Type/SgtelibModelFormulationType.cpp
Type/StepType.cpp
)

#
# Util
#
set(UTIL_HEADERS
Util/AllStopReasons.hpp
Util/ArrayOfString.hpp
Util/Clock.hpp
Util/defines.hpp
Util/Exception.hpp
Util/fileutils.hpp
Util/MicroSleep.hpp
Util/StopReason.hpp
Util/Uncopyable.hpp
Util/utils.hpp
)

set(UTIL_SOURCES
Util/AllStopReasons.cpp
Util/ArrayOfString.cpp
Util/Clock.cpp
Util/defines.cpp
Util/Exception.cpp
Util/fileutils.cpp
Util/StopReason.cpp
Util/Uncopyable.cpp
Util/utils.cpp)

#
# CS
#
set(CS_HEADERS
Algos/CoordinateSearch/CS.hpp
Algos/CoordinateSearch/CSInitialization.hpp
Algos/CoordinateSearch/CSIteration.hpp
Algos/CoordinateSearch/CSMegaIteration.hpp
Algos/CoordinateSearch/CSMesh.hpp
Algos/CoordinateSearch/CSPoll.hpp
Algos/CoordinateSearch/CSPollMethod.hpp
Algos/CoordinateSearch/CSUpdate.hpp
)

set(CS_SOURCES
Algos/CoordinateSearch/CS.cpp
Algos/CoordinateSearch/CSInitialization.cpp
Algos/CoordinateSearch/CSIteration.cpp
Algos/CoordinateSearch/CSMegaIteration.cpp
Algos/CoordinateSearch/CSMesh.cpp
Algos/CoordinateSearch/CSPoll.cpp
Algos/CoordinateSearch/CSPollMethod.cpp
Algos/CoordinateSearch/CSUpdate.cpp
)

#
# LH
#
set(LH_HEADERS
Algos/LatinHypercubeSampling/LH.hpp
)

set(LH_SOURCES
Algos/LatinHypercubeSampling/LH.cpp
)

#
# DMultiMads
#
set(DMULTIMADS_HEADERS
Algos/DMultiMads/DMultiMads.hpp
Algos/DMultiMads/DMultiMadsBarrier.hpp
Algos/DMultiMads/DMultiMadsIteration.hpp
Algos/DMultiMads/DMultiMadsMegaIteration.hpp
Algos/DMultiMads/DMultiMadsUpdate.hpp
)

set(DMULTIMADS_SOURCES
Algos/DMultiMads/DMultiMads.cpp
Algos/DMultiMads/DMultiMadsBarrier.cpp
Algos/DMultiMads/DMultiMadsIteration.cpp
Algos/DMultiMads/DMultiMadsMegaIteration.cpp
Algos/DMultiMads/DMultiMadsUpdate.cpp
)

#
# DiscoMads
#
set(DISCOMADS_HEADERS
Algos/DiscoMads/DiscoMads.hpp
Algos/DiscoMads/DiscoMadsMegaIteration.hpp
Algos/DiscoMads/DiscoMadsIteration.hpp
Algos/DiscoMads/RevealingPoll.hpp
Algos/DiscoMads/DiscoMadsBarrier.hpp
Algos/DiscoMads/DiscoMadsUpdate.hpp
)

set(DISCOMADS_SOURCES
Algos/DiscoMads/DiscoMads.cpp
Algos/DiscoMads/DiscoMadsMegaIteration.cpp
Algos/DiscoMads/DiscoMadsIteration.cpp
Algos/DiscoMads/RevealingPoll.cpp
Algos/DiscoMads/DiscoMadsBarrier.cpp
Algos/DiscoMads/DiscoMadsUpdate.cpp
)
#
# Mads
#
set(MADS_HEADERS
Algos/Mads/DoublePollMethod.hpp
Algos/Mads/GMesh.hpp
Algos/Mads/LHSearchMethod.hpp
Algos/Mads/Mads.hpp
Algos/Mads/MadsInitialization.hpp
Algos/Mads/MadsIteration.hpp
Algos/Mads/MadsMegaIteration.hpp
Algos/Mads/MadsUpdate.hpp
Algos/Mads/MegaSearchPoll.hpp
Algos/Mads/NMSearchMethod.hpp
Algos/Mads/NP1UniPollMethod.hpp
Algos/Mads/Ortho2NPollMethod.hpp
Algos/Mads/OrthoNPlus1PollMethod.hpp
Algos/Mads/Poll.hpp
Algos/Mads/PollMethodBase.hpp
Algos/Mads/QPSolverAlgoSearchMethod.hpp
Algos/Mads/QuadSearchMethod.hpp
Algos/Mads/Search.hpp
Algos/Mads/SearchMethodAlgo.hpp
Algos/Mads/SearchMethodBase.hpp
Algos/Mads/SearchMethodSimple.hpp
Algos/Mads/SgtelibSearchMethod.hpp
Algos/Mads/SimpleLineSearchMethod.hpp
Algos/Mads/SinglePollMethod.hpp
Algos/Mads/SpeculativeSearchMethod.hpp
Algos/Mads/TemplateAlgoSearchMethod.hpp
Algos/Mads/TemplateSimpleSearchMethod.hpp
Algos/Mads/UserSearchMethod.hpp
Algos/Mads/VNSmartSearchMethod.hpp
Algos/Mads/VNSSearchMethod.hpp
)

set(MADS_SOURCES
Algos/Mads/DoublePollMethod.cpp
Algos/Mads/GMesh.cpp
Algos/Mads/LHSearchMethod.cpp
Algos/Mads/Mads.cpp
Algos/Mads/MadsInitialization.cpp
Algos/Mads/MadsIteration.cpp
Algos/Mads/MadsMegaIteration.cpp
Algos/Mads/MadsUpdate.cpp
Algos/Mads/MegaSearchPoll.cpp
Algos/Mads/NMSearchMethod.cpp
Algos/Mads/NP1UniPollMethod.cpp
Algos/Mads/Ortho2NPollMethod.cpp
Algos/Mads/OrthoNPlus1PollMethod.cpp
Algos/Mads/Poll.cpp
Algos/Mads/PollMethodBase.cpp
Algos/Mads/QPSolverAlgoSearchMethod.cpp
Algos/Mads/QuadSearchMethod.cpp
Algos/Mads/Search.cpp
Algos/Mads/SearchMethodAlgo.cpp
Algos/Mads/SearchMethodBase.cpp
Algos/Mads/SearchMethodSimple.cpp
Algos/Mads/SgtelibSearchMethod.cpp
Algos/Mads/SimpleLineSearchMethod.cpp
Algos/Mads/SinglePollMethod.cpp
Algos/Mads/SpeculativeSearchMethod.cpp
Algos/Mads/TemplateAlgoSearchMethod.cpp
Algos/Mads/TemplateSimpleSearchMethod.cpp
Algos/Mads/UserSearchMethod.cpp
Algos/Mads/VNSmartSearchMethod.cpp
Algos/Mads/VNSSearchMethod.cpp
)

#
# NM
#
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
Algos/NelderMead/NMUpdate.hpp
)

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

#
# PhaseOne
#
set(PHASEONE_HEADERS
Algos/PhaseOne/PhaseOne.hpp
)

set(PHASEONE_SOURCES
Algos/PhaseOne/PhaseOne.cpp
)

#
# PSD Mads
#
if(OpenMP_CXX_FOUND)
   set(PSDMADS_HEADERS
   Algos/PSDMads/PSDMads.hpp
   Algos/PSDMads/PSDMadsMegaIteration.hpp
   )
   set(PSDMADS_SOURCES
   Algos/PSDMads/PSDMads.cpp
   Algos/PSDMads/PSDMadsMegaIteration.cpp
   )
endif()


#
# QuadModel
#
set(QUADMODEL_HEADERS
Algos/QuadModel/QuadModelAlgo.hpp
Algos/QuadModel/QuadModelEvaluator.hpp
Algos/QuadModel/QuadModelInitialization.hpp
Algos/QuadModel/QuadModelIteration.hpp
Algos/QuadModel/QuadModelIterationUtils.hpp
Algos/QuadModel/QuadModelMegaIteration.hpp
Algos/QuadModel/QuadModelOptimize.hpp
Algos/QuadModel/QuadModelSinglePass.hpp
Algos/QuadModel/QuadModelUpdate.hpp
)

set(QUADMODEL_SOURCES
Algos/QuadModel/QuadModelAlgo.cpp
Algos/QuadModel/QuadModelEvaluator.cpp
Algos/QuadModel/QuadModelInitialization.cpp
Algos/QuadModel/QuadModelIteration.cpp
Algos/QuadModel/QuadModelIterationUtils.cpp
Algos/QuadModel/QuadModelMegaIteration.cpp
Algos/QuadModel/QuadModelOptimize.cpp
Algos/QuadModel/QuadModelSinglePass.cpp
Algos/QuadModel/QuadModelUpdate.cpp
)

#
# QP Solvers
#
set(QPSOLVER_HEADERS
Algos/QPSolverAlgo/QPSolverAlgo.hpp
Algos/QPSolverAlgo/QPSolverAlgoIteration.hpp
Algos/QPSolverAlgo/QPSolverAlgoMegaIteration.hpp
Algos/QPSolverAlgo/QPSolverOptimize.hpp
Algos/QPSolverAlgo/QPSolverAlgoSinglePass.hpp
)

set(QPSOLVER_SOURCES
Algos/QPSolverAlgo/QPSolverAlgo.cpp
Algos/QPSolverAlgo/QPSolverAlgoIteration.cpp
Algos/QPSolverAlgo/QPSolverAlgoMegaIteration.cpp
Algos/QPSolverAlgo/QPSolverOptimize.cpp
Algos/QPSolverAlgo/QPSolverAlgoSinglePass.cpp
)

#
# SgtelibModel
#
set(SGTELIBMODEL_HEADERS
Algos/SgtelibModel/SgtelibModel.hpp
Algos/SgtelibModel/SgtelibModelEvaluator.hpp
Algos/SgtelibModel/SgtelibModelFilterCache.hpp
Algos/SgtelibModel/SgtelibModelInitialization.hpp
Algos/SgtelibModel/SgtelibModelIteration.hpp
Algos/SgtelibModel/SgtelibModelMegaIteration.hpp
Algos/SgtelibModel/SgtelibModelOptimize.hpp
Algos/SgtelibModel/SgtelibModelUpdate.hpp
)

set(SGTELIBMODEL_SOURCES
Algos/SgtelibModel/SgtelibModel.cpp
Algos/SgtelibModel/SgtelibModelEvaluator.cpp
Algos/SgtelibModel/SgtelibModelFilterCache.cpp
Algos/SgtelibModel/SgtelibModelInitialization.cpp
Algos/SgtelibModel/SgtelibModelIteration.cpp
Algos/SgtelibModel/SgtelibModelMegaIteration.cpp
Algos/SgtelibModel/SgtelibModelOptimize.cpp
Algos/SgtelibModel/SgtelibModelUpdate.cpp
)

#
# SSD Mads
#
set(SSDMADS_HEADERS
Algos/SSDMads/SSDMads.hpp
Algos/SSDMads/SSDMadsMegaIteration.hpp
)

set(SSDMADS_SOURCES
Algos/SSDMads/SSDMads.cpp
Algos/SSDMads/SSDMadsMegaIteration.cpp
)

#
# Template Algorithm
# Can be used as an example to create a new algorithm
#
set(TEMPLATE_ALGO_HEADERS
Algos/TemplateAlgo/TemplateAlgo.hpp
Algos/TemplateAlgo/TemplateAlgoInitialization.hpp
Algos/TemplateAlgo/TemplateAlgoIteration.hpp
Algos/TemplateAlgo/TemplateAlgoMegaIteration.hpp
Algos/TemplateAlgo/TemplateAlgoRandom.hpp
Algos/TemplateAlgo/TemplateAlgoSinglePass.hpp
Algos/TemplateAlgo/TemplateAlgoUpdate.hpp
)

set(TEMPLATE_ALGO_SOURCES
Algos/TemplateAlgo/TemplateAlgo.cpp
Algos/TemplateAlgo/TemplateAlgoInitialization.cpp
Algos/TemplateAlgo/TemplateAlgoIteration.cpp
Algos/TemplateAlgo/TemplateAlgoMegaIteration.cpp
Algos/TemplateAlgo/TemplateAlgoRandom.cpp
Algos/TemplateAlgo/TemplateAlgoSinglePass.cpp
Algos/TemplateAlgo/TemplateAlgoUpdate.cpp)

#
# Simple Line Search
#
set(SIMPLELINESEARCH_ALGO_HEADERS
Algos/SimpleLineSearch/SimpleLineSearch.hpp
Algos/SimpleLineSearch/SimpleLineSearchMegaIteration.hpp
)

set(SIMPLELINESEARCH_ALGO_SOURCES
Algos/SimpleLineSearch/SimpleLineSearch.cpp
Algos/SimpleLineSearch/SimpleLineSearchMegaIteration.cpp
)


#
# VNS
#
set(VNSMADS_HEADERS
Algos/VNSMads/VNS.hpp
#Algos/VNSMads/VNSMegaIteration.cpp
)

set(VNSMADS_SOURCES
Algos/VNSMads/VNS.cpp
)


#
# Generate Attribute headers
#

# Build attribute executable
add_executable(
  WriteAttributeDefinitionFile ${ATTRIBUTE_HEADERS_GENERATOR}
)

# Command to create the attribute headers
add_custom_command(
  OUTPUT ${ATTRIBUTE_HEADERS}
  COMMAND WriteAttributeDefinitionFile
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Attribute
  DEPENDS ${ATTRIBUTE_TEXT}
)

#
# Build libraries and executable
# Install all
#

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# export all symbols on Windows
if (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Build shared libraries (nomadUtils, nomadEval, nomadAlgos)
include(./CMakeListsLibs.txt)

# build statically linked library (nomadStatic) for PyNomad interface
if(BUILD_INTERFACE_PYTHON MATCHES ON)
   message(STATUS "  Warning: build nomad static library for Python interface.")
   include(./CMakeListsLibsStatic.txt)
endif()

# Build dynamically linked executables (nomad)
include(./CMakeListsExes.txt)


# Install header files
install(
  FILES ${ALGOS_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos
)

install(
  FILES ${ATTRIBUTE_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Attribute
)

install(
  FILES ${CACHE_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Cache
)

install(
  FILES ${EVAL_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Eval
)

install(
  FILES ${MATH_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Math
)

install(
  FILES ${NOMAD_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Nomad
)

install(
  FILES ${OUTPUT_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Output
)

install(
  FILES ${PARAM_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Param
)

install(
  FILES ${TYPE_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Type
)

install(
  FILES ${UTIL_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Util
)

install(
  FILES ${CS_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/CoordinateSearch
)

install(
  FILES ${DMULTIMADS_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/DMultiMads
)

install(
   FILES ${DISCOMADS_HEADERS}
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/DiscoMads
)

install(
  FILES ${LH_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/LatinHypercubeSampling
)

install(
  FILES ${MADS_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/Mads
)

install(
  FILES ${NM_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/NelderMead
)

install(
  FILES ${PHASEONE_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/PhaseOne
)

install(
  FILES ${PSDMADS_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/PSDMads
)

install(
  FILES ${QPSOLVER_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/QPSolverAlgo
)

install(
  FILES ${QUADMODEL_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/QuadModel
)

install(
  FILES ${SGTELIB_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/SgtelibModel
)

install(
  FILES ${SSDMADS_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/SSDMads
)

install(
  FILES ${TEMPLATE_ALGO_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/TemplateAlgo
)

install(
  FILES ${SIMPLELINESEARCH_ALGO_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/SimpleLineSearch
)

install(
  FILES ${VNSMADS_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/Algos/VNSMads
)

install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/nomad_nsbegin.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/nomad_nsend.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/nomad_platform.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/nomad_version.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nomad/
)
