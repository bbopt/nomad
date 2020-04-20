
ifndef ($(VARIANT))
VARIANT             = release
endif

ifndef PAR
PAR                 = 8
endif


UNAME := $(shell uname)

TOP                 = $(abspath ..)
BUILD_DIR           = $(TOP)/build/$(VARIANT)
SRC_DIR             = $(TOP)/src
EXT_DIR             = $(TOP)/ext
SGTE_DIR            = $(EXT_DIR)/sgtelib
SGTE_SRC_DIR        = $(SGTE_DIR)/src
SGTE_OBJ_DIR        = $(SGTE_DIR)/build
SGTE_LIB_DIR        = $(SGTE_DIR)/lib
OBJ_DIR             = $(BUILD_DIR)/obj
INCLUDE_DIR         = $(BUILD_DIR)/include
SGTE_INCLUDE_DIR    = $(INCLUDE_DIR)/sgtelib
LIB_DIR             = $(BUILD_DIR)/lib
BIN_DIR             = $(BUILD_DIR)/bin
EXE                 = $(BIN_DIR)/nomad

# Evaluate if SGTE_DIR exists or not.
ifeq "$(wildcard $(SGTE_DIR))" ""
    ifeq ($(USE_SGTELIB), 0)
        $(info Compile without sgtelib)
    else
        $(info Warning: Directory $(SGTE_DIR) does not exist.)
        $(info Here are your options:)
        $(info    1- Get the Sgtelib library at: https://github.com/bbopt/sgtelib/zipball/master)
        $(info       and unzip it as directory $(SGTE_DIR) )
        $(info )
        $(info    2- Developer: Follow git module instructions in $(NOMAD_HOME)/README_GITHUB.txt )
        $(info )
        $(info    3- Compile without the Sgtelib library, by using make argument USE_SGTELIB=0 )
        $(info )
        ifdef USE_SGTELIB
            $(error USE_SGTELIB is set but directory $(SGTE_DIR) does not exist)
        endif
    endif
else
    USE_SGTELIB = 1
endif


# Main modules
ALGO_DIRNAME        = Algos
ATTR_DIRNAME        = Attribute
CACHE_DIRNAME       = Cache
EVAL_DIRNAME        = Eval
MATH_DIRNAME        = Math
NOMAD_DIRNAME       = Nomad
OUTPUT_DIRNAME      = Output
PARAM_DIRNAME       = Param
TYPE_DIRNAME        = Type
UTIL_DIRNAME        = Util

# Sub modules
LH_DIRNAME          = $(ALGO_DIRNAME)/LatinHypercubeSampling
MADS_DIRNAME        = $(ALGO_DIRNAME)/Mads
NM_DIRNAME          = $(ALGO_DIRNAME)/NelderMead
PHASEONE_DIRNAME    = $(ALGO_DIRNAME)/PhaseOne
SGTELIBMODEL_DIRNAME = $(ALGO_DIRNAME)/SgtelibModel
QUADMODEL_DIRNAME   = $(ALGO_DIRNAME)/QuadModel

ALGO_COMPONENT      = $(SRC_DIR)/$(ALGO_DIRNAME)
ATTR_COMPONENT      = $(SRC_DIR)/$(ATTR_DIRNAME)
CACHE_COMPONENT     = $(SRC_DIR)/$(CACHE_DIRNAME)
EVAL_COMPONENT      = $(SRC_DIR)/$(EVAL_DIRNAME)
MADS_COMPONENT      = $(SRC_DIR)/$(MADS_DIRNAME)
MATH_COMPONENT      = $(SRC_DIR)/$(MATH_DIRNAME)
LH_COMPONENT        = $(SRC_DIR)/$(LH_DIRNAME)
NM_COMPONENT        = $(SRC_DIR)/$(NM_DIRNAME)
NOMAD_COMPONENT     = $(SRC_DIR)/$(NOMAD_DIRNAME)
OUTPUT_COMPONENT    = $(SRC_DIR)/$(OUTPUT_DIRNAME)
PARAM_COMPONENT     = $(SRC_DIR)/$(PARAM_DIRNAME)
PHASEONE_COMPONENT  = $(SRC_DIR)/$(PHASEONE_DIRNAME)
SGTELIBMODEL_COMPONENT = $(SRC_DIR)/$(SGTELIBMODEL_DIRNAME)
QUADMODEL_COMPONENT = $(SRC_DIR)/$(QUADMODEL_DIRNAME)
TYPE_COMPONENT      = $(SRC_DIR)/$(TYPE_DIRNAME)
UTIL_COMPONENT      = $(SRC_DIR)/$(UTIL_DIRNAME)
COMPONENTS          = $(ATTR_DIRNAME) $(CACHE_DIRNAME) $(EVAL_DIRNAME) \
                      $(LH_DIRNAME) $(MADS_DIRNAME) $(MATH_DIRNAME) \
                      $(NM_DIRNAME) $(NOMAD_DIRNAME) $(OUTPUT_DIRNAME) \
                      $(PARAM_DIRNAME) $(SGTELIBMODEL_DIRNAME) $(QUADMODEL_DIRNAME) \
                      $(TYPE_DIRNAME) $(UTIL_DIRNAME)
INCLUDE_COMPONENTS  := $(addprefix $(INCLUDE_DIR)/,$(COMPONENTS))

# TODO: Remove obj file names from here
ATTR_EXE            = WriteAttributeDefinitionFile
ATTR_EXE            := $(addprefix $(BIN_DIR)/,$(ATTR_EXE))
CACHE_OBJ           = CacheBase.o CacheSet.o
CACHE_OBJ           := $(addprefix $(OBJ_DIR)/,$(CACHE_OBJ))
EVAL_OBJ            = Barrier.o BBInput.o BBOutput.o CallbackType.o \
                      Eval.o EvalPoint.o \
                      EvalQueuePoint.o EvaluatorControl.o Evaluator.o
EVAL_OBJ            := $(addprefix $(OBJ_DIR)/,$(EVAL_OBJ))
ALGO_OBJ            = Algorithm.o CacheInterface.o EvcInterface.o Initialization.o  \
                      Iteration.o IterationUtils.o MainStep.o MegaIteration.o MeshBase.o \
                      Projection.o Step.o Subproblem.o \
                      Termination.o
ALGO_OBJ            := $(addprefix $(OBJ_DIR)/,$(ALGO_OBJ))
LH_OBJ              = LH.o
LH_OBJ              := $(addprefix $(OBJ_DIR)/,$(LH_OBJ))
MADS_OBJ            = GMesh.o LHSearchMethod.o \
                      Mads.o MadsInitialization.o MadsIteration.o \
                      MadsMegaIteration.o MegaSearchPoll.o MadsUpdate.o \
                      NMSearchMethod.o Poll.o Search.o SearchMethodAlgo.o SearchMethodBase.o \
                      SearchMethodSimple.o SgtelibSearchMethod.o QuadSearchMethod.o \
                      SpeculativeSearchMethod.o UserSearchMethod.o
MADS_OBJ            := $(addprefix $(OBJ_DIR)/,$(MADS_OBJ))
MATH_OBJ            = ArrayOfDouble.o ArrayOfPoint.o Direction.o Double.o LHS.o \
                      MatrixUtils.o Point.o RNG.o
MATH_OBJ            := $(addprefix $(OBJ_DIR)/,$(MATH_OBJ))
NM_OBJ              = NM.o NMAllReflective.o NMInitialization.o \
                      NMInitializeSimplex.o NMIteration.o NMIterationUtils.o \
                      NMMegaIteration.o NMReflective.o NMShrink.o \
                      NMSimplexEvalPoint.o NMUpdate.o
NM_OBJ              := $(addprefix $(OBJ_DIR)/,$(NM_OBJ))
OUTPUT_OBJ          = OutputInfo.o OutputQueue.o StatsInfo.o
OUTPUT_OBJ          := $(addprefix $(OBJ_DIR)/,$(OUTPUT_OBJ))
PARAM_OBJ           = AllParameters.o \
                      CacheParameters.o DisplayParameters.o \
                      EvaluatorControlParameters.o EvalParameters.o \
                      ParameterEntry.o ParameterEntries.o \
                      Parameters.o ParametersNomad3.o PbParameters.o \
                      RunParameters.o
PARAM_OBJ           := $(addprefix $(OBJ_DIR)/,$(PARAM_OBJ))
PHASEONE_OBJ        = PhaseOne.o 
PHASEONE_OBJ        := $(addprefix $(OBJ_DIR)/,$(PHASEONE_OBJ)) 
SGTELIBMODEL_OBJ    = SgtelibModel.o SgtelibModelEvaluator.o \
                      SgtelibModelFilterCache.o \
                      SgtelibModelInitialization.o \
                      SgtelibModelIteration.o \
                      SgtelibModelMegaIteration.o \
                      SgtelibModelOptimize.o \
                      SgtelibModelUpdate.o
SGTELIBMODEL_OBJ    := $(addprefix $(OBJ_DIR)/,$(SGTELIBMODEL_OBJ))
QUADMODEL_OBJ       = QuadModelAlgo.o QuadModelEvaluator.o \
                      QuadModelInitialization.o \
                      QuadModelIteration.o \
                      QuadModelIterationUtils.o \
                      QuadModelMegaIteration.o \
	              QuadModelOptimize.o \
                      QuadModelSinglePass.o \
                      QuadModelUpdate.o
QUADMODEL_OBJ       := $(addprefix $(OBJ_DIR)/,$(QUADMODEL_OBJ))
TYPE_OBJ            = BBInputType.o BBOutputType.o EvalType.o LHSearchType.o \
                      SgtelibModelFeasibilityType.o SgtelibModelFormulationType.o
TYPE_OBJ            := $(addprefix $(OBJ_DIR)/,$(TYPE_OBJ))
UTIL_OBJ            = ArrayOfString.o Clock.o Exception.o fileutils.o \
                      StopReason.o utils.o
UTIL_OBJ            := $(addprefix $(OBJ_DIR)/,$(UTIL_OBJ))

SGTE_OBJ            = Kernel.o Matrix.o Metrics.o sgtelib_help.o sgtelib.o \
                      Surrogate_CN.o Surrogate_Ensemble.o Surrogate_Factory.o \
                      Surrogate_Kriging.o Surrogate_KS.o Surrogate_LOWESS.o \
                      Surrogate.o Surrogate_Parameters.o Surrogate_PRS_CAT.o \
                      Surrogate_PRS_EDGE.o Surrogate_PRS.o Surrogate_RBF.o \
                      Surrogate_Utils.o Tests.o TrainingSet.o
SGTE_OBJ            := $(addprefix $(SGTE_OBJ_DIR)/,$(SGTE_OBJ))

UTILS_OBJ_LIB       = $(MATH_OBJ) $(OUTPUT_OBJ) $(PARAM_OBJ) $(TYPE_OBJ) $(UTIL_OBJ)
EVAL_OBJ_LIB        = $(CACHE_OBJ) $(EVAL_OBJ)
ALGOS_OBJ_LIB       = $(ALGO_OBJ) $(MADS_OBJ) $(LH_OBJ) $(NM_OBJ) \
                      $(PHASEONE_OBJ)
ifeq ($(USE_SGTELIB), 1)
ALGOS_OBJ_LIB       += $(SGTELIBMODEL_OBJ) $(QUADMODEL_OBJ)
endif
OBJ_LIB             = $(UTILS_OBJ_LIB) $(EVAL_OBJ_LIB) $(ALGOS_OBJ_LIB)


# TODO can we get the version number from somewhere else?
# NOMAD_MAJOR_VER := $(shell awk '$2=="MAJOR_VER"{print $$3}' nomad_version.hpp)
UTILS_LIB_CURRENT_VERSION = 4.0.0
EVAL_LIB_CURRENT_VERSION = 4.0.0
ALGOS_LIB_CURRENT_VERSION = 4.0.0
NOMAD_CURRENT_VERSION = 4.0.0

UTILS_NAME_AND_VERSION    = nomadUtils.$(UTILS_LIB_CURRENT_VERSION)
UTILS_LIB_DYNAMIC_NAME    = lib$(UTILS_NAME_AND_VERSION).so
UTILS_LIB_DYNAMIC         = $(LIB_DIR)/$(UTILS_LIB_DYNAMIC_NAME)

EVAL_NAME_AND_VERSION     = nomadEval.$(UTILS_LIB_CURRENT_VERSION)
EVAL_LIB_DYNAMIC_NAME     = lib$(EVAL_NAME_AND_VERSION).so
EVAL_LIB_DYNAMIC          = $(LIB_DIR)/$(EVAL_LIB_DYNAMIC_NAME)

ALGOS_NAME_AND_VERSION    = nomadAlgos.$(UTILS_LIB_CURRENT_VERSION)
ALGOS_LIB_DYNAMIC_NAME    = lib$(ALGOS_NAME_AND_VERSION).so
ALGOS_LIB_DYNAMIC         = $(LIB_DIR)/$(ALGOS_LIB_DYNAMIC_NAME)

ifeq ($(USE_SGTELIB), 1)
SGTE_NAME                 = sgtelib
SGTE_LIB_DYNAMIC_NAME     = lib$(SGTE_NAME).so
SGTE_LIB_DYNAMIC          = $(LIB_DIR)/$(SGTE_LIB_DYNAMIC_NAME)
LIB_DYNAMIC_SGTE          = -l$(SGTE_NAME)
else
SGTE_NAME                 = 
SGTE_LIB_DYNAMIC_NAME     = 
SGTE_LIB_DYNAMIC          = 
LIB_DYNAMIC_SGTE          =
endif

LIB_DYNAMIC_UTILS         = -l$(UTILS_NAME_AND_VERSION)
LIB_DYNAMIC_EVAL          = -l$(EVAL_NAME_AND_VERSION)
LIB_DYNAMIC_ALGOS         = -l$(ALGOS_NAME_AND_VERSION)

LIB_DYNAMIC               = $(UTILS_LIB_DYNAMIC) $(EVAL_LIB_DYNAMIC) \
                            $(ALGOS_LIB_DYNAMIC) $(SGTE_LIB_DYNAMIC)

LIB_DYNAMIC_LINK          = $(LIB_DYNAMIC_UTILS) $(LIB_DYNAMIC_EVAL) $(LIB_DYNAMIC_ALGOS)

ifeq ($(VARIANT), release)
CXXFLAGS            = -O2
else
CXXFLAGS            = -g -DDEBUG
endif

ifdef USE_UNORDEREDSET
CXXFLAGS            += -DUSE_UNORDEREDSET
endif
ifeq ($(USE_SGTELIB), 1)
CXXFLAGS            += -DUSE_SGTELIB
endif
ifeq ($(TIME_STATS), 1)
CXXFLAGS            += -DTIME_STATS
endif

CXXFLAGS            += -std=c++14 -Wall -Wextra -fpic
CXXFLAGS_UTILS_LIBS = -L$(LIB_DIR)
CXXFLAGS_EVAL_LIBS  = -L$(LIB_DIR)
CXXFLAGS_ALGOS_LIBS = -L$(LIB_DIR)

ifeq ($(UNAME), Darwin)
CXXFLAGS_UTILS_LIBS += -current_version $(UTILS_LIB_CURRENT_VERSION) -compatibility_version $(UTILS_LIB_CURRENT_VERSION) -install_name $(UTILS_LIB_DYNAMIC) -headerpad_max_install_names
CXXFLAGS_EVAL_LIBS  += -current_version $(EVAL_LIB_CURRENT_VERSION) -compatibility_version $(EVAL_LIB_CURRENT_VERSION) -install_name $(EVAL_LIB_DYNAMIC) -headerpad_max_install_names
CXXFLAGS_ALGOS_LIBS += -current_version $(ALGOS_LIB_CURRENT_VERSION) -compatibility_version $(ALGOS_LIB_CURRENT_VERSION) -install_name $(ALGOS_LIB_DYNAMIC) -headerpad_max_install_names
endif

ifeq ($(UNAME), Linux)
CXXFLAGS_UTILS_LIBS += -Wl,-soname,'$(UTILS_LIB_DYNAMIC_NAME)' -Wl,-rpath-link,'.' -Wl,-rpath,'$$ORIGIN'
CXXFLAGS_EVAL_LIBS  += -Wl,-soname,'$(EVAL_LIB_DYNAMIC_NAME)' -Wl,-rpath-link,'.' -Wl,-rpath,'$$ORIGIN'
CXXFLAGS_ALGOS_LIBS += -Wl,-soname,'$(ALGOS_LIB_DYNAMIC_NAME)' -Wl,-rpath-link,'.' -Wl,-rpath,'$$ORIGIN'
endif

OBJFLAGS            = -c

COMPILE             = g++ $(CXXFLAGS)

# Export all variables
export

all: check_sgte $(INCLUDE_DIR) $(OBJ_DIR) $(ATTR_EXE) $(OBJ_LIB) $(LIB_DYNAMIC) $(EXE)
ifeq ($(USE_SGTELIB), 1)
	@cd $(SGTE_SRC_DIR) && $(MAKE) -j$(PAR) all
endif
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all
	@cd $(CACHE_COMPONENT) && $(MAKE) -j$(PAR) all
	@cd $(EVAL_COMPONENT) && $(MAKE) -j$(PAR) all
	@cd $(MATH_COMPONENT) && $(MAKE) -j$(PAR) all
	@cd $(OUTPUT_COMPONENT) && $(MAKE) -j$(PAR) all
	@cd $(PARAM_COMPONENT) && $(MAKE) -j$(PAR) all
	@cd $(TYPE_COMPONENT) && $(MAKE) -j$(PAR) all
	@cd $(UTIL_COMPONENT) && $(MAKE) -j$(PAR) all


check_sgte:
ifndef USE_SGTELIB
	$(error Cannot compile)
endif

lib_only: $(INCLUDE_DIR) $(OBJ_DIR) $(OBJ_LIB) $(LIB_DYNAMIC)

$(OBJ_DIR): 
	@mkdir -p $(OBJ_DIR)

$(ATTR_EXE) :
	@cd $(ATTR_COMPONENT) && $(MAKE) all
	@echo

$(ALGO_OBJ) :
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all

$(CACHE_OBJ) :
	@cd $(CACHE_COMPONENT) && $(MAKE) -j$(PAR) all
	@echo

$(EVAL_OBJ) :
	@cd $(EVAL_COMPONENT) && $(MAKE) -j$(PAR) all
	@echo

$(MATH_OBJ) :
	@cd $(MATH_COMPONENT) && $(MAKE) -j$(PAR) all
	@echo

$(OUTPUT_OBJ) :
	@cd $(OUTPUT_COMPONENT) && $(MAKE) -j$(PAR) all
	@echo

$(PARAM_OBJ) :
	@cd $(PARAM_COMPONENT) && $(MAKE) -j$(PAR) all
	@echo

$(TYPE_OBJ) :
	@cd $(TYPE_COMPONENT) && $(MAKE) -j$(PAR) all
	@echo

$(UTIL_OBJ) :
	@cd $(UTIL_COMPONENT) && $(MAKE) -j$(PAR) all
	@echo

$(SGTE_OBJ) :
	@cd $(SGTE_SRC_DIR) && $(MAKE) -j$(PAR) all
	@echo

$(LH_OBJ) : $(ALGO_OBJ)
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all

$(MADS_OBJ) : $(ALGO_OBJ)
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all

$(NM_OBJ) : $(ALGO_OBJ)
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all

$(PHASEONE_OBJ) : $(ALGO_OBJ)
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all

$(SGTELIBMODEL_OBJ) : $(ALGO_OBJ)
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all

$(QUADMODEL_OBJ) : $(ALGO_OBJ)
	@cd $(ALGO_COMPONENT) && $(MAKE) -j$(PAR) all

$(EXE) : $(INCLUDE_DIR) $(ATTR_EXE) $(LIB_DYNAMIC) $(NOMAD_COMPONENT)/nomad.hpp $(NOMAD_COMPONENT)/nomad.cpp
	@cd $(NOMAD_COMPONENT) && $(MAKE) all LIB_DIR=$(LIB_DIR) LIB_DYNAMIC_LINK='$(LIB_DYNAMIC_LINK)'
	@echo


$(INCLUDE_DIR): $(SGTE_INCLUDE_DIR) nomad_nsbegin.hpp nomad_nsend.hpp nomad_version.hpp
	@mkdir -p $(INCLUDE_DIR)
	@cp -f *.hpp $(INCLUDE_DIR)

$(SGTE_INCLUDE_DIR):
ifeq ($(USE_SGTELIB), 1)
	@mkdir -p $(SGTE_INCLUDE_DIR)
	@cp $(SGTE_SRC_DIR)/*.hpp $(SGTE_INCLUDE_DIR)
endif


#
# Make the dynamic libraries
#

$(UTILS_LIB_DYNAMIC): $(ATTR_EXE) $(UTILS_OBJ_LIB) $(INCLUDE_DIR)
	@mkdir -p $(LIB_DIR)
	$(COMPILE) -shared -o $(UTILS_LIB_DYNAMIC) $(UTILS_OBJ_LIB) $(CXXFLAGS_UTILS_LIBS)

$(EVAL_LIB_DYNAMIC): $(ATTR_EXE) $(EVAL_OBJ_LIB) $(INCLUDE_DIR)
	@mkdir -p $(LIB_DIR)
	$(COMPILE) -shared -o $(EVAL_LIB_DYNAMIC) $(EVAL_OBJ_LIB) $(CXXFLAGS_EVAL_LIBS) $(LIB_DYNAMIC_UTILS)

$(ALGOS_LIB_DYNAMIC): $(ATTR_EXE) $(ALGOS_OBJ_LIB) $(INCLUDE_DIR) $(SGTE_LIB_DYNAMIC)
	@mkdir -p $(LIB_DIR)
	$(COMPILE) -shared -o $(ALGOS_LIB_DYNAMIC) $(ALGOS_OBJ_LIB) $(CXXFLAGS_ALGOS_LIBS) $(LIB_DYNAMIC_UTILS) $(LIB_DYNAMIC_EVAL) $(LIB_DYNAMIC_SGTE)
# Need to install the correct name for libsgtelib.so because it is moved from its original position
ifeq ($(UNAME), Darwin)
ifeq ($(USE_SGTELIB), 1)
	install_name_tool -change $(SGTE_LIB_DYNAMIC_NAME) $(LIB_DIR)/$(SGTE_LIB_DYNAMIC_NAME) $(ALGOS_LIB_DYNAMIC)
endif
endif

$(SGTE_LIB_DYNAMIC): $(SGTE_OBJ) $(SGTE_LIB_DIR)/$(SGTE_LIB_DYNAMIC_NAME)
	@mkdir -p $(LIB_DIR)
	@rm -f $(LIB_DIR)/$(SGTE_LIB_DYNAMIC_NAME)
	@ln -s $(SGTE_LIB_DIR)/$(SGTE_LIB_DYNAMIC_NAME) $(LIB_DIR)/$(SGTE_LIB_DYNAMIC_NAME)


clean:
	@cd $(ATTR_COMPONENT) && $(MAKE) clean
	@cd $(SGTE_SRC_DIR) && $(MAKE) clean
	@rm -rf $(BUILD_DIR)
