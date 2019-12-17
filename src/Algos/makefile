
ifndef ($(VARIANT))
VARIANT             = release
endif

# Note: if NOOMP is not set, flag -fopenmp will be used, and that defines _OPENMP.
ifdef NOOMP
$(info Not using OpenMP)
endif


UNAME := $(shell uname)

TOP                 = $(abspath ../..)
BUILD_DIR           = $(TOP)/build/$(VARIANT)
SRC_DIR             = $(TOP)/src
OBJ_DIR             = $(BUILD_DIR)/obj
INCLUDE_DIR         = $(BUILD_DIR)/include
BIN_DIR             = $(BUILD_DIR)/bin

COMPONENT_DIRNAME   = Algos
LH_DIRNAME          = $(COMPONENT_DIRNAME)/LatinHypercubeSampling
MADS_DIRNAME        = $(COMPONENT_DIRNAME)/Mads
NM_DIRNAME          = $(COMPONENT_DIRNAME)/NelderMead
PHASEONE_DIRNAME    = $(COMPONENT_DIRNAME)/PhaseOne
SGTELIBMODEL_DIRNAME = $(COMPONENT_DIRNAME)/SgtelibModel

LH_COMPONENT        = $(SRC_DIR)/$(LH_DIRNAME)
MADS_COMPONENT      = $(SRC_DIR)/$(MADS_DIRNAME)
NM_COMPONENT        = $(SRC_DIR)/$(NM_DIRNAME)
PHASEONE_COMPONENT  = $(SRC_DIR)/$(PHASEONE_DIRNAME)
SGTELIBMODEL_COMPONENT = $(SRC_DIR)/$(SGTELIBMODEL_DIRNAME)

ALL_FILES           = Algorithm AlgoStopReasons AllStopReasons CacheInterface \
                      EvcInterface Initialization Iteration IterationUtils \
                      MainStep MegaIteration MeshBase Projection \
                      Step Subproblem Termination
ALL_OBJ             = $(addsuffix .o, $(addprefix $(OBJ_DIR)/,$(ALL_FILES)))
ALL_HEADERS         = $(addsuffix .hpp, $(ALL_FILES))
ALL_INCLUDE_HEADERS = $(addprefix $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/,$(ALL_HEADERS))


ifeq ($(VARIANT), release)
CXXFLAGS            = -O2
else
CXXFLAGS            = -g -DDEBUG
endif
CXXFLAGS            += -std=c++14 -Wall -fpic
ifndef NOOMP
# Use OpenMP for parallelism (threads)
CXXFLAGS            += -fopenmp
endif
ifdef USE_PRIORITY_QUEUE
CXXFLAGS            += -DUSE_PRIORITY_QUEUE
endif
ifeq ($(USE_SGTELIB), 1)
CXXFLAGS            += -DUSE_SGTELIB
endif


OBJFLAGS            = -c
INCLFLAGS           = -I$(INCLUDE_DIR)

COMPILE             = g++ $(CXXFLAGS)


all: $(ALL_INCLUDE_HEADERS) $(ALL_OBJ)
	@cd $(LH_COMPONENT) && $(MAKE) all
	@cd $(MADS_COMPONENT) && $(MAKE) all
	@cd $(NM_COMPONENT) && $(MAKE) all
	@cd $(PHASEONE_COMPONENT) && $(MAKE) all
ifeq ($(USE_SGTELIB), 1)
	@cd $(SGTELIBMODEL_COMPONENT) && $(MAKE) all
endif
	@echo

$(OBJ_DIR)/%.o: %.cpp %.hpp
	$(COMPILE) $(OBJFLAGS) $(INCLFLAGS) $< -o $@


$(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/%.hpp: %.hpp
	@mkdir -p $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)
	@cp -f $<  $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/$<


clean:
	@rm -f $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/*.hpp
	@rm -f $(ALL_OBJ)
	@cd $(LH_COMPONENT) && $(MAKE) clean
	@cd $(MADS_COMPONENT) && $(MAKE) clean
	@cd $(NM_COMPONENT) && $(MAKE) clean
	@cd $(PHASEONE_COMPONENT) && $(MAKE) clean
	@cd $(SGTELIBMODEL_COMPONENT) && $(MAKE) clean
