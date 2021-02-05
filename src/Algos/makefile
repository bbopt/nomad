
ifndef ($(VARIANT))
VARIANT             = release
endif

# Note: if NOOMP is not set, flag -fopenmp will be used, and that defines _OPENMP.
ifdef NOOMP
$(info Not using OpenMP)
endif
ifneq ($(USE_SGTELIB), 0)
    USE_SGTELIB = 1
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
PSDMADS_DIRNAME     = $(COMPONENT_DIRNAME)/PSDMads
QUADMODEL_DIRNAME   = $(COMPONENT_DIRNAME)/QuadModel
SGTELIBMODEL_DIRNAME = $(COMPONENT_DIRNAME)/SgtelibModel
SSDMADS_DIRNAME     = $(COMPONENT_DIRNAME)/SSDMads
SGTEENSEMBLE_DIRNAME = $(COMPONENT_DIRNAME)/SgteEnsembleRenaud # added by renaud

LH_COMPONENT        = $(SRC_DIR)/$(LH_DIRNAME)
MADS_COMPONENT      = $(SRC_DIR)/$(MADS_DIRNAME)
NM_COMPONENT        = $(SRC_DIR)/$(NM_DIRNAME)
PHASEONE_COMPONENT  = $(SRC_DIR)/$(PHASEONE_DIRNAME)
PSDMADS_COMPONENT   = $(SRC_DIR)/$(PSDMADS_DIRNAME)
QUADMODEL_COMPONENT = $(SRC_DIR)/$(QUADMODEL_DIRNAME)
SGTELIBMODEL_COMPONENT = $(SRC_DIR)/$(SGTELIBMODEL_DIRNAME)
SSDMADS_COMPONENT   = $(SRC_DIR)/$(SSDMADS_DIRNAME)
SGTEENSEMBLE_COMPONENT = $(SRC_DIR)/$(SGTEENSEMBLE_DIRNAME) # added by renaud

ALL_FILES           = Algorithm AlgoStopReasons CacheInterface \
                      EvcInterface Initialization Iteration IterationUtils \
                      MainStep MegaIteration MeshBase Projection \
                      Step Subproblem SubproblemManager Termination
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
ifeq ($(USE_SGTELIB), 1)
CXXFLAGS            += -DUSE_SGTELIB
endif
ifeq ($(TIME_STATS), 1)
CXXFLAGS            += -DTIME_STATS
endif


OBJFLAGS            = -c
INCLFLAGS           = -I$(INCLUDE_DIR)

COMPILE             = g++ $(CXXFLAGS)


all: $(ALL_INCLUDE_HEADERS) $(ALL_OBJ)
	@cd $(LH_COMPONENT) && $(MAKE) all
	@cd $(MADS_COMPONENT) && $(MAKE) all
	@cd $(NM_COMPONENT) && $(MAKE) all
	@cd $(PHASEONE_COMPONENT) && $(MAKE) all
ifndef NOOMP
	@cd $(PSDMADS_COMPONENT) && $(MAKE) all
endif
ifeq ($(USE_SGTELIB), 1)
	@cd $(QUADMODEL_COMPONENT) && $(MAKE) all
	@cd $(SGTELIBMODEL_COMPONENT) && $(MAKE) all
	@cd $(SGTEENSEMBLE_COMPONENT) && $(MAKE) all
endif  #  @cd $(SGTEENSEMBLE_COMPONENT) && $(MAKE) all  added by renaud
	@cd $(SSDMADS_COMPONENT) && $(MAKE) all
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
	@cd $(PSDMADS_COMPONENT) && $(MAKE) clean
	@cd $(QUADMODEL_COMPONENT) && $(MAKE) clean
	@cd $(SGTELIBMODEL_COMPONENT) && $(MAKE) clean
	@cd $(SSDMADS_COMPONENT) && $(MAKE) clean
	@cd $(SGTEENSEMBLE_COMPONENT) && $(MAKE) clean

#  @cd $(SGTEENSEMBLE_COMPONENT) && $(MAKE) clean  added by renaud
