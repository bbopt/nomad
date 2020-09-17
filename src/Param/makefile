
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

COMPONENT_DIRNAME   = Param

ALL_FILES           = AllParameters AttributeFactory Attribute CacheParameters \
                      DisplayParameters EvalParameters EvaluatorControlParameters \
                      EvaluatorControlGlobalParameters \
                      ParameterEntries ParameterEntry Parameters ParametersNomad3 \
                      PbParameters RunParameters TypeAttribute
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


OBJFLAGS            = -c

COMPILE             = g++ $(CXXFLAGS)


all: $(ALL_INCLUDE_HEADERS) $(ALL_OBJ)

$(OBJ_DIR)/%.o: %.cpp %.hpp
	$(COMPILE) $(OBJFLAGS) $< -o $@

$(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/%.hpp: %.hpp
	@mkdir -p $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)
	@cp -f $<  $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/$<

clean:
	@rm -f $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/*.hpp
	@rm -f $(ALL_OBJ)
