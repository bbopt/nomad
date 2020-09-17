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

COMPONENT_DIRNAME   = Attribute
ALL_FILES           = cacheAttributesDefinition displayAttributesDefinition \
                      evalAttributesDefinition evaluatorControlAttributesDefinition \
                      evaluatorControlGlobalAttributesDefinition \
                      pbAttributesDefinition runAttributesDefinition

ALL_HEADERS         = $(addsuffix .hpp, $(ALL_FILES))
ALL_TXT             = $(addsuffix .txt, $(ALL_FILES))
ALL_INCLUDE_HEADERS = $(addprefix $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/,$(ALL_HEADERS))

EXE = $(BIN_DIR)/WriteAttributeDefinitionFile
OBJ = $(OBJ_DIR)/WriteAttributeDefinitionFile.o


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

OBJFLAGS            = -c

COMPILE             = g++ $(CXXFLAGS)


all: $(EXE)

$(OBJ): WriteAttributeDefinitionFile.cpp $(ALL_TXT)
	@mkdir -p $(OBJ_DIR)
	$(COMPILE) $(OBJFLAGS) $< -o $@

$(EXE) : $(OBJ) $(ALL_TXT)
	@mkdir -p $(BIN_DIR)
	$(COMPILE) $(OBJ) -o $@
	$(EXE)
	@mkdir -p $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)
	@cp -f $(ALL_HEADERS)  $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/


clean:
	@rm -f $(EXE) $(OBJ)
	@rm -f *.hpp $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/*.hpp
