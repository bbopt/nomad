
ifndef ($(VARIANT))
VARIANT             = release
endif

ifdef NOOMP
$(info Not using OpenMP)
endif

UNAME := $(shell uname)

TOP                 = $(abspath ../..)
BUILD_DIR           = $(TOP)/build/$(VARIANT)
SRC_DIR             = $(TOP)/src
OBJ_DIR             = $(BUILD_DIR)/obj
INCLUDE_DIR         = $(BUILD_DIR)/include
LIB_DIR             = $(BUILD_DIR)/lib
BIN_DIR             = $(BUILD_DIR)/bin
EXE                 = $(BIN_DIR)/nomad

COMPONENT_DIRNAME   = Nomad

ALL_FILES           = nomad
ALL_OBJ             = $(addsuffix .o, $(addprefix $(OBJ_DIR)/,$(ALL_FILES)))
ALL_HEADERS         = $(addsuffix .hpp, $(ALL_FILES))
ALL_INCLUDE_HEADERS = $(addprefix $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/,$(ALL_HEADERS))


NOMAD_SRC           = nomad.cpp
LIB_DYNAMIC         =
SGTE_LIB            =

CXXFLAGS            = -L$(LIB_DIR)
ifeq ($(VARIANT), release)
CXXFLAGS           += -O2
else
CXXFLAGS           += -g -DDEBUG
endif
ifdef USE_UNORDEREDSET
CXXFLAGS            += -DUSE_UNORDEREDSET
endif
ifdef USE_PRIORITY_QUEUE
CXXFLAGS            += -DUSE_PRIORITY_QUEUE
endif

CXXFLAGS            += -std=c++14 -Wall -fpic
# Use OpenMP for parallelism (threads)
ifndef NOOMP
CXXFLAGS            += -fopenmp
endif

ifeq ($(UNAME), Linux)
CXXFLAGS_LIBS = -Wl,-rpath,$(LIB_DIR)
endif


OBJFLAGS            = -c

INCLFLAGS           = -I$(INCLUDE_DIR)

COMPILE             = g++ $(CXXFLAGS) $(INCLFLAGS)


all: $(ALL_INCLUDE_HEADERS) $(EXE)

$(EXE): $(NOMAD_HEADER) $(NOMAD_SRC)
	@mkdir -p $(BIN_DIR)
	$(COMPILE) $(NOMAD_SRC) -o $@ $(CXXFLAGS_LIBS) $(LIB_DYNAMIC_LINK)

$(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/%.hpp: %.hpp
	@mkdir -p $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)
	@cp -f $<  $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/$<

clean:
	@rm -f $(INCLUDE_DIR)/$(COMPONENT_DIRNAME)/*.hpp
	@rm -f $(ALL_OBJ)
