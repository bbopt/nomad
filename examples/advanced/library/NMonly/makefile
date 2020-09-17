ifndef ($(VARIANT))
VARIANT             = release
endif

UNAME := $(shell uname)

TOP                 = $(shell pwd | sed 's/\/examples.*//')
BUILD_DIR           = $(TOP)/build/$(VARIANT)
OBJ_DIR             = $(BUILD_DIR)/obj
INCLUDE_DIR         = $(BUILD_DIR)/include
LIB_DIR             = $(BUILD_DIR)/lib
BIN_DIR             = $(BUILD_DIR)/bin
EXE                 = fixedVariable.exe


UTILS_LIB_CURRENT_VERSION = 4.0.0
EVAL_LIB_CURRENT_VERSION = 4.0.0
ALGOS_LIB_CURRENT_VERSION = 4.0.0

UTILS_NAME_AND_VERSION    = nomadUtils.$(UTILS_LIB_CURRENT_VERSION)
EVAL_NAME_AND_VERSION     = nomadEval.$(EVAL_LIB_CURRENT_VERSION)
ALGOS_NAME_AND_VERSION    = nomadAlgos.$(ALGOS_LIB_CURRENT_VERSION)

SGTE_LIB_NAME             = sgtelib

LIB_DYNAMIC               = -l$(UTILS_NAME_AND_VERSION) -l$(EVAL_NAME_AND_VERSION) -l$(ALGOS_NAME_AND_VERSION) -l$(SGTE_LIB_NAME)


ifeq ($(VARIANT), release)
CXXFLAGS            = -O2
else
CXXFLAGS            = -g -DDEBUG
endif

CXXFLAGS            += -std=c++14 -Wall -fpic
# Use OpenMP for parallelism (threads)
ifndef NOOMP
CXXFLAGS            += -fopenmp
endif

CXXFLAGS           += -L$(LIB_DIR)

ifeq ($(UNAME), Linux)
CXXFLAGS_LIBS = -Wl,-rpath,$(LIB_DIR)
endif

INCLFLAGS           = -I$(INCLUDE_DIR)

COMPILE             = $(CXX) $(CXXFLAGS) $(INCLFLAGS) $(CXXFLAGS_LIBS)

all: u.exe NMOpt.exe

NMOpt.exe: $(INCLUDE_DIR) $(OBJ_DIR) NMOpt.cpp
	$(COMPILE) -o $@ NMOpt.cpp $(LIB_DYNAMIC)


u.exe: u.cpp
	g++ -o u.exe u.cpp

clean:
	rm -f *.o *.exe cache.txt

