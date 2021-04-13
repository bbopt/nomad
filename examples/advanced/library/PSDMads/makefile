ifndef ($(VARIANT))
VARIANT             = release
endif

UNAME := $(shell uname)

TOP                 = $(shell pwd | sed 's/\/examples.*//')
BUILD_DIR           = $(TOP)/build/$(VARIANT)
OBJ_DIR             = $(BUILD_DIR)/obj
INCLUDE_DIR         = $(BUILD_DIR)/include/nomad/
LIB_DIR             = $(BUILD_DIR)/lib
BIN_DIR             = $(BUILD_DIR)/bin
EXE                 = rosenbrock.exe



UTILS_LIB_NAME    = nomadUtils
EVAL_LIB_NAME     = nomadEval
ALGOS_LIB_NAME    = nomadAlgos

SGTELIB_NAME              = sgtelib

LIB_DYNAMIC               = -l$(UTILS_LIB_NAME) -l$(EVAL_LIB_NAME) -l$(ALGOS_LIB_NAME) -l$(SGTELIB_NAME)

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

all: rosenbrock.exe

rosenbrock.exe: $(INCLUDE_DIR) rosenbrock.cpp
	$(COMPILE) -o $@ rosenbrock.cpp $(LIB_DYNAMIC)

clean:
	rm -f *.exe *.o

