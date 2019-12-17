EXE                = ../bin/sgtelib.exe
LIB_NAME           = libsgtelib.so
LIB_DIR            = ../lib
LIB                = $(LIB_DIR)/$(LIB_NAME)
COMPILATOR         = g++
COMPILATOR_OPTIONS = -O3 -pedantic -Wall -Wextra -fpic

LIB_CURRENT_VERSION= 1.0.0

LIBS               = -lc -lm 
LDLIB_EXE          = -lsgtelib 
LDFLAGS            = -L$(LIB_DIR)

UNAME             := $(shell uname)

BUILD_DIR          = ../build

ifeq ($(UNAME), Darwin)
SHARED_OPTIONS_LIB = -current_version $(LIB_CURRENT_VERSION) -compatibility_version $(LIB_CURRENT_VERSION) -install_name $(LIB_NAME) -headerpad_max_install_names
endif

ifeq ($(UNAME), Linux)
SHARED_OPTIONS_EXE = -Wl,-rpath-link,'.' -Wl,-rpath,'$$ORIGIN/../lib'
COMPILATOR_OPTIONS+= -ansi
LIBS              += -ldl
endif

INCLUDE            = -I.
COMPILE            = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS_LIB           = TrainingSet.o Surrogate_Parameters.o Surrogate_KS.o Surrogate_RBF.o \
                     Surrogate_PRS.o Surrogate_PRS_EDGE.o Surrogate_LOWESS.o Surrogate_Kriging.o\
                     Surrogate_PRS_CAT.o Surrogate_Ensemble.o Surrogate_CN.o \
                     Surrogate.o Matrix.o Kernel.o Surrogate_Utils.o Surrogate_Factory.o \
                     Tests.o sgtelib_help.o Metrics.o
OBJS_MAIN          = sgtelib.o 

OBJS_MAIN         := $(addprefix $(BUILD_DIR)/,$(OBJS_MAIN))
OBJS_LIB          := $(addprefix $(BUILD_DIR)/,$(OBJS_LIB))

.SECONDARY: %.o

all: $(OBJS_LIB) $(LIB) $(OBJS_MAIN) $(EXE)

$(EXE): $(OBJS_MAIN) $(LIB)
	@echo
	@echo "== Making $@ =="
	$(call DisplayCommand,$(COMPILATOR) -o $@ $< $(LIBS) $(LDLIB_EXE) $(COMPILATOR_OPTIONS) $(LDFLAGS) $(SHARED_OPTIONS_EXE))
ifeq ($(UNAME), Darwin)
	@install_name_tool -change $(LIB_NAME) @loader_path/$(LIB) $@
endif
	@echo

$(LIB): $(OBJS_LIB)
	@echo
	@echo "== Making $@ =="
	$(call DisplayCommand,$(COMPILATOR) -shared -o $@ $(OBJS_LIB) $(LIBS) $(COMPILATOR_OPTIONS) $(SHARED_OPTIONS_LIB))

$(BUILD_DIR)/%.o: %.cpp %.hpp
	@mkdir -p $(BUILD_DIR) $(LIB_DIR)	../bin
	$(call DisplayCommand,$(COMPILE) $< -o $@)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(BUILD_DIR) $(LIB_DIR) ../bin
	$(call DisplayCommand,$(COMPILE) $< -o $@)

.PHONY : clean	
clean:
	@echo " --> cleaning obj files" 
	$(call DisplayCommand0,rm -f $(OBJS_LIB) $(OBJS_MAIN) )

.PHONY : del
del:
	@echo " --> cleaning trash files" 
	$(call DisplayCommand0,rm -f core *~)
	@echo " --> cleaning obj files"
	$(call DisplayCommand0,rm -f $(OBJS_LIB) $(OBJS_MAIN) )
	@echo " --> cleaning exe file" 
	$(call DisplayCommand0,rm -f $(EXE))
	@echo " --> cleaning library file" 
	$(call DisplayCommand0,rm -f $(LIB))


ifdef SHORT_DISPLAY
define DisplayCommand	
	@if $(1) ; then echo "---> $@ OK" ; else echo "compilation failed for $@" ; fi
endef
define DisplayCommand0
	@if $(1) ; then echo " ----> OK" ; else echo " ----> FAILED" ; fi
endef
else
define DisplayCommand	
	$(1)
endef
define DisplayCommand0
	$(1)
endef
endif
