SRCS=$(filter-out vibes.cpp, $(wildcard *.cpp))
BINS=$(SRCS:.cpp=)

CXXFLAGS := $(shell pkg-config --cflags ibex) 
LIBS	 := $(shell pkg-config --libs  ibex)

ifeq ($(DEBUG), yes)
CXXFLAGS := $(CXXFLAGS) -O0 -g -pg -Wall 
else
CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG 
endif

CXXFLAGS := $(CXXFLAGS) -std=c++11 -DIBEX_BENCHS_DIR=\"${ROOT_DIR}../benchs/solver\" -U__STRICT_ANSI__

DOC_OUT_TXT=doc-arithmetic.txt doc-modeling.txt doc-separator.txt doc-set.txt 
DOC_SET_OUT=set-sep set-example set-inter set-union set-contract set-interval

EX_WITH_VIBES="doc-set lab1 lab2 lab3 lab4 lab5 lab6 lab7 lab8"
VIBES_MSG="INFO: The example $@ uses VIBes to plot data, you need to launch VIBes-viewer before running this test"

all: $(BINS)

% :	%.cpp
	@case $(EX_WITH_VIBES) in *$@*) echo "$(VIBES_MSG)" ;; esac
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(BINS) $(DOC_OUT_TXT) $(DOC_SET_OUT)
	
#plugins/optim/examples/doc-contractor.cpp
