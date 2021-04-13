
# Verify Anaconda is loaded. Needed for Cython.Build.
conda_version=$(shell conda --version)
ifeq ("$(conda_version)","")
$(error Anaconda is needed to run this makefile)
else
$(info Anaconda Version $(conda_version))
endif

ifdef USE_OPENMP
FLAG_OPENMP=1
$(info USE_OPENMP is set)
else
$(info USE_OPENMP is not set)
FLAG_OPENMP=0
endif

TOP = $(shell pwd | sed 's/\/interfaces.*//')


all:
	python setup_PyNomad.py $(FLAG_OPENMP) $(TOP)/build/release build_ext --inplace

lib:
	python setup_PyNomad.py $(FLAG_OPENMP) $(TOP)/build/release install --user

clean:
	@echo Cleaning PyNomad directory...
	@rm -rf build *.so PyNomad.cpp
	@echo done.
