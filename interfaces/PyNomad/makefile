
# Verify Anaconda is loaded. Needed for Cython.Build.
conda_version=$(shell conda --version)
ifeq ("$(conda_version)","")
$(error Anaconda is needed to run this makefile)
else
$(info Anaconda Version $(conda_version))
endif

ifdef USE_OPENMP
$(info USE_OPENMP is set)
else
$(info USE_OPENMP is not set)
endif

all:
	python setup_PyNomad.py $(USE_OPENMP) build_ext --inplace

lib:
	python setup_PyNomad.py install --user

clean:
	@echo Cleaning PyNomad directory...
	@rm -rf build *.so PyNomad.cpp
	@echo done.
