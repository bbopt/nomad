Creating the Python Interface to NOMAD (PyNomad) requires to build 
source codes. The build procedure relies on Python 3.6 and Cython 0.24 
or higher. A simple way to make it work is to first install the Anaconda 
package manager.

************
KNOWN ISSUES
************
If NOMAD has been compiled with gcc 9.1, running PyNomad may fail.

************************
HOW TO BUILD AND INSTALL
************************
The interface build is managed by CMake that can be run at NOMAD root. 

The configuration command:
   cmake -DBUILD_INTERFACE_PYTHON=ON -S . -B build/release 
must be performed with Cython available (that can be done within a Conda 
environment: conda activate ... or activate ...).

For Windows, the default Anaconda is Win64. Visual Studio can support both 
Win32 and Win64 compilations. The configuration must be forced to use Win64:
   cmake -DBUILD_INTERFACE_PYTHON=ON -S . -B build/release -G"Visual Studio 15 2017 Win64". 
The Visual Studio version must be adapted.

The command 
   cmake --build build/release 
or cmake --build build/release --config Release (for Windows)
is used for building the selected configuration.

The command 
   cmake --install build/release
must be run before using the PyNomad module.


The next release of NOMAD will provide installation through pip using 
a wheel.

**********
HOW TO USE
**********
Some tests are proposed in the directory to check that everything 
is up and running.

Import PyNomad as a module and run PyNomad.info() to obtain the interface 
usage. To obtain help on a Nomad parameter, run PyNomad.help("keyword").

NOMAD parameters are provided in a list of strings using the same syntax 
as used in the NOMAD parameter files. Several tests and examples are 
proposed in the PyNomad directory to check that everything is up and 
running:
        python runTest.py
        python runTest_BlockEval.py
        Python runTest_PbWithConst.py

