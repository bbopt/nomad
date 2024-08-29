**********************
**********************
NOMAD and Python 
**********************
**********************

NOMAD and Python can be used in combination. There are two ways to perform
optimization using objective and constraint function evaluated by a Python script.

The simplest way is to run the Python script as a blackbox in batch mode to evaluate each 
given point. An example is provided in $NOMAD_HOME/examples/basic/batch/PythonBB.

Another way is to obtain the Nomad interface for Python (PyNomad) as described in what 
follows. 


**********************
**********************
Installation from PyPi
**********************
**********************

The simplest way to install PyNomad package is to proceed from PyPI.

A python package installation guide is available at
https://packaging.python.org/en/latest/tutorials/installing-packages/

To install the last version, from a shell command line perform:

pip install PyNomadBBO 

PyNomadBBO from PyPI relies on Python3 version 3.8 and above.

**********
HOW TO USE
**********
Some tests are proposed in the directory to check that everything 
is up and running. From the command line, if pytest is available, 
simply run the command:

pytest

To have more info, start python in a shell, import PyNomad as a module 
and run PyNomad.info() to obtain the interface usage. To obtain help 
on a Nomad parameter, run PyNomad.help("keyword"). To list all 
attributes and functions, execute dir(PyNomad).

NOMAD parameters are provided in a list of strings using the same syntax 
as used in the NOMAD parameter files. Several tests and examples are 
proposed in the $NOMAD_HOME/examples/advanced/library/PyNomad directory.

        python3 simpleExample_basic.py
        python3 simpleExample_BlockEval.py
        Python3 simpleExample_PbWithConst.py



****************
****************
Building PyNomad 
****************
****************
Alternatively, one can build PyNomad (and Nomad binaries) from Nomad 
source code as described below.  

The build procedure relies on Python 3.6 and Cython 0.24 
or higher. 

The procedure is often straightforward but several issues that can prevent
a successful building of the binaries can arise.

************
KNOWN ISSUES
************
If NOMAD has been compiled with gcc 9.1, running PyNomad may fail.

On OSX with ARM64 and X86 architectures.
There might be incompatibility between the Python API binaries and
the binaries obtained when building Nomad. When importing PyNomad in
a python shell, this can result in a message like :
ImportError: dlopen(xxxx/PyNomad.cpython-39-darwin.so, 0x0002) .... 
(mach-o file, but is an incompatible architecture (have 'arm64', need 'x86_64'))
or
ImportError: dlopen(xxxx/PyNomad.cpython-39-darwin.so, 0x0002): symbol not found 
in flat namespace '__ZN9NOMAD_4_210Parameters17_typeOfAttributesE'

This issue can be resolved by forcing the architecture when configuring 
the build with a flag like -DCMAKE_OSX_ARCHITECTURES=x86_64. The
build directory and the PyNomad.cpython-39-darwin.so in 
$NOMAD_HOME/interfaces/PyNomad MUST BE REMOVED BEFORE REBUILD.

************************
HOW TO BUILD AND INSTALL
************************
The interface build is managed by CMake that can be run at NOMAD root. 
For now, PyNomad cannot work with OpenMP enabled.

The configuration command:
   cmake -DBUILD_INTERFACE_PYTHON=ON -DTEST_OPENMP=OFF -S . -B build/release 
must be performed with Cython available (that can be done within a Conda 
environment: conda activate ... or activate ... OR with a virtual environment
containing cython and wheel).  

For Windows, the default Anaconda is Win64. Older Visual Studio versions can 
support both Win32 and Win64 compilations. The configuration must be forced
 to use Win64:
   cmake -DBUILD_INTERFACE_PYTHON=ON -DTEST_OPENMP=OFF -S . -B build/release -G"Visual Studio 15 2017 Win64". 

The Visual Studio version must be adapted.

The command 
   cmake --build build/release 
or cmake --build build/release --config Release (for Windows)
is used for building the selected configuration.

The command 
   cmake --install build/release
must also be run to install libraries.

IMPORTANT:
To install PyNomad wheel in your Python environment you must do
    pip install --user --force-reinstall dist/*whl 
Or
    pip install dist/*whl 
In the the PyNomad directory

**********
HOW TO USE
**********
Some tests are proposed in the directory to check that everything 
is up and running. From the command line, if pytest is available, 
simply run the command:

pytest

To have more info, start python in a shell, import PyNomad as a module 
and run PyNomad.info() to obtain the interface usage. To obtain help 
on a Nomad parameter, run PyNomad.help("keyword"). To list all 
attributes and functions, execute dir(PyNomad).

NOMAD parameters are provided in a list of strings using the same syntax 
as used in the NOMAD parameter files. Several tests and examples are 
proposed in the $NOMAD_HOME/examples/advanced/library/PyNomad directory.

        python3 simpleExample_basic.py
        python3 simpleExample_BlockEval.py
        Python3 simpleExample_PbWithConst.py
