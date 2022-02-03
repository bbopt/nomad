Creating the Python Interface to NOMAD (PyNomad) is based on Cython 0.24
(or above) and Python 3.7.
A simple way to have Cython and Python is to install the Anaconda package.

KNOWN ISSUES
If NOMAD has been compiled with gcc 9.1, running PyNomad may fail.
In that case, recompiling NOMAD with gcc 8.2 should fix the issue.
Building is disabled when OpenMP is enabled because of Python GIL
problems when using threads in Nomad library.  

Not all functionalities of Nomad are available in PyNomad.

HOW TO BUILD
The interface build is managed by CMake. See README.txt at NOMAD root. 
The next release will provide installation through pip using a wheel.


HOW TO USE
Some tests are proposed in the directory to check that everything is up and running.

Import PyNomad as a module and run PyNomad.info() to obtain the interface usage.
To obtain help on a Nomad parameter, run PyNomad.help("keyword").

Basic tests can be performed:
        python runTest.py
        python runTest_BlockEval.py
        Python runTest_PbWithConst.py
