Creating the Python Interface to NOMAD (PyNomad) is based on Cython 0.24
(or above) and Python 3.7.
A simple way to have Cython and Python is to install the Anaconda package.

KNOWN ISSUES
If NOMAD has been compiled with gcc 9.1, running PyNomad may fail.
In that case, recompiling NOMAD with gcc 8.2 should fix the issue.

Not all functionalities of Nomad are available in PyNomad.

HOW TO BUILD
The interface build is managed by CMake. See README.txt at NOMAD root.


HOW TO USE
Some tests are proposed in the directory to check that everything is up and running.
Run the following in a command line:
        python runTestInfoHelp.py
        python runTest.py
        python runTest_BlockEval.py
