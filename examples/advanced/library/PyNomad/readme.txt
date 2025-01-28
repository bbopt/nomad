# Project Title

Simple examples to illustrate PyNomad utilization.

## Installation and limitations

To run these examples, PyNomad must have been installed. Refer to the user guide or $NOMAD_HOME/README.txt for building PyNomad from source. PyNomad is also available from PyPi (pip install PyNomadBBO) for different OS and Python versions.

PyNomad built for OSX using Clang does not support OpenMP. Example 'simpleExample_basic_parallelEval.py' will trigger an exception because the parameter 'NB_THREADS_PARALLEL_EVAL' is set to a value greater than one.
This concern also the versions for OSX obtained from PyPi. 

For OSX, it is possible to build PyNomad using gcc which supports OpenMP.


