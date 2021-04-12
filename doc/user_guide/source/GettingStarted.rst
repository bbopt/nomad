.. _getting_started:

Getting started
===============

NOMAD is an efficient tool for simulation-based design optimizations provided in the form

...

How to create blackbox programs
===============================

To conduct optimization in batch mode the users must define their separate blackbox program coded as a stand-alone program. Blackbox program executions are managed by NOMAD with system calls.

How to provide parameters
=========================

In batch mode, the parameters are provided in a text file using predefined keywords followed by one or more argument. Here are some of the most important parameters defining an optimization problem (without brackets):

* The number of variables (DIMENSION n).
* The name of the blackbox executable that outputs the objective and the constraints (BB_EXE bb_name).
* Bounds on variables are defined with the LOWER_BOUND lb and UPPER_BOUND ub parameters.
* The output types of the blackbox executable: objective and constraints (BB_OUTPUT_TYPE obj cons1...consM).
* A starting point (X0 x0).
* An optional stopping criterion (MAX_BB_EVAL max_bb_eval, for example). If no stopping criterion is specified, the algorithm will stop as soon as the mesh size reaches a given tolerance.
* Any entry on a line is ignored after the character ‘#’.

The order in which the parameters appear in the file or their case is unimportant.

To continue

How to conduct optimization
===========================

Optimization is conducted by starting NOMAD executable in a command window with the parameter file name given as argument. To illustrate the execution, the example provided in $NOMAD_HOME/examples/basic/batch/single_obj/ is considered:

Todo
