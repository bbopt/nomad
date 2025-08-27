.. _getting_started:

Getting started
===============

NOMAD is an efficient tool for simulation-based design optimizations provided in the form:

.. math::

   \min_{x \in \Omega} f(x)

where the feasible set :math:`\Omega = \{ x \in X : c_j(x) \leq 0, j \in J\} \subset \mathbb{R}^n`, :math:`f, c_j : X \rightarrow \mathbb{R} \cup \{ \infty \}`
for  all :math:`j \in J= \{ 1,2,\ldots,m \}`, and where :math:`X` is a subset of :math:`\mathbb{R}^n`.
The functions :math:`f` and :math:`c_j`, :math:`j ∈ J`, are typically blackbox functions whose evaluations require computer simulation.

NOMAD can be used in two different modes: batch mode and library mode.
The batch mode is intended for a basic usage and is briefly presented below (more details will be provided in :ref:`basic_nomad_usage`).
This chapter explains how to get started with NOMAD in batch mode. The following topics will be covered:

* :ref:`create_blackbox_program`
* :ref:`provide_parameters` for defining the problem and displaying optimization results
* :ref:`conduct_optimization`

.. note::
   Building NOMAD binaries and running the examples provided during the installation requires to have a ``C++`` compiler installed on your machine.

   Compilation instructions rely on **CMake** and have been tested with **GCC** (GNU Compiler Collection), Clang and Visual Studio.


.. _create_blackbox_program:

Create blackbox programs
^^^^^^^^^^^^^^^^^^^^^^^^

To conduct optimization in batch mode the users must define their separate blackbox program coded as a standalone program.
Blackbox program executions are managed by NOMAD with system calls.

A valid blackbox program:
    `-` takes an input vector file as single argument,
    `-` reads space-separated values in input vector file,
    `-` returns evaluation values on standard output or file,
    `-` returns an evaluation status.

In what follows we use the example in the ``$NOMAD_HOME/examples/basic/batch/single_obj``.
This example optimization problem has a single objective, 5 variables, 2 nonlinear constraints and 8 bound constraints:


.. image:: ../figs/example1.png
    :align: left
    :width: 650
    :alt: example 1

.. note:: The blackbox programs may be coded in any language (even scripts) but must respect **NOMAD format**:

    1. The blackbox program must be callable in a terminal window at the command prompt and take the input vector file name as a single argument.
    For the example above, the blackbox executable is ``bb.exe``, one can execute it with the command  ``./bb.exe x.txt``. Here ``x.txt`` is a text file containing a total of 5 values.

    2. NOMAD will manage the creation of the **input file consisting of one value for each variable separated by space** and the execution of the blackbox program.

    3. The blackbox program must return the evaluation values by displaying them in the **standard output** (default) or by writing them in an output file (using ``BB_REDIRECTION no``, see example ``$NOMAD_HOME\examples\advanced\batch\BBOutputRedirection``).
    It must also **return an evaluation status of 0** to indicate that the evaluation went well. Otherwise NOMAD considers that the evaluation has failed.

    4. The minimum number of values displayed by the blackbox program corresponds to the number of constraints plus one (or two for bi-objective problems) representing the objective function(s) that one seeks to minimize.
    The constraints values correspond to left-hand side of constraints of the form :math:`c_j \leq 0` (for example, the constraint :math:`0 \leq x_1 + x_2 \leq 10` must be displayed with the two quantities :math:`c_1(x)=-x_1-x_2` and :math:`c_2(x)=x_1+x_2-10`).

The blackbox ``C++`` program of the previous example to evaluate the objective and the two constraints for a given design vector is given as:

.. code-block:: c++
   :linenos:

    #include <cmath>
    #include <iostream>
    #include <fstream>
    #include <cstdlib>
    using namespace std;

    int main ( int argc , char ** argv ) {

    double f = 1e20, c1 = 1e20 , c2 = 1e20;
    double x[5];

    if ( argc >= 2 ) {
        c1 = 0.0 , c2 = 0.0;
        ifstream in ( argv[1] );
        for ( int i = 0 ; i < 5 ; i++ ) {
            in >> x[i];
            c1 += pow ( x[i]-1 , 2 );
            c2 += pow ( x[i]+1 , 2 );
        }
        f = x[4];
        if ( in.fail() )
            f = c1 = c2 = 1e20;
        else {
            c1 = c1 - 25;
            c2 = 25 - c2;
        }
        in.close();
    }
    cout << f << " " << c1 << " " << c2 << endl;
    return 0;
    }

The blackbox compilation and test are as follows:

1. Change directory to ``$NOMAD_HOME/examples/basic/batch/single_obj``.

2. Optionally, compile the blackbox program with the following command ``g++ -o bb.exe bb.cpp`` (**GNU compiler**). This step is not really required because the building procedure with *CMake* normally builds the blackbox executable for this example.

3. Test the executable with the text file ``x.txt`` containing ``0 0 0 0 0`` by entering the command ``bb.exe x.txt``.

4. This test  should display ``0 -20 20``, which means that the point :math:`x = (0~0~0~0~0)^T` has an objective value of :math:`f(x)=0`,
but is not feasible, since the second constraint is not satisfied (:math:`c_2(x) = 20 > 0`).

::

  > cd $NOMAD_HOME/examples/basic/batch/single_obj
  > g++ -o bb.exe bb.cpp
  > more x.txt
  0 0 0 0 0
  > ./bb.exe x.txt
  0 -20 20

.. note::

  The order of the displayed outputs must correspond to the order defined in the parameter file (see :ref:`bb_output_type` for details).
  If variables have bound constraints, they must be defined in the parameters file and should not appear in the blackbox code.


.. _provide_parameters:

Provide parameters
^^^^^^^^^^^^^^^^^^

In batch mode, the parameters are provided in a text file using predefined keywords followed by one or more argument.

.. note::

  Help on parameters is accessible at the command prompt:
  ``$NOMAD_HOME/bin/nomad -h param_name`` (linux/osx) or ``%NOMAD_HOME%\build\release\bin\nomad.exe -h param_name``

Here are some of the most important parameters defining an optimization problem (without brackets):

* The number of variables (``DIMENSION n``).
* The name of the blackbox executable that outputs the objective and the constraints (``BB_EXE bb_name``).
* Bounds on variables are defined with the ``LOWER_BOUND lb`` and ``UPPER_BOUND ub`` parameters.
* The output types of the blackbox executable: objective and constraints (``BB_OUTPUT_TYPE obj cons1 ... consM``).
* A starting point (``X0 x0``).
* An optional stopping criterion (``MAX_BB_EVAL max_bb_eval``, for example).
  If no stopping criterion is specified, the algorithm will stop as soon as the mesh size reaches a given tolerance.
* Any entry on a line is ignored after the character ``#``.


.. note::

  The order in which the parameters appear in the file or their case is unimportant.

Example of a basic parameters file extracted from ``$NOMAD_HOME/examples/basic/batch/single_obj/param.txt``.
The comments in the file describes some of the syntactic rules to provide parameters:

::

    DIMENSION      5              # number of variables

    BB_EXE         bb.exe         # 'bb.exe' is a program that
    BB_OUTPUT_TYPE OBJ PB EB      # takes in argument the name of
                                  # a text file containing 5
                                  # values, and that displays 3
                                  # values that correspond to the
                                  # objective function value (OBJ),
                                  # and two constraints values g1
                                  # and g2 with g1 <= 0 and
                                  # g2 <= 0; 'PB' and 'EB'
                                  # correspond to constraints that
                                  # are treated by the Progressive
                                  # and Extreme Barrier approaches
                                  # (all constraint-handling
                                  #  options are described in the
                                  #  detailed parameters list)

    X0             ( 0 0 0 0 0 )  # starting point

    LOWER_BOUND    * -6           # all variables are >= -6
    UPPER_BOUND    ( 5 6 7 - - )  # x_1 <= 5, x_2 <= 6, x_3 <= 7
                                  # x_4 and x_5 have no bounds

    MAX_BB_EVAL    100            # the algorithm terminates when
                                  # 100 black-box evaluations have
                                  # been made



The constraints defined in the parameters file are of different types.
The first constraint :math:`c_1(x) \leq 0` is treated by the *Progressive Barrier* approach (*PB*), which allows constraint violations.
The second constraint, :math:`c_3(x) \leq 0`, is treated by the  *Extreme Barrier* approach (*EB*) that forbids violations.
Hence, evaluations not satisfying extreme barrier constraints are simply not considered when trying to improve the solution.

In the example above, the algorithmic parameters of NOMAD need not to be set because default
values are considered. This will provide the best results in most situations.


.. _conduct_optimization:

Conduct optimization
^^^^^^^^^^^^^^^^^^^^

Optimization is conducted by starting NOMAD executable in a command window with the parameter file name given as argument.

::

    $NOMAD_HOME/bin/nomad param.txt

To illustrate the execution, the example provided in ``$NOMAD_HOME/examples/basic/batch/single_obj/`` is considered::

  > cd $NOMAD_HOME/examples/basic/batch/single_obj
  > ls
  bb.cpp bb.exe CMakeLists.txt makefile param.txt x.txt
  >$NOMAD_HOME/bin/nomad param.txt
  BBE ( SOL ) OBJ
    1   (   0          0          0          0          0        )    0        (Phase One)
    8   (   0          4          0          0          0        )    0        (Phase One)
    28  (   1.4        5          0         -0.6       -0.4      )   -0.4
    29  (   2.6        4          0         -1.4       -0.8      )   -0.8
    33  (   1.63       3          0.92      -1.78      -0.88     )   -0.88
    37  (   2.46       3          0.97      -1.87      -0.92     )   -0.92
    41  (   3.2        3          0.16      -1.26      -1.05     )   -1.05
    42  (   4.27       2         -0.23      -1.07      -1.36     )   -1.36
    47  (   3.0        1          1.22      -1.92      -1.5      )   -1.5
    48  (   3.2        0          1.83      -2.19      -1.86     )   -1.86
    57  (   3.91      -0          1.02      -1.32      -1.95     )   -1.95
    67  (   3.61      -0          1.28      -1.83      -1.99     )   -1.99
    78  (   3.94       1          0.63      -1.14      -2.02     )   -2.02
    79  (   4.32       1          0.02      -0.61      -2.11     )   -2.11
    84  (   3.68       0          0.97      -1.23      -2.15     )   -2.15
    88  (   3.91       1          0.5       -0.6       -2.2      )   -2.2
    89  (   4.07       1          0.1        0.01      -2.31     )   -2.31
    94  (   3.67       1          0.56      -0.47      -2.36     )   -2.36
    95  (   3.35       1          0.84      -0.39      -2.48     )   -2.48
    99  (   4.15       1         -0.37       0.57      -2.49     )   -2.49
    Reached stop criterion: Max number of blackbox evaluations (Eval Global) 100
    A termination criterion is reached: Max number of blackbox evaluations (Eval Global) No more points to evaluate 100

    Best feasible solution:     #1540 ( 4.15 1 -0.37 0.57 -2.49 )   Evaluation OK    f =  -2.4900000000000002132     h =   0

    Best infeasible solution:   #1512 ( 3.79 0 1.14 -1.75 -1.97 )   Evaluation OK    f =  -1.9699999999999999734     h =   0.03500640999999999475

    Blackbox evaluations:        100
    Total model evaluations:     1348
    Cache hits:                  3
    Total number of evaluations: 103
