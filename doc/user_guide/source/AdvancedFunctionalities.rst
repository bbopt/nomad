.. _advanced_functionalities:

Advanced functionalities
========================

Advanced parameters
-------------------




.. _bloc_evaluations

Blackbox evaluation of a block of trial points
----------------------------------------------

At different phases of the MADS algorithm, different numbers of trial points are generated.
For example, having selected the direction type as ORTHO 2N, the maximum number of points generated during
the Poll step will be 2N+2. These points can be partitioned into blocks of trial points to be
submitted sequentially for evaluation to a blackbox program. The maximum size of a block of
evaluations is controlled by the BB_MAX_BLOCK_SIZE. By default, a block contains a single trial
point. This can be changed by the user but the blackbox program must support the evaluation
of a varying number of trial points, up to BB_MAX_BLOCK_SIZE.

Due to the strategy of by-block evaluation, the maximum number of evaluations requested to
NOMAD may be exceeded if BB_MAX_BLOCK_SIZE > 1. The reason for this behaviour is that
block results are analyzed only after completion and the maximum number of evaluations may
be exceeded when checking this termination criterion.
The opportunistic strategy (enabled by default) may apply after each block of trial points.
Evaluations of blocks of trial points can be performed in parallel by the blackbox program. This
strategy of parallelization must be setup by the user within the blackbox. Examples are provided
in what follows.

Batch mode
^^^^^^^^^^
In batch mode, NOMAD creates input files which can contain at most
BB_MAX_BLOCK_SIZE trial points separated by a linebreak. Each point is given as a row of values.
The user must provide a blackbox program that can read the input file, evaluate them and
output the objective and constraints functions (in the order provided by the BB_OUTPUT_TYPE
parameter) for each trial point in the same order as provided in the input file.
A blackbox program may fail to evaluate some of the trial points. When block of trial points is
submitted the content of the output file must reflect the outputs for each point.
If one value provided in the output file
cannot be read by NOMAD, then the corresponding trial point is considered as having failed.
The trial points that have failed will not be evaluated again.
An example of blackbox program written in Perl scripting language is provided in the
directory ``$NOMAD_HOME/examples/basic/batch/single_obj_parallel``. The script ``parallel_BBWrapper.pl``
calls up to 4 instances of the executable bb.exe to evaluate 4 trial points in parallel.

::

  > cd $NOMAD_HOME/examples/basic/batch/single_obj_parallel
  > more x.txt
  1 2 3 4 5
  0 0 0 0 0
  2 2 2 2 2
  5 4 3 2 1
  > perl parallel_BBWrapper.pl x.txt
  5 5 -65
   0 -20 20
   2 -20 -20
   1 5 -65

The same directory holds the parameter file that specifies this blackbox program with blocks of 4 trial points:

::

    DIMENSION      5              # number of variables

    BB_EXE "$perl parallel_BBWrapper.pl"
    BB_MAX_BLOCK_SIZE 4

    BB_OUTPUT_TYPE OBJ PB EB

    X0             ( 0 0 0 0 0 )  # starting point

    LOWER_BOUND    * -6.0         # all variables are >= -6
    UPPER_BOUND    ( 5 6 7 - - )  # x_1 <= 5, x_2 <= 6, x_3 <= 7
                                  # x_4 and x_5 have no bounds

    MAX_BLOCK_EVAL     20         # the algorithm terminates when
                                  # 20 blocks have been evaluated

    TMP_DIR /tmp
    DISPLAY_DEGREE 2
    DISPLAY_STATS BLK_EVA BLK_SIZE OBJ
    DISPLAY_ALL_EVAL true


Library mode
^^^^^^^^^^^^
Please refer to ``$NOMAD_HOME/examples/basic/library/single_obj_parallel`` for an example
on how to manage a block of evaluations in parallel using pThreads and Semaphore.

When evaluations are performed by blocks, i.e., when ``BB_MAX_BLOCK_SIZE`` is greater
than one, the opportunistic strategy applies after evaluating a block of trial points.



.. _parallel_evaluations

Parallel evaluations
--------------------

...

PSD Mads
--------

...
