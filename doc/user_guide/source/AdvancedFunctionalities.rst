.. _advanced_functionalities:

Advanced functionalities
========================

...

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
A blackbox program may fail to evaluate some of the trial points. 

a failed evaluation. However, when a block of trial points is submitted the content of the output
file must specify which points have failed by using the keyword FAIL at the corresponding position
in the output file. The keyword FAIL should be put only once per trial point independently
of the number of outputs given by BB_OUTPUT_TYPE.

If one value provided in the output file
cannot be read by NOMAD then the corresponding trial point is also considered as having failed.
The trial points that have failed will not be evaluated again.
A blackbox program can stop prematurely the evaluation of a block of trial points, for example
when a best incumbent trial point has been identified within the blackbox. However, to prevent
that the remaining trial points of the block be considered as having failed, it is required to
explicitly reject them by putting the keyword REJECT instead of their corresponding objective
and constraint functions. Contrary to trial points that have failed, the rejected ones could be
resubmitted later on.
An example of blackbox program written in Perl scripting language is provided in the example
directory (see Fig. 71). This script calls up to 4 instances of the executable bb.exe to evaluate
4 trial points in parallel.

Figure 71: Example of a blackbox program evaluating a block of 4 trial points in $NOMAD_HOME/examples/basic/batch/single_obj_parallel.
The parameter file that specifies this blackbox program with blocks of 4 trial points is given in
Fig. 72.

Figure 72: Example of parameter file with BB_MAX_BLOCK_SIZE
$NOMAD_HOME/examples/basic/batch/single_obj_parallel.

Library mode
^^^^^^^^^^^^
Please refer to $NOMAD_HOME/examples/basic/library/single_obj_parallel for an ex-
ample on how to manage a block of evaluations in parallel using pThreads and Semaphore.

When evaluations are performed by blocks (EVAL_LIST_MAX_BLOCK_SIZE greater than one) the
opportunistic strategy applies after evaluating a block of trial points.



Parallel evaluations
--------------------

...

PSD Mads
--------

...
