# Utilization of a cache file for a hot restart

This project illustrates how a first optimization run can be continued with 
a second run with more evaluations.

The example uses a cache file to resume an optimization if the initial attempt
is prematurely stopped. To simulate this situation we have two parameter files
that can be ran in sequence.

This functionality is only possible if the problem definition is not changed
between the runs.

## Usage

To run the first optimization execute the following command:

$NOMAD_HOME/bin/nomad param_firstRun.txt

This creates a cache file `cache.txt` containing 100 evaluated points.

To continue this optimization with 50 more points, execute the following command:

$NOMAD_HOME/bin/nomad param_secondRun.txt

In the second run, the `MAX_BB_EVAL` parameter is increased to 150, and the 
`USE_CACHE_FILE_FOR_RERUN` option is enabled. This means that the first 100 
evaluations of the second run will be loaded from the cache file instead of 
being re-evaluated by the blackbox function. The optimization algorithm will 
then proceed with new evaluations until it reaches the maximum number of 
evaluations specified.



## Notes

Running two times using `param_firstRun.txt` will produce different results
because the cache file is used. We can note that the cache hits increase between
two runs. 

For using the best point obtained after a first run and stored in a cache,
the user can remove initial point provide in `X0`.
 