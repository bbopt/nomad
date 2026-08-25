//////////// THIS FILE MUST BE CREATED BY EXECUTING WriteAttributeDefinitionFile ////////////
//////////// DO NOT MODIFY THIS FILE MANUALLY ///////////////////////////////////////////////

#ifndef __NOMAD_4_6_RUNATTRIBUTESDEFINITIONCOOP__
#define __NOMAD_4_6_RUNATTRIBUTESDEFINITIONCOOP__

_definition = {
{ "COOP_MADS_OPTIMIZATION",  "bool",  "false",  " COOP-MADS optimization algorithm ",  " \n  \n . Use COOP-MADS algorithm. Only available if code is built with OpenMP enabled. \n  \n . Argument: bool \n  \n . Description: Parallel concurrent Mads \n  \n . This option deactivates any other optimization strategy. \n  \n . Example: COOP_MADS_OPTIMIZATION true \n  \n . Default: false\n\n",  "  advanced coop mads parallel  "  , "true" , "false" , "true" },
{ "COOP_MADS_NB_PROBLEM",  "size_t",  "4",  " Number of COOP-MADS problems ",  " \n  \n . When using COOP MADS optimization, select the number of \n   Mads problems ran in parallel. \n    \n . In addition each Mads algorithm can perform evaluations in parallel when \n   the NB_THREADS_PARALLEL_EVAL is greater than 1. \n  \n . Argument: a positive integer. \n  \n . This attribute is used only when COOP-Mads optimization is active. \n  \n . Example: COOP_MADS_NB_PROBLEM 2 \n  \n . Default: 4\n\n",  "  advanced psd mads parallel decomposition subproblem  "  , "true" , "false" , "true" },
{ "COOP_MADS_OPTIMIZATION_CACHE_SEARCH",  "bool",  "false",  " COOP-MADS cache search for incumbent synchronization ",  " \n  \n . Perform a cache search to update the best incumbent obtained by all Mads. \n  This allows to synchronize the best solutions of the parallel Mads instances. \n  \n . This attribute is used only when COOP-Mads optimization is active. \n  \n . Argument: bool. \n  \n . Example: COOP_MADS_OPTIMIZATION_CACHE_SEARCH false \n  \n . Default: false\n\n",  "  advanced psd mads parallel decomposition subproblem  "  , "true" , "false" , "true" } };

#endif
