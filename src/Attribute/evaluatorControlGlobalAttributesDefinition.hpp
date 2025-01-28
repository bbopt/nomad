/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/
//////////// THIS FILE MUST BE CREATED BY EXECUTING WriteAttributeDefinitionFile ////////////
//////////// DO NOT MODIFY THIS FILE MANUALLY ///////////////////////////////////////////////

#ifndef __NOMAD_4_5_EVALUATORCONTROLGLOBALATTRIBUTESDEFINITION__
#define __NOMAD_4_5_EVALUATORCONTROLGLOBALATTRIBUTESDEFINITION__

_definition = {
{ "BB_MAX_BLOCK_SIZE",  "size_t",  "1",  " Size of blocks of points, to be used for parallel evaluations ",  " \n . Maximum size of a block of evaluations send to the blackbox \n   executable at once. Blackbox executable can manage parallel \n   evaluations on its own. Opportunistic strategies may apply after \n   each block of evaluations. \n  \n . Depending on the algorithm phase, the blackbox executable will \n   receive at most BB_MAX_BLOCK_SIZE points to evaluate. \n  \n . When this parameter is greater than one, the number of evaluations \n   may exceed the MAX_BB_EVAL stopping criterion. \n  \n . Argument: integer > 0. \n  \n . Example: BB_MAX_BLOCK_SIZE 3 \n            The blackbox executable receives blocks of \n            at most 3 points for evaluation. \n  \n . Default: 1\n\n",  "  advanced block parallel  "  , "true" , "true" , "true" },
{ "SURROGATE_MAX_BLOCK_SIZE",  "size_t",  "1",  " Size of blocks of points, to be used for parallel evaluations ",  " \n . Maximum size of a block of evaluations send to the surrogate \n   executable at once. Surrogate executable can manage parallel \n   evaluations on its own. \n  \n . Depending on the algorithm phase, the surrogate executable will \n   receive at most SURROGATE_MAX_BLOCK_SIZE points to evaluate. \n  \n . Argument: integer > 0. \n  \n . Example: SURROGATE_MAX_BLOCK_SIZE INF \n            The surrogate executable receives blocks with \n            all points evailable for evaluation. \n  \n . Default: 1\n\n",  "  advanced block parallel surrogate  "  , "true" , "true" , "true" },
{ "EVAL_QUEUE_CLEAR",  "bool",  "true",  " Opportunistic strategy: Flag to clear EvaluatorControl queue between each run ",  " \n  \n . Opportunistic strategy: If a success is found, clear evaluation queue of \n   other points. \n  \n . If this flag is false, the points in the evaluation queue that are not yet \n   evaluated might be evaluated later. \n  \n . If this flag is true, the points in the evaluation queue that are not yet \n   evaluated will be flushed. \n  \n . Outside of opportunistic strategy, this flag has no effect. \n  \n . Default: true\n\n",  "  advanced opportunistic oppor eval evals evaluation evaluations clear flush  "  , "true" , "true" , "true" },
{ "EVAL_SURROGATE_COST",  "size_t",  "INF",  " Cost of the surrogate function versus the true function ",  " \n   . Cost of the surrogate function relative to the true function \n  \n   . Argument: one nonnegative integer. \n  \n   . INF means there is no cost \n  \n   . Examples: \n         EVAL_SURROGATE_COST 3    # three surrogate evaluations count as one blackbox \n                                  # evaluation: the surrogate is three times faster \n         EVAL_SURROGATE_COST INF  # set to infinity: A surrogate evaluation does \n                                  # not count at all \n  \n   . See also: SURROGATE_EXE, EVAL_SURROGATE_OPTIMIZATION \n . Default: INF\n\n",  "  advanced static surrogate  "  , "true" , "false" , "true" },
{ "MAX_BB_EVAL",  "size_t",  "INF",  " Stopping criterion on the number of blackbox evaluations ",  " \n  \n . Maximum number of blackbox evaluations. When OpenMP is activated, this budget \n maybe exceeded due to parallel evaluations. \n  \n . Argument: one positive integer. \n  \n . An INF value serves to disable the stopping criterion. \n  \n . Does not consider evaluations taken in the cache (cache hits) \n  \n . Example: MAX_BB_EVAL 1000 \n  \n . Default: INF\n\n",  "  basic stop stops stopping max maximum criterion criterions blackbox blackboxes bb  "  , "false" , "true" , "true" },
{ "MAX_BLOCK_EVAL",  "size_t",  "INF",  " Stopping criterion on the number of blocks evaluations ",  " \n  \n . Maximum number of blocks evaluations \n  \n . Argument: one positive integer. \n  \n . An INF value serves to disable the stopping criterion. \n  \n . Example: MAX_BLOCK_EVAL 100 \n  \n . Default: INF\n\n",  "  advances block stop parallel  "  , "true" , "true" , "true" },
{ "MAX_EVAL",  "size_t",  "INF",  " Stopping criterion on the number of evaluations (blackbox and cache) ",  " \n  \n . Maximum number of evaluations, including evaluations taken in the cache \n   (cache hits) \n  \n . Argument: one positive integer. \n  \n . An INF value serves to disable the stopping criterion. \n  \n . Example: MAX_EVAL 1000 \n  \n . Default: INF\n\n",  "  advanced stop stops stopping max maximum criterion criterions blackbox blackboxes bb eval evals evaluation evaluations cache  "  , "false" , "true" , "true" },
{ "MAX_SURROGATE_EVAL_OPTIMIZATION",  "size_t",  "INF",  " Stopping criterion on the number of static surrogate evaluations ",  " \n  \n . Maximum number of static surrogate evaluations \n  \n . Argument: one positive integer. \n  \n . An INF value serves to disable the stopping criterion. \n  \n . Does not consider evaluations taken in the cache (cache hits) \n  \n . Only used when EVAL_SURROGATE_OPTIMIZATION is true. \n  \n . Example: MAX_SURROGATE_EVAL_OPTIMIZATION 1000 \n  \n . Default: INF\n\n",  "  basic stop max maximum surrogate  "  , "false" , "false" , "true" },
{ "MODEL_MAX_BLOCK_SIZE",  "size_t",  "INF",  " Internal parameter for QUAD_MODEL_MAX_BLOCK_SIZE and SGTELIB_MODEL_MAX_BLOCK_SIZE ",  " \n . CANNOT BE MODIFIED BY USER. Internal parameter. \n . Default: INF\n\n",  "  internal  "  , "true" , "true" , "true" },
{ "MODEL_MAX_EVAL",  "size_t",  "1000",  " Internal parameter for QUAD_MODEL_MAX_EVAL and SGTELIB_MODEL_MAX_EVAL ",  " \n . CANNOT BE MODIFIED BY USER. Internal parameter. \n . Default: 1000\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "QUAD_MODEL_MAX_BLOCK_SIZE",  "size_t",  "INF",  " Size of blocks of points, to be used for parallel evaluations ",  " \n . Maximum size of a block of evaluations send to the quad model evaluator \n   at once. Opportunistic strategies may apply after each block of evaluations. \n  \n . Depending on the algorithm phase, the quad model will \n   receive at most QUAD_MODEL_MAX_BLOCK_SIZE points to evaluate. \n  \n . When this parameter is greater than one, the number of evaluations \n   may exceed the QUAD_MODEL_MAX_EVAL stopping criterion. \n  \n . Argument: integer > 0. \n  \n . Example: QUAD_MODEL_MAX_BLOCK_SIZE 100 \n            The blackbox executable receives blocks of \n            at most 100 points for evaluation. \n  \n . Default: INF\n\n",  "  advanced block parallel  "  , "true" , "true" , "true" },
{ "QUAD_MODEL_MAX_EVAL",  "size_t",  "2000",  " Max number of model evaluations for optimization of the quad model problem ",  " \n . Max number of model evaluations for each optimization of the quad model problem. \n  \n . Argument: one integer > 0. \n  \n . The number of evaluations may exceed this parameter when parameter \n   QUAD_MODEL_MAX_BLOCK_SIZE is greater than one. \n    \n . For faster quadratic model search execution this value can be reduced to 1000. \n   However, this can effect the quality of trial points obtained by quad model \n   search. \n  \n . Example: QUAD_MODEL_MAX_EVAL 5000 \n . Default: 2000\n\n",  "  advanced quad search model model_search  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_MAX_BLOCK_SIZE",  "size_t",  "INF",  " Size of blocks of points, to be used for parallel evaluations ",  " \n . Maximum size of a block of evaluations send to the sgtelib model evaluator \n   at once. Opportunistic strategies may apply after each block of evaluations. \n  \n . Depending on the algorithm phase, the sgtelib will \n   receive at most SGTELIB_MODEL_MAX_BLOCK_SIZE points to evaluate. \n  \n . When this parameter is greater than one, the number of evaluations \n   may exceed the SGTELIB_MODEL_MAX_EVAL stopping criterion. \n  \n . Argument: integer > 0. \n  \n . Example: SGTELIB_MODEL_MAX_BLOCK_SIZE 100 \n            The blackbox executable receives blocks of \n            at most 100 points for evaluation. \n  \n . Default: INF\n\n",  "  advanced block parallel  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_MAX_EVAL",  "size_t",  "2000",  " Max number of model evaluations for each optimization of the sgtelib model problem ",  " \n . Max number of model evaluations for each optimization of the sgtelib model problem. \n  \n . Argument: one integer > 0.  \n  \n . The number of evaluations may exceed this parameter when parameter SGTELIB_MODEL_MAX_BLOCK_SIZE \n   is greater than one. \n  \n . Example: SGTELIB_MODEL_MAX_EVAL 5000 \n . Default: 2000\n\n",  "  advanced sgtelib search model model_search  "  , "true" , "true" , "true" },
{ "TMP_DIR",  "std::string",  "",  " Directory where to put temporary files ",  " \n  \n . Temporary directory for blackbox input/output files \n  \n . Argument: one string indicating a directory \n  \n . Improved performance with a local temporary directory \n  \n . Default is empty. Input/output files are put in problem directory \n  \n . Example: TMP_DIR /tmp \n  \n . Default: Empty string.\n\n",  "  advanced  "  , "false" , "false" , "true" },
{ "USE_CACHE_FILE_FOR_RERUN",  "bool",  "false",  " Cache file for rerun ",  " \n  \n . Flag to allow a cache set for rerun optimization. \n  \n . The cache for rerun is filled with evaluation points in the cache file (if  \n   it is provided, otherwise the option is inactive). \n    \n . The cache set for rerun complements the regular cache set. The regular \n   is used for testing if a point has already be evaluated  \n   (but not only for that). \n    \n . If an algorithm proposes a new trial point (not in regular cache), that is in \n   the cache for rerun, the evaluation results will be used. \n    \n . Points not in cache for rerun will be evaluated with blackbox. This allow to \n   perform a hot restart or reset the state of an algorithm to suggest new points. \n  \n . If no cache file is provided simply revert to blackbox evaluation. \n  \n . Argument: bool \n  \n . Example: USE_CACHE_FILE_FOR_RERUN true \n  \n . Default: false\n\n",  "  advanced  "  , "false" , "false" , "true" },
{ "NB_THREADS_PARALLEL_EVAL",  "int",  "1",  " Max number of threads used for parallel evaluations of each algorithm ",  " \n . A value greater than 1 is possible only if code is built with OpenMP enabled. \n  \n . Argument: One positive integer, or -1. A value of -1 means OpenMP decides \n   by itself. \n  \n . Example: NB_THREADS_PARALLEL_EVAL 4 \n . Default: 1\n\n",  "  advanced parallel openmp omp evaluator evaluation  "  , "true" , "false" , "true" } };

#endif
