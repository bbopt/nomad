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

#ifndef __NOMAD_4_5_EVALUATORCONTROLATTRIBUTESDEFINITION__
#define __NOMAD_4_5_EVALUATORCONTROLATTRIBUTESDEFINITION__

_definition = {
{ "EVAL_OPPORTUNISTIC",  "bool",  "true",  " Opportunistic strategy: Terminate evaluations as soon as a success is found ",  " \n  \n . Opportunistic strategy: Terminate evaluations as soon as a success is found \n  \n . This parameter is the default value for other OPPORTUNISTIC parameters, \n    including Search steps \n  \n . This parameter is the value used for Poll step \n  \n . In addition, a custom criterion can also be provided (library mode only) in \n  a callback function (see example in \n  $NOMAD_HOME/examples/advanced/library/CustomOpportunistic ) \n  \n . Argument: one boolean (yes or no) \n  \n . Type 'nomad -h opportunistic' to see advanced options \n  \n . Example: EVAL_OPPORTUNISTIC no  # complete evaluations \n  \n . Default: true\n\n",  "  advanced opportunistic oppor eval evals evaluation evaluations terminate list success successes  "  , "true" , "true" , "true" },
{ "EVAL_SURROGATE_OPTIMIZATION",  "bool",  "false",  " Use static surrogate as blackbox for optimization ",  " \n  . Use solely static surrogate instead of the blackbox for optimization. \n  \n  . In batch mode, SURROGATE_EXE needs to be defined. \n  \n  . In library mode, an Evaluator for SURROGATE eval type needs to be defined. \n  \n  . Argument: bool \n  \n  . Example: EVAL_SURROGATE_OPTIMIZATION true \n  \n  . See also: SURROGATE_EXE, MAX_SURROGATE_EVAL_OPTIMIZATION \n  \n . Default: false\n\n",  "  advanced static surrogate  "  , "true" , "false" , "true" },
{ "EVAL_USE_CACHE",  "bool",  "true",  " Use cache in algorithms ",  " \n . When this parameter is false, the Cache is not used at all. Points may be \n   re-evaluated. \n  \n . Recommended when DIMENSION is large and evaluations are not costly. \n  \n . Cache may be used for top algorithm, and disabled for a sub-algorithm. \n  \n . If CACHE_FILE is non-empty, cache file will still be read and written. \n  \n . Default: true\n\n",  "  advanced  "  , "true" , "false" , "true" },
{ "EVAL_QUEUE_SORT",  "NOMAD::EvalSortType",  "QUADRATIC_MODEL",  " How to sort points before evaluation ",  " \n . Argument: One of DIR_LAST_SUCCESS, LEXICOGRAPHICAL, RANDOM, SURROGATE, \n   QUADRATIC_MODEL, USER \n  \n . DIR_LAST_SUCCESS: Points that are generated in a direction similar to the \n   last direction that provided a successful point are evaluated first. \n  \n . LEXICOGRAPHICAL: Points are sorted in lexicographical order before evaluation. \n  \n . RANDOM: Mix points randomly before evaluation, instead of sorting them. \n  \n . SURROGATE: Sort points using values given by static surrogate. See parameter SURROGATE_EXE. \n  \n . QUADRATIC_MODEL: Sort points using values given by dynamic quadratic models. \n  \n . USER: In library mode only. See example, $NOMAD_HOME/examples/advanced/library/CustomCompOrdering. \n This is set automatically when providing the user compare priority class. \n  \n . Example: EVAL_QUEUE_SORT LEXICOGRAPHICAL \n  \n . Default: QUADRATIC_MODEL\n\n",  "  advanced  "  , "true" , "true" , "true" },
{ "PSD_MADS_SUBPROBLEM_MAX_BB_EVAL",  "size_t",  "INF",  " Max number of evaluations for each subproblem ",  " \n  \n . Used in the context of Parallel Space Decomposition (PSD) MADS algorithm. \n  \n . Select the max number of evaluations in each Mads subproblem. \n  \n . Argument: a positive integer. \n  \n . Example: PSD_MADS_SUBPROBLEM_MAX_BB_EVAL 10 \n  \n . Default: INF\n\n",  "  advanced psd mads parallel decomposition subproblem  "  , "true" , "false" , "true" },
{ "SUBPROBLEM_MAX_BB_EVAL",  "size_t",  "INF",  " Internal parameter for PSD_MADS_SUBPROBLEM_MAX_BB_EVAL ",  " \n  \n . CANNOT BE MODIFIED BY USER. Internal parameter. \n  \n . Default: INF\n\n",  "  internal  "  , "false" , "false" , "true" } };

#endif
