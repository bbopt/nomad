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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONDMULTI__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONDMULTI__

_definition = {
{ "DMULTIMADS_EXPANSIONINT_LINESEARCH",  "bool",  "false",  " DMultiMads Expansion integer linesearch ",  " \n  \n . DMultiMads Expansion integer linesearch. Is only active when the problem \n   possesses integer variables. \n  \n . Argument: one boolean ('yes' or 'no') \n  \n . Example: DMULTIMADS_EXPANSIONINT_LINESEARCH yes \n  \n . Default: false\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" },
{ "DMULTIMADS_QUAD_MODEL_STRATEGY",  "NOMAD::DMultiMadsQuadSearchType",  "MULTI",  " Quad Model search strategies for DMultiMads ",  " \n . Quad model search strategy used by DMultiMads to generate trial points. \n  \n . Arguments: Quad model search strategy in: \n     DMS   : uses the Direct MultiSearch strategy proposed for the Direct \n     MultiSearch algorithm. The Quad Model search tries to find new non- \n     dominated points by minimizing subsets of objectives in a local region \n     centered around the current frame center; starting from one objective, \n     until minimizing all objectives at the same time. At each level of \n     combinations, the search checks if it finds a non-dominated solution. \n     Otherwise, it increases the number of objectives considered. In the \n     worst case, this strategy has to solve 2^nobj - 1 subproblems, that \n     can be costly. \n     DOM   : uses the Dominance Move strategy. In this case, the Quad model \n     search tries to minimize a 'dominance move' single-objective function, \n     and explores the objective space towards new non-dominated solutions \n     according to the current set of solutions. \n     MULTI : uses the MultiMads strategy. When the current frame center \n     is an extreme solution of the current set of solutions, i.e., it \n     minimizes one of the objective components, the Quad model search \n     starts an expansion phase: it then tries to minimize further the \n     considering objective component. Otherwise, the Quad model search \n     tries to find new non-dominated solutions in a local region of \n     the objective space defined by the current frame center and some \n     of its 'closest' neighbors. \n  \n . This parameter is active only when QUAD_MODEL_SEARCH is set to TRUE. \n   This is the case by default. \n  \n . Example: DMULTIMADS_QUAD_MODEL_STRATEGY DMS \n  \n . Default: MULTI\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" },
{ "DMULTIMADS_MIDDLEPOINT_SEARCH",  "bool",  "false",  " DMultiMads Middle Point search ",  " \n  \n . DMultiMads Middle Point search. \n  \n . Argument: one boolean ('yes' or 'no') \n  \n . Example: DMULTIMADS_MIDDLEPOINT_SEARCH no \n  \n . Default: false\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" },
{ "DMULTIMADS_MIDDLEPOINT_SEARCH_CACHE_MAX",  "size_t",  "50",  " DMultiMads middle point search ",  " \n  \n . Number of times the Middle Point search for DMultiMads is allowed to \n   consult the cache for each objective. When the maximum number of tentatives \n   is reached, no point is generated by the Middle Point search for the current \n   objective. \n  \n . Argument: one positive integer < INF. \n  \n . Example: DMULTIMADS_MIDDLEPOINT_SEARCH_CACHE_MAX 50 \n  \n . Default: 50\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" },
{ "DMULTIMADS_NM_STRATEGY",  "NOMAD::DMultiMadsNMSearchType",  "DOM",  " Nelder-Mead search strategies for DMultiMads ",  " \n . Nelder-Mead search strategy used by DMultiMads to generate trial points. \n  \n . Arguments: Nelder-Mead search strategy in: \n     DOM   : uses the Dominance Move strategy. In this case, the NM search \n     tries to minimize a 'dominance move' single-objective function, and \n     explores the objective space towards new non-dominated solutions \n     according to the current set of solutions. \n     MULTI : uses the MultiMads strategy. When the current frame center \n     is an extreme solution of the current set of solutions, i.e., it \n     minimizes one of the objective components, the NM search starts \n     an expansion phase: it then tries to minimize further the \n     considering objective component. Otherwise, the NM search tries to \n     find new non-dominated solutions in a local region of the objective \n     space defined by the current frame center and some of its \n     'closest' neighbors. \n  \n . This parameter is active only when NM_SEARCH is set to TRUE. \n   This is the case by default. \n  \n . Example: DMULTIMADS_NM_STRATEGY MULTI \n  \n . Default: DOM\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" },
{ "DMULTIMADS_OPTIMIZATION",  "bool",  "false",  " DMultiMads solves multiobjective optimization problems ",  " \n  \n . When BB_OUTPUT_TYPE contains more than one objective, DMultiMads algorithm \n must be explicitly enabled. \n  \n . ChT. TODO modify the note. Pareto front can be output in a file. \n . Note: It is recommended to output all evaluation points into a file (use \n HISTORY_FILE). The user must do some post-processing for these results to obtain \n  the approximation of the pareto front. Future version will have the option \n  to output the pareto front in a file. \n  \n . Argument: bool \n  \n . Example: DMULTIMADS_OPTIMIZATION true \n  \n . Default: false\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" },
{ "DMULTIMADS_SELECT_INCUMBENT_THRESHOLD",  "size_t",  "1",  " Control the choice of the DMultiMads incumbent ",  " \n  \n . Is set at this value at the initialization of the DMultiMads algorithm. \n  \n . The higher, the more potential candidates belonging to the current set of \n non dominated solutions can be potentialy chosen as the current incumbent at a \n given iteration. When DMULTIMADS_SELECT_INCUMBENT_THRESHOLD is equal to 0, only \n points possessing a maximal frame size value can be selected as the current \n incumbent. \n  \n . Argument: size_t \n  \n . Example: DMULTIMADS_SELECT_INCUMBENT_THRESHOLD 5 \n  \n . Default: 1\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" },
{ "DMULTIMADS_QMS_PRIOR_COMBINE_OBJ",  "bool",  "true",  " Select compute method for objective of DMultiMads quad model search ",  " \n  \n . TODO select the correct option and get rid of the parameter \n  \n . Default: true\n\n",  "  advanced multi objective algorithm pareto front  "  , "true" , "true" , "true" } };

#endif
