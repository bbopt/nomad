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

#ifndef __NOMAD_4_5_DEPRECATEDATTRIBUTESDEFINITION__
#define __NOMAD_4_5_DEPRECATEDATTRIBUTESDEFINITION__

_definition = {
{ "ASYNCHRONOUS",  "bool",  "true",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: true\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "BB_INPUT_INCLUDE_SEED",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "BB_INPUT_INCLUDE_TAG",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "CACHE_SAVE_PERIOD",  "size_t",  "25",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 25\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "CACHE_SEARCH",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "CLOSED_BRACE",  "std::string",  "}",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: }\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "DISABLE",  "NOMAD::ArrayOfString",  "",  " Deprecated. DISABLE MODELS is replaced by QUAD_MODEL_SEARCH false and SGTELIB_MODEL_SEARCH false. DISABLE EVAL_SORT is replaced by EVAL_QUEUE_SORT LEXICOGRAPHICAL.",  " \n . Default: Empty string.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "EXTENDED_POLL_ENABLED",  "bool",  "true",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: true\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "EXTENDED_POLL_TRIGGER",  "NOMAD::Double",  "0.1",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 0.1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "F_TARGET",  "NOMAD::Double",  "0.0",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 0.0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "HAS_SGTE",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "INITIAL_MESH_INDEX",  "int",  "0",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "INITIAL_POLL_SIZE",  "NOMAD::ArrayOfDouble",  "-",  " Deprecated from Nomad 3: replaced by INITIAL_FRAME_SIZE ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "INTENSIFICATION_TYPE",  "std::string",  "POLL",  " Deprecated from Nomad 3: Intensification not implemented ",  " \n . Default: POLL\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "INT_POLL_DIR_TYPES",  "NOMAD::DirectionTypeList",  "ORTHO 1",  " Deprecated from Nomad 3: Intensification not implemented ",  " \n . Default: ORTHO 1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "L_CURVE_TARGET",  "NOMAD::Double",  "-",  " Deprecated from Nomad 3: Not implemented ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MAX_CACHE_MEMORY",  "size_t",  "2000",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 2000\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MAX_CONSECUTIVE_FAILED_ITERATIONS",  "int",  "-1",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: -1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MAX_EVAL_INTENSIFICATION",  "int",  "-1",  " Deprecated from Nomad 3: Intensification not implemented ",  " \n . Default: -1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MAX_SGTE_EVAL",  "size_t",  "1000",  " Deprecated from Nomad 3: replaced by QUAD_MODEL_MAX_EVAL and SGTELIB_MODEL_MAX_EVAL ",  " \n . Default: 1000\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MAX_SIM_BB_EVAL",  "int",  "-1",  " Deprecated from Nomad 3: not implemented ",  " \n . Default: -1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MESH_COARSENING_EXPONENT",  "size_t",  "1",  " Deprecated from Nomad 3: Only GMesh is used in Nomad 4 ",  " \n . Default: 1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MESH_REFINING_EXPONENT",  "int",  "-1",  " Deprecated from Nomad 3: Only GMesh is used in Nomad 4 ",  " \n . Default: -1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MESH_TYPE",  "std::string",  "X",  " Deprecated from Nomad 3: Only GMesh is used in Nomad 4 ",  " \n . Default: X\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MESH_UPDATE_BASIS",  "NOMAD::Double",  "4.0",  " Deprecated from Nomad 3: Only GMesh is used in Nomad 4 ",  " \n . Default: 4.0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MIN_POLL_SIZE",  "NOMAD::ArrayOfDouble",  "-",  " Deprecated from Nomad 3: replaced by MIN_FRAME_SIZE ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_EVAL_SORT",  "bool",  "true",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: true\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_EVAL_SORT_CAUTIOUS",  "bool",  "false",  " Deprecated from Nomad 3: not present in Nomad 4 ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_NP1_QUAD_EPSILON",  "NOMAD::Double",  "0.01",  " Deprecated from Nomad 3: ortho n+1 not implemented  ",  " \n . Default: 0.01\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_QUAD_MAX_Y_SIZE",  "size_t",  "500",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 500\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_QUAD_MIN_Y_SIZE",  "int",  "2",  " Deprecated from Nomad 3: not present in Nomad 4 ",  " \n . Default: 2\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_QUAD_RADIUS_FACTOR",  "NOMAD::Double",  "2.0",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 2.0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_QUAD_USE_WP",  "bool",  "false",  " Deprecated from Nomad 3: not used anymore ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_RADIUS_FACTOR",  "NOMAD::Double",  "2.0",  " Deprecated from Nomad 3: replaced by QUAD_MODEL_SEARCH_xxx_FACTOR \n  in Nomad 4 ",  " \n . Default: 2.0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_SEARCH",  "std::string",  "QUADRATIC",  " Deprecated from Nomad 3: replaced by QUAD_MODEL_SEARCH and \n   SGTELIB_MODEL_SEARCH in Nomad 4 ",  " \n . Default: QUADRATIC\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_SEARCH_MAX_TRIAL_PTS",  "size_t",  "10",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 10\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_SEARCH_OPTIMISTIC",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_SEARCH_OPPORTUNISTIC",  "bool",  "false",  " Deprecated from Nomad 3: Model search is never opportunistic in Nomad 4 ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MODEL_SEARCH_PROJ_TO_MESH",  "bool",  "true",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: true\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MULTI_F_BOUNDS",  "NOMAD::ArrayOfDouble",  "-",  " Deprecated from Nomad 3: Not implemented ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MULTI_FORMULATION",  "std::string",  "PRODUCT",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: PRODUCT\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MULTI_NB_MADS_RUNS",  "size_t",  "INF",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: INF\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MULTI_OVERALL_BB_EVAL",  "size_t",  "INF",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: INF\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "MULTI_USE_DELTA_CRIT",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NEIGHBORS_EXE",  "std::string",  "",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: Empty string.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_DELTA_E",  "NOMAD::Double",  "2",  " Deprecated from Nomad 3: replaced by NM_DELTA_E ",  " \n . Default: 2\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_DELTA_IC",  "NOMAD::Double",  "-0.5",  " Deprecated from Nomad 3: replaced by NM_DELTA_IC ",  " \n . Default: -0.5\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_DELTA_OC",  "NOMAD::Double",  "0.5",  " Deprecated from Nomad 3: replaced by NM_DELTA_OC ",  " \n . Default: 0.5\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_GAMMA",  "NOMAD::Double",  "0.5",  " Deprecated from Nomad 3: replaced by NM_GAMMA (valid for search and opt. )",  " \n . Default: 0.5\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_INCLUDE_FACTOR",  "NOMAD::Double",  "8",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 8\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_INIT_Y_BEST_VON",  "bool",  "false",  " Deprecated from Nomad 3: always true in Nomad 4 ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_INIT_Y_ITER",  "bool",  "false",  " Deprecated from Nomad 3: not used anymore ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_INTENSIVE",  "bool",  "false",  " Deprecated from Nomad 3: not used anymore ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_MAX_TRIAL_PTS",  "int",  "-1",  " Deprecated from Nomad 3: use NM_SEARCH_MAX_TRIAL_PTS_NFACTOR instead ",  " \n . Default: -1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_MIN_SIMPLEX_VOL",  "size_t",  "0",  " Deprecated from Nomad 3: stopping criterion not used anymore ",  " \n . Default: 0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_OPPORTUNISTIC",  "bool",  "false",  " Deprecated from Nomad 3: replaced by NM_SEARCH_STOP_ON_SUCCESS ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_SCALED_DZ",  "bool",  "false",  " Deprecated from Nomad 3: always true in Nomad 4 ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_USE_ONLY_Y",  "bool",  "false",  " Deprecated from Nomad 3: not used anymore ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NM_SEARCH_USE_SHORT_Y0",  "bool",  "false",  " Deprecated from Nomad 3: not used anymore ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPEN_BRACE",  "std::string",  "{",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: {\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPPORTUNISTIC_CACHE_SEARCH",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPPORTUNISTIC_EVAL",  "NOMAD::ArrayOfDouble",  "-",  " Deprecated from Nomad 3: replaced by EVAL_OPPORTUNISTIC ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPPORTUNISTIC_LH",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPPORTUNISTIC_LUCKY_EVAL",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPPORTUNISTIC_MIN_EVAL",  "size_t",  "INF",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: INF\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPPORTUNISTIC_MIN_F_IMPRVMT",  "NOMAD::Double",  "-",  " Deprecated from Nomad 3: Not implemented ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPPORTUNISTIC_MIN_NB_SUCCESS",  "size_t",  "INF",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: INF\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "OPT_ONLY_SGTE",  "bool",  "false",  " Deprecated from Nomad 3: replaced by EVAL_SURROGATE_OPTIMIZATION ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "PERIODIC_VARIABLE",  "NOMAD::ArrayOfDouble",  "-",  " Deprecated from Nomad 3: Not implemented ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "POINT_DISPLAY_LIMIT",  "size_t",  "20",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 20\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "POLL_UPDATE_BASIS",  "NOMAD::Double",  "2.0",  " Deprecated from Nomad 3: Only GMesh is used in Nomad 4 ",  " \n . Default: 2.0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "RANDOM_EVAL_SORT",  "bool",  "false",  " Deprecated from Nomad 3: replaced by EVAL_QUEUE_SORT RANDOM ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "ROBUST_MADS",  "bool",  "false",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "ROBUST_MADS_STANDARD_DEV_FACTOR",  "NOMAD::Double",  "2",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: 2\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SCALING",  "NOMAD::ArrayOfDouble",  "-",  " Deprecated from Nomad 3: Not implemented ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SEC_POLL_DIR_TYPE",  "NOMAD::DirectionTypeList",  "DOUBLE",  " Deprecated from Nomad 3: replaced by DIRECTION_TYPE_SECONDARY_POLL ",  " \n . Default: DOUBLE\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTE_CACHE_FILE",  "std::string",  "",  " Deprecated from Nomad 3: Cache for static surrogate not implemented ",  " \n . Default: Empty string.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTE_COST",  "size_t",  "INF",  " Deprecated from Nomad 3: replaced by EVAL_SURROGATE_COST ",  " \n . Default: INF\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTE_EVAL_SORT",  "bool",  "true",  " Deprecated from Nomad 3: replaced by EVAL_QUEUE_SORT SURROGATE ",  " \n . Default: true\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTE_EXE",  "std::string",  "sgte.exe",  " Deprecated from Nomad 3: replaced by SURROGATE_EXE ",  " \n . Default: sgte.exe\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTELIB_MODEL_CANDIDATES_NB",  "int",  "-1",  " Deprecated from Nomad 3: replaced by SGTELIB_MODEL_SEARCH_CANDIDATES_NB in Nomad 4 ",  " \n . Default: -1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTELIB_MODEL_EVAL_NB",  "size_t",  "10000",  " Deprecated from Nomad 3: replaced by SGTELIB_MODEL_MAX_EVAL ",  " \n . Default: 10000\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTELIB_MODEL_EXCLUSION_AREA",  "NOMAD::Double",  "0.0",  " Deprecated from Nomad 3: replaced by SGTELIB_MODEL_SEARCH_EXCLUSION_AREA in Nomad 4 ",  " \n . Default: 0.0\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTELIB_MODEL_FILTER",  "std::string",  "2345",  " Deprecated from Nomad 3: replaced by SGTELIB_MODEL_SEARCH_FILTER in Nomad 4 ",  " \n . Default: 2345\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SGTELIB_MODEL_TRIALS",  "int",  "1",  " Deprecated from Nomad 3: replaced by SGTELIB_MODEL_SEARCH_TRIALS in Nomad 4 ",  " \n . Default: 1\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "SNAP_TO_BOUNDS",  "bool",  "true",  " Deprecated from Nomad 3: Not implemented ",  " \n . Default: true\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "TREND_MATRIX",  "NOMAD::ArrayOfDouble",  "-",  " Deprecated from Nomad 3: Trend matrix not implemented  ",  " \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "TREND_MATRIX_BASIC_LINE_SEARCH",  "bool",  "false",  " Deprecated from Nomad 3: Trend matrix not implemented  ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "TREND_MATRIX_EVAL_SORT",  "bool",  "false",  " Deprecated from Nomad 3: Trend matrix not implemented  ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "VNS_SEARCH",  "bool",  "false",  " Replaced by VNS_MADS_SEARCH ",  " \n . Default: false\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "NB_THREADS_OPENMP",  "int",  "1",  " Deprecated from Nomad 4.4. Check NB_THREADS_PARALLEL_EVAL.",  " \n      Replaced by NB_THREADS_PARALLEL_EVAL to control explictlty the number of \n      threads used for parallel evaluations in evaluator control. \n . Default: 1\n\n",  "  internal  "  , "false" , "false" , "true" } };

#endif
