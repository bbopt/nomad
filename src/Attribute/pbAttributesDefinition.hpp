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

#ifndef __NOMAD_4_5_PBATTRIBUTESDEFINITION__
#define __NOMAD_4_5_PBATTRIBUTESDEFINITION__

_definition = {
{ "BB_INPUT_TYPE",  "NOMAD::BBInputTypeList",  "* R",  " The variable blackbox input types ",  " \n  \n . Blackbox input types \n  \n . List of types for each variable \n  \n . Available types: \n     . B: binary \n     . I: integer \n     . R: continuous \n  \n . Examples: \n     . BB_INPUT_TYPE * I       # all variables are integers \n     . BB_INPUT_TYPE ( R I B ) # for all 3 variables \n     . BB_INPUT_TYPE 1-3 B     # NOT YET SUPPORTED ( variables 1 to 3 are binary ) \n     . BB_INPUT_TYPE 0 I       # NOT YET SUPPORTED ( first variable is integer ) \n  \n . Default: * R\n\n",  "  basic blackbox blackboxes input inputs type types int integer integers binary bin continuous \n categorical  "  , "false" , "false" , "true" },
{ "DIMENSION",  "size_t",  "0",  " Dimension of the optimization problem (required) ",  " \n  \n . Number of variables \n  \n . Argument: one positive integer. \n  \n . Example: DIMENSION 3 \n  \n . Default: 0\n\n",  "  basic dimension dimensions dim dims problem problems prob pb pbs optimization size  "  , "false" , "false" , "true" },
{ "FIXED_VARIABLE",  "NOMAD::Point",  "-",  " Fix some variables to some specific values ",  " \n  \n . Fix some variables to some specific values \n  \n . Arguments: variable indexes and values \n  \n . Values for fixed variables are optional. Values of X0 will be used. \n  \n . Examples: \n     . FIXED_VARIABLE ( 0.0 - 0.0 )  # Variables 0 and 2 are fixed to value 0.0. \n                                     # Variable 1 is not fixed. \n  \n     . FIXED_VARIABLE 0              # Variable 0 is fixed to its X0 value. \n     . FIXED_VARIABLE 2-4            # Variables 2, 3 and 4 are fixed \n                                     # to their X0 values. \n . No default value.\n\n",  "  advanced fixed variable variables  "  , "false" , "false" , "true" },
{ "GRANULARITY",  "NOMAD::ArrayOfDouble",  "-",  " The granularity of the variables ",  " \n  \n . Set the granularity of variables to some specific values \n  \n . Arguments: granularity indexes and values (positive) \n  \n . Examples: \n     . GRANULARITY ( 0.01 0.0 0.01 ) # granularity of variables 0 and 2 are \n                                     # set to 0.01. \n       Variable 1 is real. \n     . GRANULARITY 0-1 0.01          # 2 first variables set granularity to 0.01 \n     . GRANULARITY * 0.01            # all variables set to granularity 0.01 \n  \n . No default value.\n\n",  "  advanced granular integer integers variable variables step  "  , "false" , "false" , "true" },
{ "INITIAL_FRAME_SIZE",  "NOMAD::ArrayOfDouble",  "-",  " The initial frame size of MADS ",  " \n  \n . Initial frame size \n  \n . Arguments: one or DIMENSION positive real(s) \n  \n . Reinterpreted empty default: \n     10% of the range if bounds are defined, |x0|/10 otherwise \n  \n . NOMAD uses one frame size per variable to achieve scaling \n  \n . The initial mesh size is determined from initial frame size when provided, but \n providing both is not allowed. \n  \n . Examples \n . INITIAL_FRAME_SIZE * 1.0          # for all variables \n . INITIAL_FRAME_SIZE 1 0.5        # for variable 1 only \n . INITIAL_FRAME_SIZE 2-4 0.25     # for variable 2 to 4 \n  \n . No default value.\n\n",  "  advanced intial poll frame mesh size mads gmesh  "  , "false" , "false" , "true" },
{ "INITIAL_MESH_SIZE",  "NOMAD::ArrayOfDouble",  "-",  " The initial mesh size of MADS ",  " \n  \n . Initial mesh size \n  \n . Arguments: one or DIMENSION positive real(s) \n  \n . NOMAD uses one mesh size per variable. \n  \n . Initial frame size is determined from initial mesh size when provided \n  \n . Examples: \n     . INITIAL_MESH_SIZE * 1.0          # for all variables \n     . INITIAL_MESH_SIZE 1 0.5        # for variable 1 only \n     . INITIAL_MESH_SIZE 2-4 0.25     # for variables 2 to 4 \n  \n . No default value.\n\n",  "  advanced initial mesh size mads gmesh  "  , "false" , "false" , "true" },
{ "LOWER_BOUND",  "NOMAD::ArrayOfDouble",  "-",  " The optimization problem lower bounds for each variable ",  " \n  \n . Lower bounds for each variable \n  \n . Arguments: DIMENSION reals \n  \n . Examples: \n     LOWER_BOUND * 0.0   # all variables are nonnegative \n     LOWER_BOUND 0-2 0.0 # the 3 first variables are nonnegative \n     LOWER_BOUND 0 0.0   # the first variable is nonnegative \n  \n . No default value.\n\n",  "  basic bound bounds lower variable variables constraint constraints  "  , "false" , "false" , "true" },
{ "MIN_FRAME_SIZE",  "NOMAD::ArrayOfDouble",  "-",  " Termination criterion on minimal frame size of MADS ",  " \n  \n . Minimum frame size. Can be set explicitely or automatically to 1 for \n   integer or binary variables (during check). \n  \n . Arguments: same logic as INITIAL_FRAME_SIZE ('r' can be used) \n  \n . Example: MIN_FRAME_SIZE r1E-5 \n  \n . No default value.\n\n",  "  advanced min minimum poll frame size stop stopping terminate terminates \n termination terminations mads  "  , "false" , "false" , "true" },
{ "MIN_MESH_SIZE",  "NOMAD::ArrayOfDouble",  "-",  " Termination criterion on minimal mesh size of MADS ",  " \n  \n . Minimum mesh size \n  \n . Arguments: same logic as INITIAL_MESH_SIZE ('*' can be used) \n  \n . Example: MIN_MESH_SIZE * 1E-5 \n  \n . No default value.\n\n",  "  advanced min minimum frame mesh size stop stopping terminate terminates termination terminations \n mads  "  , "false" , "false" , "true" },
{ "UPPER_BOUND",  "NOMAD::ArrayOfDouble",  "-",  " The optimization problem upper bounds for each variable ",  " \n  \n . Upper bounds for each variable \n  \n . Arguments: DIMENSION reals \n  \n . Examples: \n     UPPER_BOUND * 10.0   # all variables are less than or equal to 10.0 \n     UPPER_BOUND 0-2 10.0 # the 3 first variables are less than or equal to 10.0 \n     UPPER_BOUND 0 10.0   # the first variable is less than or equal to 10.0 \n  \n . No default value.\n\n",  "  basic bound bounds upper variable variables constraint constraints  "  , "false" , "false" , "true" },
{ "VARIABLE_GROUP",  "NOMAD::ListOfVariableGroup",  "-",  " The groups of variables) ",  " \n  \n . When defined, Mads poll generates trial points in a separate subspace for \n   each group of variables. \n  \n . Arguments: a list of variable indices or a range of indices \n  \n . More than one group of variables can be defined. Providing a single group of \n   variables puts the remaining indices into another variable group. \n  \n . The indices in the groups of variables must be unique. \n  \n . Examples: \n     . VARIABLE_GROUP 0 1 2 6 \n     . VARIABLE_GROUP 2-6 \n  \n . No default value.\n\n",  "  advanced group groups variable variables var vars poll  "  , "true" , "false" , "false" },
{ "X0",  "NOMAD::ArrayOfPoint",  "-",  " The initial point(s) ",  " \n  \n . Vector of starting point(s) \n  \n . Arguments: text file name or DIMENSION reals \n  \n . More than one starting point can be defined in a separate text file \n   (X0 x0s.txt) with one point per line. A single point can be provided in the \n   parameter file: X0 (0 0 0). \n  \n . All points are evaluated: X0 evaluations are not opportunistic. \n   Initial LH_SEARCH points are handled as X0s (no opportunism for evaluation). \n  \n . May be infeasible \n  \n . Cannot be outside bounds \n  \n . Must respect fixed variables (parameter FIXED_VARIABLE) \n  \n . Examples: \n     . X0 x0.txt \n  \n     . X0   * 0.0    # First starting point \n       X0 1 * 1.0    # Second starting point \n  \n     . X0 ( 0 1 2 )  # if DIMENSION = 3 \n  \n . No default value.\n\n",  "  basic initial variable variables var vars init point points bound bounds lower upper cache start starting  "  , "false" , "false" , "false" },
{ "POINT_FORMAT",  "NOMAD::ArrayOfDouble",  "-",  " Format of the doubles for trial points",  " \n  \n . POINT_FORMAT and BB_EVAL_FORMAT have the same values \n  \n . POINT_FORMAT and BB_EVAL_FORMAT are computed from the BB_INPUT_TYPE parameter. \n  \n . Gives the format precision for doubles of trial points. \n  \n . CANNOT BE MODIFIED BY USER. Internal parameter. \n  \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" } };

#endif
