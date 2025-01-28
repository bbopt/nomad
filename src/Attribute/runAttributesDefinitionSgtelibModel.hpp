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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONSGTELIBMODEL__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONSGTELIBMODEL__

_definition = {
{ "SGTELIB_MODEL_EVAL",  "bool",  "0",  " Sgtelib Model Sampling of points ",  " \n  \n . Sgtelib Model sampling \n  \n . Argument: bool \n  \n . Best points are taken from the cache \n  \n . A model is computed \n  \n . The most promising points according to that model are evaluated \n  \n . No opportunism \n  \n . This option deactivates other optimization strategies. \n  \n . Example: SGTELIB_MODEL_EVAL true \n  \n . Default: 0\n\n",  "  advanced sgtelib model sampling  "  , "true" , "false" , "true" },
{ "SGTELIB_MODEL_SEARCH",  "bool",  "false",  " Model search using Sgtelib ",  " \n  \n . MADS model search, using Bastien Talgorn's Sgtelib using a model definition \n  \n . Argument: one boolean ('yes' or 'no') \n  \n . Disabled for more than 50 variables \n  \n . See SGTELIB_MODEL_DEFINITION for the definition of the model \n  \n . Example: SGTELIB_MODEL_SEARCH yes \n  \n . Default: false\n\n",  "  basic mads model search sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_DISPLAY",  "std::string",  "",  " Display of a model ",  " \n . Control the display of the sgtelib model search and sgtelib model optimization. \n   These details are only shown if DISPLAY_DEGREE is FULL (3) or more. \n  \n . Arguments: a string containing one or several of the following letters \n  \n . \"S\": General information on the model search or optimization \n  \n . \"F\": Details of the filter step \n  \n . \"O\": Details of the models optimization \n  \n . \"P\": Details of the projection \n  \n . \"U\": Details of the model update \n  \n . \"I\": Advancement of the model optimization \n  \n . \"X\": Display of all of the model evaluations \n  \n . Example: SGTELIB_MODEL_DISPLAY SPF # display the general information on the search \n                                        and on the filter and projection steps \n . Default: Empty string.\n\n",  "  developer advanced model quad sgtelib  "  , "false" , "false" , "true" },
{ "SGTELIB_MODEL_DEFINITION",  "NOMAD::ArrayOfString",  "",  " Definition of the Sgtelib model ",  " \n . Argument: Array of string that represent the Sgtelib model definition. See sgtelib manual. \n  \n . See SGTELIB_MODEL_SEARCH or SGTELIB_MODEL_EVAL to enable model use. \n  \n . Example: TYPE PRS DEGREE 1 # builds a linear model \n .          TYPE PRS DEGREE 2 # builds a quadratic model \n .          TYPE RBF          # builds an RBF model \n .          TYPE ENSEMBLE     # builds an ensemble of models \n            # builds a lowess model with local linear regression \n            # and optimized kernel shape: \n .          TYPE LOWESS DEGREE 1 KERNEL_COEF OPTIM \n .          # Variation that gives good results: \n            TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM KERNEL_COEF OPTIM RIDGE 0 METRIC AOECV \n . Default: Empty string.\n\n",  "  advanced sgtelib search model model_search interpolation regression  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_SEARCH_TRIALS",  "size_t",  "1",  " Max number of sgtelib model search failures before going to the poll step ",  " \n . Max number of sgtelib model search failures before going to the poll step. \n  \n . Argument: a positive integer. \n  \n . Note: The minimum between this parameter and MAX_ITERATION_PER_MEGAITERATION \n   will be used. \n  \n . Example: SGTELIB_MODEL_SEARCH_TRIALS 5 \n . Default: 1\n\n",  "  developer trials sgtelib model search  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_FORMULATION",  "NOMAD::SgtelibModelFormulationType",  "FS",  " Formulation of the sgtelib model problem ",  " \n . Formulation of the sgtelib model problem. \n  \n . Argument: one string in {'FS', 'EIS', 'FSP', \n                            'EFI', 'EFIS','EFIM','EFIC', \n                            'PFI', \n                            'D', \n                            'EXTERN'} \n  \n . Description of the sgtelib problem formulations : \n     (FS)   min f    -d.sigma_f \n            st  c_j  -d.sigma_j <= 0 \n  \n     (EIS)  min -EI  -d.sigma_f \n            st  c_j  -d.sigma_j <= 0 \n  \n     (FSP)  min f    -d.sigma_f \n            st  P >= 1/2 \n  \n     (EFI)  min -EFI \n  \n     (EFIS) min -EFI -d.sigma_f \n  \n     (EFIM) min -EFI -d.sigma_f.mu \n  \n     (EFIM) min -EFI -d.(EI.mu+P.sigma_f) \n  \n     (PFI)  min -PFI \n  \n     (D)    min -distance_to_closest \n  \n . Example: SGTELIB_MODEL_FORMULATION EFI \n . Default: FS\n\n",  "  developer advanced problem expected improvement diversification model sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_FEASIBILITY",  "NOMAD::SgtelibModelFeasibilityType",  "C",  " Method used to model the feasibility of a point ",  " \n . Method used to model the feasibility of a point. \n  \n . Arguments: one character in {'C', 'H', 'M', 'P'} \n  \n . Example: SGTELIB_MODEL_FEASIBILITY C  # 1 model per constraint \n            SGTELIB_MODEL_FEASIBILITY H  # 1 model of the aggregate constraint \n            SGTELIB_MODEL_FEASIBILITY M  # 1 model of the max of the constraints \n            SGTELIB_MODEL_FEASIBILITY B  # 1 binary model of the feasibility \n . Default: C\n\n",  "  developer advanced feasibility constraints interpolation regression model sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_DIVERSIFICATION",  "NOMAD::Double",  "0.01",  " Coefficient of the exploration term in the sgtelib model problem ",  " \n . Coefficient of the exploration term in the sgtelib model problem. \n  \n . Argument: one positive real \n  \n . Example: SGTELIB_MODEL_DIVERSIFICATION 0    # no exploration \n            SGTELIB_MODEL_DIVERSIFICATION 0.01 # light exploration \n            SGTELIB_MODEL_DIVERSIFICATION 0.1  # medium exploration \n            SGTELIB_MODEL_DIVERSIFICATION 1    # strong exploration \n . Default: 0.01\n\n",  "  developer advanced model sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_SEARCH_EXCLUSION_AREA",  "NOMAD::Double",  "0.0",  " Exclusion area for the sgtelib model search around points of the cache ",  " \n . Defines an exclusion area for the sgtelib model search around points of the cache \n  \n . Arguments: one real number in [0, 0.5] \n  \n . Example: SGTELIB_MODEL_SEARCH_EXCLUSION_AREA 0 # no exclusion area \n            SGTELIB_MODEL_SEARCH_EXCLUSION_AREA 0.1 # small exclusion area \n            SGTELIB_MODEL_SEARCH_EXCLUSION_AREA 0.5 # large exclusion area \n  \n . Default: 0.0\n\n",  "  developer advanced model sgtelib search exclusion  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_SEARCH_CANDIDATES_NB",  "int",  "-1",  " Number of candidates returned by the sgtelib model search ",  " \n . Number of candidates returned by the sgtelib model search. \n  \n . Argument: one integer \n  \n . If smaller or equal to 0, then the number of candidates \n   will be the largest value between BB_MAX_BLOCK_SIZE and \n   2 * DIMENSION \n  \n . Example: SGTELIB_MODEL_SEARCH_CANDIDATES_NB 8 \n . Default: -1\n\n",  "  developer advanced model sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MIN_POINTS_FOR_MODEL",  "size_t",  "1",  " Minimum number of valid points necessary to build a model ",  " \n . Defines the minimum number of valid points beyond which no model will \n   be build \n  \n . Argument: a positive integer < INF. \n  \n . Example: SGTELIB_MIN_POINTS_FOR_MODEL 5 \n  \n . Default: 1\n\n",  "  developer advanced model sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MAX_POINTS_FOR_MODEL",  "size_t",  "500",  " Maximum number of valid points used to build a model ",  " \n . Defines the maximum number of valid points kept to build a model. \n   Extra points too far from center are ignored. \n  \n . Arguments: one positive integer.  \n  \n . Example: SGTELIB_MAX_POINTS_FOR_MODEL 96 \n  \n . Default: 500\n\n",  "  developer advanced model sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_SEARCH_FILTER",  "std::string",  "2345",  " Methods used in the sgtelib search filter to return several search candidates ",  " \n . Methods used in the sgtelib search filter to return several search candidates \n  \n . Arguments: a string containing several integers from 0 to 5 \n  \n . Method 0: Select the best candidate \n  \n . Method 1: Select the most remote candidate \n  \n . Method 2: Select the best candidate, with minimal distance to the cache \n  \n . Method 3: Select the best candidate, with minimal margin in feasibility \n  \n . Method 4: Select the candidate with the best isolation number \n  \n . Method 5: Select the candidate with the best density number \n  \n . Examples: SGTELIB_MODEL_SEARCH_FILTER 0    # Only method 0 will be used \n             SGTELIB_MODEL_SEARCH_FILTER 01   # Alternate between method 0 and 1 \n             SGTELIB_MODEL_SEARCH_FILTER 2345 # Cycle through methods 2, 3, 4 and 5 \n . Default: 2345\n\n",  "  developer advanced model search sgtelib  "  , "true" , "true" , "true" },
{ "SGTELIB_MODEL_RADIUS_FACTOR",  "NOMAD::Double",  "2.0",  " Sgtelib model radius factor ",  " \n . Sgtelib model radius factor \n  \n . This parameter is used to select points to build the sgtelib model \n  \n . Frame size is multiplied by this factor to get the search radius \n  \n . Points inside a circle centered on the poll center, within this radius, \n   are selected to build the sgtelib model \n  \n . Arguments: one strictly positive real \n  \n . Example: SGTELIB_MODEL_RADIUS_FACTOR 1.0 \n . Default: 2.0\n\n",  "  developer sgtelib model radius  "  , "true" , "true" , "true" } };

#endif
