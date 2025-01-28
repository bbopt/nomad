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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONQUADMODEL__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONQUADMODEL__

_definition = {
{ "QUAD_MODEL_SEARCH",  "bool",  "true",  " Quad model search ",  " \n  \n . MADS model search, using Bastien Talgorn's Sgtelib with quad models \n  \n . Argument: one boolean ('yes' or 'no') \n  \n . Disabled for more than 50 variables \n  \n . Example: QUAD_MODEL_SEARCH yes \n  \n . Default: true\n\n",  "  basic mads quad model search model_search  "  , "true" , "true" , "true" },
{ "QUAD_MODEL_SEARCH_SIMPLE_MADS",  "bool",  "false",  " Quad model search using a simpler version of Mads ",  " \n  \n . SimpleMads has only one poll method and no search method. \n  \n . The evaluator control is by-passed. The model is directly evaluated. \n  \n . SimpleMads use a simpler progressive barrier approach to handle constraints. \n  \n . This approach is more straightforward but maybe less efficient. \n  \n . Example: QUAD_MODEL_SEARCH_SIMPLE_MADS yes \n  \n . Default: false\n\n",  "  basic mads quad model search model_search  "  , "true" , "true" , "true" },
{ "QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR",  "NOMAD::Double",  "1",  " Scale the bounds for the quad model search  ",  " \n  \n . The quad model is built on evaluation points around a frame center. This \n defines min and max bounds, from which we define a model center. The model \n search is limited to tighter (can be wider with the parameter set to less than 1) \n bounds by reducing the distance of the optimization bounds to the model center. \n We use a reduction factor for that. If the reduction factor equals one, the min \n and max bounds are used as optimization bounds for the search. The greater the \n reduction factor, the tighter the bounds. \n  \n . Argument: one Double greater than 0 \n  \n . Example: QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR 3.0 \n  \n . Default: 1\n\n",  "  develop mads quad model search sgtelib model_search bounds  "  , "true" , "true" , "true" },
{ "QUAD_MODEL_DISPLAY",  "std::string",  "",  " Display of a model ",  " \n . Control the display of the quad model search and quad model optimization. \n   These details are only shown if DISPLAY_DEGREE is FULL (3) or more. \n  \n . Arguments: a string containing one or several of the following letters \n  \n . \"S\": General information on the model search or optimization \n  \n . \"F\": Details of the filter step \n  \n . \"O\": Details of the models optimization \n  \n . \"P\": Details of the projection \n  \n . \"U\": Details of the model update \n  \n . \"I\": Advancement of the model optimization \n  \n . \"X\": Display of all of the model evaluations \n  \n . Example: QUAD_MODEL_DISPLAY SPF # display the general information on the search \n                                        and on the filter and projection steps \n . Default: Empty string.\n\n",  "  developer advanced model quad sgtelib  "  , "false" , "false" , "true" },
{ "QUAD_MODEL_OPTIMIZATION",  "bool",  "false",  " Quad model stand alone optimization for constrained and unconstrained pbs ",  " \n  \n . Quadratic model optimization for constrained and unconstrained \n   optimization. \n    \n . LH_EVAL can be used to sample points and evaluate them before running \n   quad model optimization. \n  \n . Argument: bool \n  \n . Stand alone quadratic model optimization will deactivate any optimization \n   strategy. \n  \n . Example: QUAD_MODEL_OPTIMIZATION true \n  \n . Default: false\n\n",  "  advanced sgtelib quadratic quad optimization simplex  "  , "true" , "false" , "true" },
{ "QUAD_MODEL_SEARCH_BOX_FACTOR",  "NOMAD::Double",  "4.0",  " Quadratic model search point selection factor ",  " \n . Quadratic model search point selection factor \n  \n . This parameter is used to select points to build the quadratic model for \n   the search method \n  \n . The max, min and average of trial points coordinates are used to identify a \n   box. This box contains all trial points to sort. \n    \n . The box is enlarged by this factor. Cache evaluation points inside this box \n   are selected to build the quadratic model for search. \n  \n . Arguments: one strictly positive real \n  \n . Example: QUAD_MODEL_SEARCH_BOX_FACTOR 1.0 \n . Default: 4.0\n\n",  "  developer quadratic model search radius  "  , "true" , "true" , "true" },
{ "QUAD_MODEL_BOX_FACTOR",  "NOMAD::Double",  "4.0",  " Quadratic model points selection box factor ",  " \n  \n . This parameter is used to select points to build the quadratic model. \n  \n . The max, min and average of trial points coordinates are used to identify a \n   box. \n    \n . The box is enlarged by this factor. Cache evaluation points inside this box \n   are selected to build the quadratic model. \n  \n . Arguments: one strictly positive real \n  \n . Example: QUAD_MODEL_BOX_FACTOR 1.0 \n . Default: 4.0\n\n",  "  developer quadratic model sort search  "  , "true" , "true" , "true" },
{ "QUAD_MODEL_SEARCH_FORCE_EB",  "bool",  "false",  " Quadratic model search optimization using extreme barrier ",  " \n  \n . By default (flag=false) quadratic model search optimization problem transforms \n   all constraints into progressive barrier constraints. \n    \n . Handling constraints with the progressive barrier is more costly (time) than \n   with the extreme barrier because updating the progressive barrier is more \n   complex. \n  \n . On costly blackbox optimization problem the extra time for using the progressive \n   barrier during the quadratic model search optimization is minimal. But it can \n   be more important on simple analytical blacbox optimization problems with many \n   quad model searches. \n    \n . Forcing EB constraints for quad model search is faster but can result in less \n   successful trial points. \n    \n . Arguments: boolean \n  \n . Example: QUAD_MODEL_FORCE_EB true \n  \n . Default: false\n\n",  "  developer quadratic model search  "  , "true" , "true" , "true" } };

#endif
