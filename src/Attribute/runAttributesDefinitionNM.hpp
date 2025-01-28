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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONNM__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONNM__

_definition = {
{ "NM_OPTIMIZATION",  "bool",  "false",  " Nelder Mead stand alone optimization for constrained and unconstrained pbs ",  " \n  \n . Nelder Mead optimization for constrained and unconstrained optimization \n  \n . Argument: bool \n  \n . Stand alone Nelder Mead optimization will deactivate any optimization strategy. \n  \n . Example: NM_OPTIMIZATION true \n  \n . Default: false\n\n",  "  advanced nelder mead simplex  "  , "true" , "false" , "true" },
{ "NM_SEARCH",  "bool",  "true",  " Nelder Mead optimization used as a search step for Mads ",  " \n  \n . Nelder Mead optimization as a search step for Mads \n  \n . Argument: bool \n  \n . Example: NM_SEARCH false \n  \n . Default: true\n\n",  "  advanced nelder mead simplex mads search "  , "true" , "true" , "true" },
{ "NM_SIMPLEX_INCLUDE_LENGTH",  "NOMAD::Double",  "INF",  " Construct NM simplex using points in cache.",  " \n  \n . Construct NM simplex using points in cache within a given distance of poll \n   center in absolute value. See also NM_SIMPLEX_INCLUDE_FACTOR. \n  \n . Argument: Positive double. INF means all points are considered. \n  \n . Example: NM_SIMPLEX_INCLUDE_LENGTH 0.2 \n  \n . Default: INF\n\n",  "  advanced nelder mead simplex length "  , "true" , "true" , "true" },
{ "NM_SIMPLEX_INCLUDE_FACTOR",  "size_t",  "8",  " Construct NM simplex using points in cache.",  " \n  \n . Construct NM simplex using points in cache within a given length of frame center \n   relative. The length equals the include factor multiplied by the frame size. \n   Used only if the mesh is defined. See also NM_SIMPLEX_INCLUDE_LENGTH. \n  \n . Argument: Positive integer.  \n  \n . Example: NM_SIMPLEX_INCLUDE_FACTOR 10 \n  \n . Default: 8\n\n",  "  advanced nelder mead simplex include factor length poll  "  , "true" , "true" , "true" },
{ "NM_DELTA_E",  "NOMAD::Double",  "2",  " NM expansion parameter delta_e.",  " \n  \n . Nelder Mead expansion parameter \n  \n . Argument: Positive NOMAD::Double > 1 \n  \n . Example: NM_DELTA_E 2.5 \n  \n . Default: 2\n\n",  "  advanced nelder mead simplex expansion  "  , "true" , "true" , "true" },
{ "NM_DELTA_IC",  "NOMAD::Double",  "-0.5",  " NM inside contraction parameter delta_ic.",  " \n  \n . Nelder Mead inside contraction parameter \n  \n . Argument: Negative NOMAD::Double \n  \n . Example: NM_DELTA_IC -1 \n  \n . Default: -0.5\n\n",  "  advanced nelder mead simplex inside contraction  "  , "true" , "true" , "true" },
{ "NM_DELTA_OC",  "NOMAD::Double",  "0.5",  " NM outside contraction parameter delta_oc.",  " \n  \n . Nelder Mead outside contraction parameter \n  \n . Argument: Positive NOMAD::Double <= 1 \n  \n . Example: NM_DELTA_OC 0.8 \n  \n . Default: 0.5\n\n",  "  advanced nelder mead simplex outside contraction  "  , "true" , "true" , "true" },
{ "NM_GAMMA",  "NOMAD::Double",  "0.5",  " NM shrink parameter gamma.",  " \n  \n . Nelder Mead shrink parameter \n  \n . Argument: Positive NOMAD::Double <= 1 \n  \n . Example: NM_GAMMA 0.8 \n  \n . Default: 0.5\n\n",  "  advanced nelder mead simplex shrink  "  , "true" , "true" , "true" },
{ "NM_SEARCH_MAX_TRIAL_PTS_NFACTOR",  "size_t",  "80",  " NM-Mads search stopping criterion.",  " \n  \n . NM-Mads stopping criterion. Max number of trial pts < dimension * NFactor \n  \n . Argument: Positive integer. INF disables this criterion. \n  \n . Example: NM_SEARCH_MAX_TRIAL_PTS_NFACTOR 100 \n  \n . Default: 80\n\n",  "  advanced nelder mead mads search stop trial  "  , "true" , "true" , "true" },
{ "NM_SEARCH_RANK_EPS",  "NOMAD::Double",  "0.01",  " NM-Mads epsilon for the rank of DZ.",  " \n  \n . Precision to detect when a vector increases the rank or not. \n  \n . Argument: Positive double. \n  \n . Example: NM_SEARCH_RANK_EPS 1E-4 \n  \n . Default: 0.01\n\n",  "  advanced nelder mead mads search rank DZ  "  , "true" , "true" , "true" },
{ "NM_SEARCH_STOP_ON_SUCCESS",  "bool",  "false",  " NM-Mads search stops on success.",  " \n  \n . NM-Mads search opportunistically stops on success. \n  \n . Argument: boolean. \n  \n . Example: NM_SEARCH_STOP_ON_SUCCESS false \n  \n . Default: false\n\n",  "  advanced nelder mead mads search opportunistic success  "  , "true" , "true" , "true" } };

#endif
