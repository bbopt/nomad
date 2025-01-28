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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONLH__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONLH__

_definition = {
{ "LH_EVAL",  "size_t",  "0",  " Latin Hypercube Sampling of points (no optimization) ",  " \n  \n . Latin-Hypercube sampling (evaluations) \n  \n . Argument: A positive integer p < INF.  \n  \n . p: number of LH points \n  \n . All points will be evaluated (no opportunism). \n  \n . This option will not work with Mads but can be combined with quadratic \n   model optimization to have enough points to construct models. \n  \n . The LH sampling requires to have both lower and upper bounds defined. \n  \n . Example: LH_EVAL 100 \n  \n . Default: 0\n\n",  "  basic latin hypercube sampling  "  , "true" , "true" , "true" },
{ "LH_SEARCH",  "NOMAD::LHSearchType",  "-",  " Latin Hypercube Sampling Search method ",  " \n  \n . Latin-Hypercube sampling (search) \n  \n . Arguments: two size_t p0 and pi < INF. \n  \n . p0: number of initial LH search points. These points are handled as X0s \n   (no opportunism for evaluation). \n  \n . pi: LH search points at each iteration. The iteration search can be \n   opportunistic or not (parameter EVAL_OPPORTUNISTIC). \n  \n . Example: LH_SEARCH 100 0 \n  \n . No default value.\n\n",  "  basic search latin hypercube sampling  "  , "true" , "true" , "true" } };

#endif
