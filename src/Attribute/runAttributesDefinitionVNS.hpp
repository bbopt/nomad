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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONVNS__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONVNS__

_definition = {
{ "VNS_MADS_OPTIMIZATION",  "bool",  "false",  " VNS MADS stand alone optimization for constrained and unconstrained pbs ",  " \n  \n . Shaking + optimization for constrained and unconstrained optimization \n  \n . Argument: bool \n  \n . Stand alone VNS Mads optimization will deactivate any optimization strategy. \n  \n . Example: VNS_MADS_OPTIMIZATION true \n  \n . Default: false\n\n",  "  advanced global optimization vns neighborhood  "  , "true" , "false" , "true" },
{ "VNSMART_MADS_SEARCH",  "bool",  "false",  " VNS Mads search under condition of consecutive fails ",  " \n . VNS Mads optimization used as a search step for Mads under condition \n   on the number of consecutive fails \n . Variable Neighborhood Search + Mads optimization as a search step for Mads \n . Criterion: threshold on the number of consecutive fails \n . Argument: bool \n . Example: VNSMART_MADS_SEARCH true \n . Default: false\n\n",  "  advanced global mads search vns neighborhood "  , "true" , "false" , "true" },
{ "VNSMART_MADS_SEARCH_THRESHOLD",  "int",  "3",  " Threshold for VNS (SMART) Mads search ",  " \n . Number of consecutive fails to enable VNS Mads search step \n . Variable Neighborhood Search + Mads optimization as a search step for Mads \n . Argument: int \n . Example: VNSMART_MADS_SEARCH_THRESHOLD 5 \n . Default: 3\n\n",  "  advanced global mads search vns neighborhood "  , "true" , "false" , "true" },
{ "VNS_MADS_SEARCH",  "bool",  "false",  " VNS Mads optimization used as a search step for Mads ",  " \n  \n . Variable Neighborhood Search + Mads optimization as a search step for Mads \n  \n . Argument: bool \n  \n . Example: VNS_MADS_SEARCH false \n  \n . Default: false\n\n",  "  advanced global mads search vns neighborhood "  , "true" , "true" , "true" },
{ "VNS_MADS_SEARCH_TRIGGER",  "NOMAD::Double",  "0.75",  " VNS Mads search trigger",  " \n  \n . The VNS trigger is the maximum desired ratio of VNS blackbox evaluations \n   over the total number of blackbox evaluations. \n    \n . The VNS search is never executed with a null trigger while a value of 1 \n   allows the search at every iteration \n    \n . If \"VNS_MADS_SEARCH yes\", the default value of 0.75 is taken for the trigger \n  \n . Argument: Double \n  \n . Example: VNS_MADS_SEARCH_TRIGGER 0.9 \n  \n . Default: 0.75\n\n",  "  advanced global mads search vns neighborhood ratio  "  , "true" , "true" , "true" },
{ "VNS_MADS_SEARCH_WITH_SURROGATE",  "bool",  "false",  " VNS Mads search with surrogate",  " \n  \n . VNS search can use static surrogate evaluations for optimization instead of \n   blackbox evaluation. \n    \n . If enabled and a static surrogate (batch mode or library mode) is not \n   available, an exception is triggered. \n    \n . Argument: bool \n  \n . Example: VNS_MADS_SEARCH_WITH_SURROGATE true \n  \n . Default: false\n\n",  "  advanced global mads search vns neighborhood ratio surrogate  "  , "true" , "true" , "true" },
{ "VNS_MADS_SEARCH_MAX_TRIAL_PTS_NFACTOR",  "size_t",  "100",  " VNS-Mads search stopping criterion.",  " \n  \n . VNS Mads stopping criterion. Max number of trial pts < dimension * NFactor \n  \n . Argument: Positive integer. INF disables this criterion. \n  \n . Example: VNS_MADS_SEARCH_MAX_TRIAL_PTS_NFACTOR 10 \n  \n . Default: 100\n\n",  "  advanced global vns neighborhood mads search stop trial  "  , "true" , "true" , "true" } };

#endif
