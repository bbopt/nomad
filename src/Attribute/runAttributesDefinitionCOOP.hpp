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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONCOOP__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONCOOP__

_definition = {
{ "COOP_MADS_OPTIMIZATION",  "bool",  "false",  " COOP-MADS optimization algorithm ",  " \n  \n . Use COOP-MADS algorithm. Only available if code is built with OpenMP enabled. \n  \n . Argument: bool \n  \n . Description: Parallel concurrent Mads \n  \n . This option deactivates any other optimization strategy. \n  \n . Example: COOP_MADS_OPTIMIZATION true \n  \n . Default: false\n\n",  "  advanced coop mads parallel  "  , "true" , "false" , "true" },
{ "COOP_MADS_NB_PROBLEM",  "size_t",  "4",  " Number of COOP-MADS problems ",  " \n  \n . When using COOP MADS optimization, select the number of \n   Mads problems ran in parallel. \n    \n . In addition each Mads algorithm can perform evaluations in parallel when \n   the NB_THREADS_PARALLEL_EVAL is greater than 1. \n  \n . Argument: a positive integer. \n  \n . This attribute is used only when COOP-Mads optimization is active. \n  \n . Example: COOP_MADS_NB_PROBLEM 2 \n  \n . Default: 4\n\n",  "  advanced psd mads parallel decomposition subproblem  "  , "true" , "false" , "true" },
{ "COOP_MADS_OPTIMIZATION_CACHE_SEARCH",  "bool",  "true",  " COOP-MADS cache search for incumbent synchronization ",  " \n  \n . Perform a cache search to update the best incumbent obtained by all Mads. \n  This allows to synchronize the best solutions of the parallel Mads instances. \n  \n . This attribute is used only when COOP-Mads optimization is active. \n  \n . Argument: bool. \n  \n . Example: COOP_MADS_OPTIMIZATION_CACHE_SEARCH false \n  \n . Default: true\n\n",  "  advanced psd mads parallel decomposition subproblem  "  , "true" , "false" , "true" } };

#endif
