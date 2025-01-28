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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONPSDSSD__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONPSDSSD__

_definition = {
{ "PSD_MADS_OPTIMIZATION",  "bool",  "0",  " PSD-MADS optimization algorithm ",  " \n  \n . Use PSD-MADS algorithm. \n  \n . Argument: bool \n  \n . Description: Parallel Space Decomposition with Mads (no parallelism) \n  \n . This option deactivates any other optimization strategy. \n  \n . Example: PSD_MADS_OPTIMIZATION true \n  \n . Default: 0\n\n",  "  advanced psd mads parallel decomposition  "  , "true" , "false" , "true" },
{ "PSD_MADS_NB_VAR_IN_SUBPROBLEM",  "size_t",  "2",  " Number of variables in PSD-MADS subproblems ",  " \n  \n . When using Parallel Space Decomposition (PSD) MADS algorithm, select the \n   number of variables in Mads subproblems. \n  \n . Argument: a positive integer < INF. \n  \n . Description: Size of subproblems in PSD-Mads. \n  \n . This attribute is used only when PSD-Mads optimization is active. \n  \n . Example: PSD_MADS_NB_VAR_IN_SUBPROBLEM 3 \n  \n . Default: 2\n\n",  "  advanced psd mads parallel decomposition subproblem  "  , "true" , "false" , "true" },
{ "PSD_MADS_NB_SUBPROBLEM",  "size_t",  "INF",  " Number of PSD-MADS subproblems ",  " \n  \n . When using Parallel Space Decomposition (PSD) MADS algorithm, select the number of \n   Mads subproblems. By default (INF), the number of subproblems is adjusted to \n   cover all variables after one pass. \n  \n . Argument: a positive integer. \n  \n . This attribute is used only when PSD-Mads optimization is active. \n  \n . Example: PSD_MADS_NB_SUBPROBLEM 2 \n  \n . Default: INF\n\n",  "  advanced psd mads parallel decomposition subproblem  "  , "true" , "false" , "true" },
{ "PSD_MADS_ITER_OPPORTUNISTIC",  "bool",  "true",  " Opportunistic strategy between the Mads subproblems in PSD-MADS ",  " \n  \n . When using Parallel Space Decomposition (PSD) MADS algorithm, the launch \n   of Mads subproblems during an iteration can be opportunistically stopped when \n   a success is obtained by a Mads subproblem. \n  \n . Argument: bool \n  \n . This attribute is used only when PSD-Mads optimization is active. \n  \n . Example: PSD_MADS_OPPORTUNISTIC false \n  \n . Default: true\n\n",  "  advanced parallel space mads parallel decomposition subproblem opportunistic  "  , "true" , "false" , "true" },
{ "PSD_MADS_ORIGINAL",  "bool",  "false",  " Use NOMAD 3 strategy for mesh update in PSD-MADS ",  " \n  \n . When using Parallel Space Decomposition (PSD) MADS algorithm, \n   NOMAD 3 strategy is to always update the mesh whenever a new pollster is launched. \n   NOMAD 4 strategy is more defined as for which conditions must be met for \n   the mesh to be updated. \n  \n . Argument: bool \n  \n . This attribute is used only when PSD-Mads optimization is active. \n  \n . Example: PSD_MADS_ORIGINAL false \n  \n . Default: false\n\n",  "  advanced parallel space mads parallel decomposition subproblem original  "  , "true" , "false" , "true" },
{ "PSD_MADS_SUBPROBLEM_PERCENT_COVERAGE",  "NOMAD::Double",  "70",  " Percentage of variables that must be covered in subproblems before updating mesh ",  " \n  \n . When using Parallel Space Decomposition (PSD) MADS algorithm, \n   update (enlarge or refine) the mesh when this percentage of variables is \n   covered by subproblems. \n  \n . A lower value makes for more frequent updates. A larger value makes \n   mesh updates less frequent. \n  \n . Argument: Double between 0 and 100 \n  \n . This attribute is used only when PSD-Mads optimization is active. \n  \n . Example: PSD_MADS_SUBPROBLEM_PERCENT_COVERAGE 80 \n  \n . Default: 70\n\n",  "  advanced parallel space mads parallel subproblem  "  , "true" , "false" , "true" } };

#endif
