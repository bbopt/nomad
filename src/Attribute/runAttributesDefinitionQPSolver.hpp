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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONQPSOLVER__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONQPSOLVER__

_definition = {
{ "QP_OPTIMIZATION",  "bool",  "false",  " Quad model stand alone QP optimization for constrained and unconstrained pbs ",  " \n  \n . QP standalone optimization on quadratic models for constrained and \n   unconstrained problems \n  \n . Argument: bool \n  \n . Standalone optimization will deactivate any other optimization strategy. \n  \n . Example: QP_OPTIMIZATION true \n  \n . Default: false\n\n",  "  advanced sgtelib quadratic quad optimization programming qp  "  , "true" , "false" , "true" },
{ "QP_SEARCH",  "bool",  "false",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . A QP search step for Mads that make several trial points from the poll center \n  \n . Argument: bool \n  \n . Example: QP_SEARCH true \n  \n . Default: false\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_SelectAlgo",  "size_t",  "0",  " Select the algorithm for QP solver ",  " \n  \n . Select the algorithm for QPSolver. \n  \n . #0: Augmented lagrangian \n  \n . #1: Interior point method \n  \n . #2: L1 penalty augmented lagrangian \n  \n . Argument: size_t \n  \n . Example: QP_SelectAlgo 0 \n  \n . Default: 0\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_maxIter",  "size_t",  "20",  " QPSolver outter loop iteration limit ",  " \n  \n . QPSolver outter loop iteration limit. \n  \n . Common to algos #0, AugLag and #2, L1AugLag. \n  \n . Argument: size_t \n  \n . Example: maxIter 100 \n  \n . Default: 20\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_tolDistDX",  "NOMAD::Double",  "-1.0",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver smallest progress accepted \n  \n . Argument: Double \n  \n . Example: QP_tolDistDX 1e-15 \n  \n . Default: -1.0\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_absoluteTol",  "NOMAD::Double",  "1e-3",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver absolute tolerance \n  \n . Argument: Double \n  \n . Example: QP_absoluteTol 1e-3 \n  \n . Default: 1e-3\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_tolCond",  "NOMAD::Double",  "1e-15",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver relative tolerance w.r.t to condition number \n  \n . Argument: Double \n  \n . Example: QP_tolCond 1e-15 \n  \n . Default: 1e-15\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_tolMesh",  "NOMAD::Double",  "1.0",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver relative tolerance w.r.t to mesh size \n  \n . Argument: Double \n  \n . Example: QP_tolMesh 1.0 \n  \n . Default: 1.0\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_relativeTol",  "NOMAD::Double",  "1e-3",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver relative tolerance \n  \n . Argument: Double \n  \n . Example: QP_relativeTol 1e-3 \n  \n . Default: 1e-3\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_verbose",  "bool",  "false",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver verbosity \n  \n . Argument: bool \n  \n . Example: verbose false \n  \n . Default: false\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_verboseFull",  "bool",  "false",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver verbosity of all subiterations \n  \n . Argument: bool \n  \n . Example: QP_verboseFull false \n  \n . Default: false\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_mu0",  "NOMAD::Double",  "0.5",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver initial penalty parameter value for Augmented Lagrangian solver \n  \n . Argument: Double \n  \n . Example: QP_AugLag_mu0 0.5 \n  \n . Default: 0.5\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_muDecrease",  "NOMAD::Double",  "2",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver penalty parameter decrease factor for Augmented Lagrangian solver \n  \n . Argument: Double \n  \n . Example: QP_AugLag_muDecrease 2 \n  \n . Default: 2\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_eta0",  "NOMAD::Double",  "1.0",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver initial constraint tolerance value for Augmented Lagrangian solver \n  \n . Argument: Double \n  \n . Example: QP_AugLag_eta0 1.0 \n  \n . Default: 1.0\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_omega0",  "NOMAD::Double",  "1.0",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver initial subproblem's tolerance value for Augmented Lagrangian solver \n  \n . Argument: Double \n  \n . Example: QP_AugLag_omega0 1.0 \n  \n . Default: 1.0\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_maxIterInner",  "size_t",  "50",  " QPSolver inner iteration limit for the subproblem ",  " \n  \n . QPSolver augmented lagrangian inner iteration limit for the subproblem \n  \n . Argument: size_t \n  \n . Example: QP_AugLag_maxIterInner 10 \n  \n . Default: 50\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_tolDistDXInner",  "NOMAD::Double",  "1e-15",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver smallest progress accepted in the subproblem \n  \n . Argument: Double \n  \n . Example: QP_AugLag_tolDistDXInner 1e-15 \n  \n . Default: 1e-15\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_maxSuccessivFail",  "size_t",  "3",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver maximum number of successiv failure for the subproblem \n  \n . Argument: size_t \n  \n . Example: QP_AugLag_maxSuccessivFail 3 \n  \n . Default: 3\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_AugLag_successRatio",  "NOMAD::Double",  "0.99",  " A quad model based search step for Mads using a QP solver ",  " \n  \n . QPSolver ratio of expected decrease, so that f_{k+1} < ratio * f_k \n  \n . Argument: NOMAD::Double (between 0 and 1) \n  \n . Example: QP_AugLag_successRatio 0.99 \n  \n . Default: 0.99\n\n",  "  advanced algorithm search quadratic quad model programming qp  "  , "true" , "true" , "true" },
{ "QP_SEARCH_MODEL_BOX_SIZE_LIMIT",  "NOMAD::Double",  "0",  " QP solver generates trial points if bounds box size is above limit  ",  " \n  \n . QP solver generates trial points if bounds box size is above limit. \n  \n . Default: 0, always goes into trial point generation step. \n  \n . Default: 0\n\n",  "  advanced algorithm search quadratic quad model programming qp box  "  , "true" , "true" , "true" } };

#endif
