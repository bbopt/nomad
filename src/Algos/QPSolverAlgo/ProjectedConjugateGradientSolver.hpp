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
/**
 \file   ProjectedConjugateGradientSolver.hpp
 \brief  Projected Conjugate Gradient algorithm
 \author Tangi Migot and Ludovic Salomon
 \see    ProjectedConjugateGradientSolver.cpp
 */
#ifndef __NOMAD_4_5_PROJECTED_CONJUGATE_GRADIENT_SOLVER__
#define __NOMAD_4_5_PROJECTED_CONJUGATE_GRADIENT_SOLVER__

#include "../../../ext/sgtelib/src/Matrix.hpp"

#include "../../nomad_nsbegin.hpp"

enum class ProjectedConjugateGradientSolverStatus
{
    MATRIX_DIMENSIONS_FAILURE, ///< Problem with matrix dimensions
    FACTORIZATION_FAILURE, ///< Factorization for solving linear systems has failed
    QUAD_ROOTS_ERROR, ///< Resolution of quadratic roots subproblem has failed
    TR_NUM_ERROR, ///< Trust-region numerical error
    TR_PARAM_ERROR, ///< Trust-region parameter error
    MAX_ITER_REACHED, ///< Maximum number of iterations reached
    BOUNDARY_REACHED, ///< The algorithm stops because boundary of the trust-region is reached
    NEGATIVE_CURVATURE, ///< The G matrix is not symmetric definite-positive
    NO_INIT_SOLUTION, ///< No initial solution satisfying constraints
    SOLVED ///< Solved: usually when minimum tolerance reached.
};

/// Projected conjugate gradient solver
/// Solve:
/// min 1/2 x' G x + c' x
/// s.t. A x = b
///      |x| <= delta
class ProjectedConjugateGradientSolver
{
public:
    static ProjectedConjugateGradientSolverStatus solve(SGTELIB::Matrix& x,
                                                        const SGTELIB::Matrix& G,
                                                        const SGTELIB::Matrix& c,
                                                        const SGTELIB::Matrix& A,
                                                        const SGTELIB::Matrix& b,
                                                        const double delta,
                                                        const bool verbose = false);
};

#include "../../nomad_nsend.hpp"

#endif //__NOMAD_4_5_PROJECTED_CONJUGATE_GRADIENT_SOLVER__
