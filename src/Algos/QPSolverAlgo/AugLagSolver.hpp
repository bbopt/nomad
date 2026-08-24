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

#ifndef __NOMAD_AUGLAG_SOLVER__
#define __NOMAD_AUGLAG_SOLVER__

#include "../../../ext/sgtelib/src/Matrix.hpp"

#include "../../nomad_nsbegin.hpp"

enum class AugLagSolverStatus
{
    BOUNDS_ERROR, ///< Problem with lower bounds and upper bounds
    MATRIX_DIMENSIONS_FAILURE, ///< Problem with matrix dimensions
    MAX_ITER_REACHED, ///< Maximum number of iterations reached
    LM_FAILURE, ///< Levenberg-Marquardt failure
    NUM_ERROR, ///< Conjugate gradient numerical error
    TIGHT_VAR_BOUNDS, ///< Bounds on variables are too tight
    STAGNATION_ITERATES, ///< Distance between successive iterates are too low
    SOLVED, ///< Problem solved
    UNDEFINED ///< Undefined status
};

// Augmented Lagrangian solver
// Solve the QCQP
class AugLagSolver
{
public:
    AugLagSolverStatus solve(SGTELIB::Matrix& x,
                             SGTELIB::Matrix& QPModel,
                             SGTELIB::Matrix& lb,
                             SGTELIB::Matrix& ub) const;

    // Parameters
    double mu_init; // Initial augmented Lagrangian penalty parameter
    double mu_decrease; // Augmented Lagrangian penalty decrease coefficient
    double omega_init; // Initial tolerance parameter to solve the inner subproblem
    double eta_init; // Initial tolerance parameter on the constraints
    double tol_dist_successive_x;

    size_t max_iter_outer;
    size_t max_iter_inner;
    size_t max_iter_bcqp;

    size_t verbose_level;

private:
    bool checkParameters() const;

    static bool checkDimensions(const SGTELIB::Matrix& x,
                                const SGTELIB::Matrix& QPModel,
                                const SGTELIB::Matrix& lb,
                                const SGTELIB::Matrix& ub);

    static bool checkBoundsCompatibilities(const SGTELIB::Matrix& lb,
                                           const SGTELIB::Matrix& ub);

    enum class BoundAugLagSolverStatus
    {
        MAX_ITER_REACHED, ///< Maximum iterations reached
        NUM_ERROR, ///< Numerical error
        FAILURE, ///< Has failed
        SOLVED, ///< Solved
        STAGNATION_ITERATES, ///< Solver has stagnated
        ONE_STEP_MADE, ///< At least one step has been made
        UNDEFINED
    };

    BoundAugLagSolverStatus solveBoundAugLag(SGTELIB::Matrix& XS,
                                             const SGTELIB::Matrix& QPModel,
                                             const SGTELIB::Matrix& lvar,
                                             const SGTELIB::Matrix& uvar,
                                             const SGTELIB::Matrix& lambda,
                                             const double mu,
                                             const double omega) const;

    static void computeBoundAugLagQPModel(SGTELIB::Matrix& model,
                                          const SGTELIB::Matrix& QPModel,
                                          const SGTELIB::Matrix& XS,
                                          const SGTELIB::Matrix& lambda,
                                          const double mu);

    static double computeAugLagObj(const SGTELIB::Matrix& QPModel,
                                   const SGTELIB::Matrix& XS,
                                   const SGTELIB::Matrix& lambda,
                                   const double mu);

    static void computeAugLagGrad(SGTELIB::Matrix& gLag,
                                  const SGTELIB::Matrix& QPModel,
                                  const SGTELIB::Matrix& XS,
                                  const SGTELIB::Matrix& lambda,
                                  const double mu);

    static void computeAugLagHessian(SGTELIB::Matrix& HLag,
                                     const SGTELIB::Matrix& QPModel,
                                     const SGTELIB::Matrix& XS,
                                     const SGTELIB::Matrix& lambda,
                                     const double mu);

    // Return || x - P[x - grad L(x)] ||_inf, where L is the lagrangian at x
    // and P[X] the projection of X on [lb, ub]
    static double computeFirstOrderError(const SGTELIB::Matrix& x,
                                         const SGTELIB::Matrix& gradL,
                                         const SGTELIB::Matrix& lb,
                                         const SGTELIB::Matrix& ub);

    static void projectOnBounds(SGTELIB::Matrix& x,
                                const SGTELIB::Matrix& lb,
                                const SGTELIB::Matrix& ub);
};

#include "../../nomad_nsend.hpp"

#endif //__NOMAD_AUGLAG_SOLVER__
