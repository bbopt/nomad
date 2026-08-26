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

#ifndef __NOMAD_BCQP_SOLVER__
#define __NOMAD_BCQP_SOLVER__

#include "../../../ext/sgtelib/src/Matrix.hpp"

#include "../../nomad_nsbegin.hpp"

enum class BCQPSolverStatus
{
    BOUNDS_ERROR, ///< Problem with lower bounds and upper bounds
    MATRIX_DIMENSIONS_FAILURE, ///< Problem with matrix dimensions
    MAX_ITER_REACHED, ///< Maximum number of iterations reached
    NUM_ERROR, ///< Conjugate gradient numerical error
    TIGHT_VAR_BOUNDS, ///< Bounds on variables are too tight
    STAGNATION_ITERATES, ///< Distance between successive iterates are too low
    SOLVED, ///< Problem solved
    UNDEFINED ///< Undefined status
};

/// BCQP solver
/// Solve the following problem
/// min Q(x) = g0 + g'x + (1/2) x'H x
///  s.t. lb <= x <= ub
class BCQPSolver
{
public:
    BCQPSolverStatus solve(SGTELIB::Matrix &x,
                           const SGTELIB::Matrix &QPModel,
                           const SGTELIB::Matrix &lb,
                           const SGTELIB::Matrix &ub) const;

    // Parameters
    double tol_dist_successive_x;

    size_t max_iter;
    size_t verbose_level;

private:
    static bool checkDimensions(const SGTELIB::Matrix& x,
                                const SGTELIB::Matrix& QPModel,
                                const SGTELIB::Matrix& lb,
                                const SGTELIB::Matrix& ub);

    static bool checkBoundsCompatibilities(const SGTELIB::Matrix& lb,
                                           const SGTELIB::Matrix& ub);

    static void projectedGradientStep(SGTELIB::Matrix& x,
                                      SGTELIB::Matrix& gradQ,
                                      std::vector<bool>& activeLowerVars,
                                      std::vector<bool>& activeUpperVars,
                                      const SGTELIB::Matrix& QPModel,
                                      const SGTELIB::Matrix& lb,
                                      const SGTELIB::Matrix& ub,
                                      const double kappa,
                                      const bool verbose);

    static double projectedArmijoLineSearch(SGTELIB::Matrix& x,
                                            SGTELIB::Matrix& gradQ,
                                            const SGTELIB::Matrix& x_start,
                                            const SGTELIB::Matrix& QPModel,
                                            const SGTELIB::Matrix& lb,
                                            const SGTELIB::Matrix& ub,
                                            const SGTELIB::Matrix& d,
                                            const double fvalue_start,
                                            const double slope,
                                            const double t_max);

    static bool conjugateGradientStep(SGTELIB::Matrix& x,
                                      const SGTELIB::Matrix& H,
                                      const SGTELIB::Matrix& g,
                                      const double xi,
                                      const bool verbose);

    // Return the maximum step size allowed and the step size chosen
    // Useful for logging
    static std::pair<double, double> lineSearchStep(SGTELIB::Matrix& x,
                                                    SGTELIB::Matrix& gradQ,
                                                    SGTELIB::Matrix& xp,
                                                    const SGTELIB::Matrix& d,
                                                    const SGTELIB::Matrix& QPModel,
                                                    const SGTELIB::Matrix& lb,
                                                    const SGTELIB::Matrix& ub,
                                                    std::vector<bool>& activeLowerVars,
                                                    std::vector<bool>& activeUpperVars,
                                                    const bool hasNegativeCurvature);

    // Return || x - P[x - grad Q(x)] ||_inf,
    // and P[X] the projection of X on [lb, ub]
    static double computeFirstOrderError(const SGTELIB::Matrix& x,
                                         const SGTELIB::Matrix& gradQ,
                                         const SGTELIB::Matrix& lb,
                                         const SGTELIB::Matrix& ub);

    static void projectOnBounds(SGTELIB::Matrix& x,
                                const SGTELIB::Matrix& lb,
                                const SGTELIB::Matrix& ub);


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_BCQP_SOLVER__