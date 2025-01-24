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
#ifndef __NOMAD_4_5_LEVENBERG_MARQUARDT_SOLVER__
#define __NOMAD_4_5_LEVENBERG_MARQUARDT_SOLVER__

#include "../../../ext/sgtelib/src/Matrix.hpp"

#include "../../nomad_nsbegin.hpp"

enum class LMSolverStatus
{
    BOUNDS_ERROR, ///< Problem with lower bounds and upper bounds
    MATRIX_DIMENSIONS_FAILURE, ///< Problem with matrix dimensions
    MAX_ITER_REACHED, ///< Maximum number of iterations reached
    FAILURE, ///< Levenberg-Marquardt algorithm has failed
    STRICT_PT_FAILURE, ///< Computation of a strict solution has failed
    NUM_ERROR, ///< Trust-region numerical error
    PARAM_ERROR, ///< Parameter error
    TIGHT_VAR_BOUNDS, ///< Bounds on variables are too tight
    STAGNATION_ITERATES, ///< Distance between successive iterates are too low
    IMPROVED, ///< Has improved the solution
    SOLVED ///< Problem solved
};

/// Levenberg-Marquardt solver.
/// Solve:
///  min   || c(x) + s ||^2
/// (x, s)
/// s.t. lb <= x <= ub
///      s >= 0
/// NB: (x, s) is encoded as a single COLUMN vector
class LevenbergMarquardtSolver
{
public:
    LMSolverStatus solve(SGTELIB::Matrix& x,
                         SGTELIB::Matrix& XS,
                         const SGTELIB::Matrix& QPModel,
                         const SGTELIB::Matrix& lvar,
                         const SGTELIB::Matrix& uvar,
                         SGTELIB::Matrix& cX) const;

    // Parameters
    // Tolerances
    double feasibility_tol; ///< When || c(x) + s || below feasibility_tol, it is solved
    double tol; ///< The tolerance
    double tol_dist_successive_x;

    size_t max_iter; ///< The number of maximum iterations
    bool sol_be_strict; ///< Ensures the point is strictly feasible

    size_t verbose_level;

private:

    static bool checkDimensions(const SGTELIB::Matrix& x,
                                const SGTELIB::Matrix& XS,
                                const SGTELIB::Matrix& QPModel,
                                const SGTELIB::Matrix& lvar,
                                const SGTELIB::Matrix& uvar,
                                const SGTELIB::Matrix& cX);

    static bool checkBoundsCompatibilities(const SGTELIB::Matrix& lvar,
                                           const SGTELIB::Matrix& uvar,
                                           const int n);

    static bool checkBoundsTightness(const SGTELIB::Matrix& lvar,
                                     const SGTELIB::Matrix& uvar,
                                     const int n);

    bool checkStartingPointInBounds(const SGTELIB::Matrix& XS,
                                    const SGTELIB::Matrix& lvar,
                                    const SGTELIB::Matrix& uvar,
                                    const int n) const;

};

#include "../../nomad_nsend.hpp"

#endif //__NOMAD_4_5_LEVENBERG_MARQUARDT_SOLVER__
