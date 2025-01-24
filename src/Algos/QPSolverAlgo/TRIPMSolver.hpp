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
#ifndef __NOMAD_4_5_TRIPM_SOLVER__
#define __NOMAD_4_5_TRIPM_SOLVER__

#include "../../../ext/sgtelib/src/Matrix.hpp"

#include "../../nomad_nsbegin.hpp"

enum class TRIPMSolverStatus
{
    BOUNDS_ERROR, ///< Problem with lower bounds and upper bounds
    MATRIX_DIMENSIONS_FAILURE, ///< Problem with matrix dimensions
    MAX_ITER_REACHED, ///< Maximum number of iterations reached
    LM_FAILURE, ///< Levenberg-Marquardt algorithm has failed
    STRICT_PT_FAILURE, ///< Computation of a strict solution has failed
    NUM_ERROR, ///< Trust-region numerical error
    PARAM_ERROR, ///< Parameter error
    TIGHT_VAR_BOUNDS, ///< Bounds on variables are too tight
    STAGNATION_ITERATES, ///< Distance between successive iterates are too low
    SOLVED ///< Problem solved
};

/// Trust-Region interior-point method solver.
/// Solve the QCQP.
class TRIPMSolver
{
public:
    TRIPMSolverStatus solve(SGTELIB::Matrix& x,
                            const SGTELIB::Matrix& QPModel,
                            const SGTELIB::Matrix& lb,
                            const SGTELIB::Matrix& ub) const;

    // Parameters
    double mu_init;
    double mu_decrease;
    double tol_dist_successive_x;

    size_t max_iter_outer;
    size_t max_iter_inner;

    size_t verbose_level;

private:
    TRIPMSolverStatus solveReducedPb(SGTELIB::Matrix& x,
                                     const SGTELIB::Matrix& QPModel,
                                     const SGTELIB::Matrix& lb,
                                     const SGTELIB::Matrix& ub) const;

    static bool checkDimensions(const SGTELIB::Matrix& x,
                                const SGTELIB::Matrix& QPModel,
                                const SGTELIB::Matrix& lb,
                                const SGTELIB::Matrix& ub);

    static bool checkBoundsCompatibilities(const SGTELIB::Matrix& lb,
                                           const SGTELIB::Matrix& ub);

    bool checkParams() const;

    static bool computeStrictlyFeasiblePoint(SGTELIB::Matrix& x,
                                             SGTELIB::Matrix& cons,
                                             const SGTELIB::Matrix& QPModel,
                                             const SGTELIB::Matrix& lvar,
                                             const SGTELIB::Matrix& uvar);

    struct TRIPMErrorMetric
    {
        double projlagGradNorm; // || x - P[x - grad L(x)] ||_inf , where P[x] is the projection of x on [lvar, uvar]
        double projObjGrad; // || x - P[x - grad f(x)] ||_inf, where P[x] is the projection of x on [lvar, uvar]
        double slackLambdaMuNorm; // || - S y - mu e ||_inf
        double cxNorm; // || max(0, c(x)) ||_inf
        double cxInitNorm; // || max(0, c(x0)) ||_inf
    };

    static void computeErrorFunctionMetric(TRIPMErrorMetric& errMetric,
                                           const SGTELIB::Matrix& XS,
                                           const SGTELIB::Matrix& QPModel,
                                           const SGTELIB::Matrix& lvar,
                                           const SGTELIB::Matrix& uvar,
                                           const SGTELIB::Matrix& lambda,
                                           const double mu,
                                           const bool isBarrierProblem);

    static bool computeSlackMultipliers(SGTELIB::Matrix& slackMultipliers,
                                        const SGTELIB::Matrix& XS,
                                        const SGTELIB::Matrix& Jx,
                                        const SGTELIB::Matrix& Gx,
                                        const double mu);

    enum class BarrierSolverStatus
    {
        NUM_ERROR, ///< Numerical error
        FAILURE, ///< Has failed
        SOLVED, ///< Solved
        STAGNATION_ITERATES, ///< Solver has stagnated
        ONE_STEP_MADE, ///< At least one step has been made
        UNDEFINED ///< Undefined (used for initialization)
    };

    BarrierSolverStatus solveBarrierSubproblem(SGTELIB::Matrix& x,
                                               SGTELIB::Matrix& XSp,
                                               SGTELIB::Matrix& cslack,
                                               SGTELIB::Matrix& lambda,
                                               SGTELIB::Matrix& Gx,
                                               SGTELIB::Matrix& cons,
                                               SGTELIB::Matrix& Jx,
                                               TRIPMErrorMetric& errMetric,
                                               const SGTELIB::Matrix& QPModel,
                                               const SGTELIB::Matrix& XS,
                                               const SGTELIB::Matrix& lvar,
                                               const SGTELIB::Matrix& uvar,
                                               const double mu,
                                               const double atolOpt,
                                               const double atolFeas,
                                               const int verbose_degree) const;

    static double computeMeritFctBarrier(const SGTELIB::Matrix& QPModel,
                                         const SGTELIB::Matrix& lvar,
                                         const SGTELIB::Matrix& uvar,
                                         const SGTELIB::Matrix& x,
                                         const SGTELIB::Matrix& XS,
                                         const double mu,
                                         const double nu);

    // Compute Second Order Correction Step y.
    // This procedure is adapted from procedure SOC, described in:
    //
    // "An interior point algorithm for large-scale nonlinear programming"
    // by R.H. Byrd, M.E. Hribar, and J. Nocedal
    // SIAM Journal on Optimization, 9(4), 877-900 (1999)
    // https://doi.org/10.1137/S1052623497325107
    //
    static bool computeSecondOrderCorrectionStep(SGTELIB::Matrix& y,
                                                 const SGTELIB::Matrix& XS,
                                                 const SGTELIB::Matrix& cslackXcan,
                                                 const SGTELIB::Matrix& Jx,
                                                 const SGTELIB::Matrix& pxs,
                                                 const SGTELIB::Matrix& vxs);
};

#include "../../nomad_nsend.hpp"

#endif //__NOMAD_4_5_TRIPM_SOLVER__
