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

#include "BCQPSolver.hpp"

#include "../../Util/utils.hpp"
#include "QPModelUtils.hpp"

NOMAD::BCQPSolverStatus NOMAD::BCQPSolver::solve(SGTELIB::Matrix& x,
                                                 const SGTELIB::Matrix& QPModel,
                                                 const SGTELIB::Matrix& lb,
                                                 const SGTELIB::Matrix& ub) const
{
    if (!checkDimensions(x, QPModel, lb, ub)) {
        return NOMAD::BCQPSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }

    if (!checkBoundsCompatibilities(lb, ub)) {
        return NOMAD::BCQPSolverStatus::BOUNDS_ERROR;
    }

    const int n = x.get_nb_rows();

    // Particular case : when the gap between lb and ub is too small, we immediately leave the procedure.
    // NB: A scaling could be applied to prevent this behavior. But at extremely low tolerance, this is not
    // really useful. The solver could have difficulty to move in the decision space.
    for (int i = 0; i < n; ++i)
    {
        if (std::fabs(lb.get(i, 0) - ub.get(i, 0)) <= 1e-8)
        {
            return NOMAD::BCQPSolverStatus::TIGHT_VAR_BOUNDS;
        }
    }

    // Compute a feasible point
    projectOnBounds(x, lb , ub);

    // Allocate memory for the following vectors
    SGTELIB::Matrix gradQ("gradQ", n, 1);
    SGTELIB::Matrix H("H", n, n);
    SGTELIB::Matrix xm("xm", n, 1); // previous iterate
    SGTELIB::Matrix xp("xp", n, 1); // next iterate
    SGTELIB::Matrix d("d", n, 1); // descent direction computed by conjugate gradient
    d.fill(0);

    std::vector<bool> activeLowerVars(n);
    std::vector<bool> previousActiveLowerVars(n);
    std::vector<bool> activeUpperVars(n);
    std::vector<bool> previousActiveUpperVars(n);
    std::vector<bool> activeVars(n);
    std::vector<bool> bindingVars(n);

    // Algorithm parameters
    constexpr double kappa = 0.1;
    constexpr double tau = 1e-7;

    // Utility procedures
    auto computeBindingVarsSet = [](std::vector<bool>& bindingVars,
                                    const SGTELIB::Matrix& gradQ,
                                    const std::vector<bool>& activeLowerVars,
                                    const std::vector<bool>& activeUpperVars)
    {
        const int n = gradQ.get_nb_rows();
        for (int i = 0; i < n; ++i)
        {
            const double gradQi = gradQ.get(i, 0);
            bindingVars[i] = (activeLowerVars[i] && (gradQi >= 0)) || (activeUpperVars[i] && (gradQi <= 0));
        }
    };

    // Compute gradient of the initial point
    QPModelUtils::getModelObjGrad(gradQ, QPModel, x);
    const double ng0 = gradQ.norm_inf();

    // Initialize some parameters for detection algorithm's stagnation
    double distXLoop = NOMAD::INF; // distance between two successive iterates
    double distQQm = NOMAD::INF; // absolute difference between Q(x) and Q(xm)
    bool areActiveVarsSetsChanged = true; // flag to detect if active set of variables has changed

    const bool verbose = verbose_level > 0;
    if (verbose)
    {
        std::printf("\nBound-constrained quadratic algorithm\n");
        std::printf("Number of total variables: %d\n", n);
        std::printf("Stopping criterion tolerance for optimality: %e\n", tau * std::max(1.0, ng0));
        std::printf("Stopping tolerance for projected gradient step: %e\n\n", kappa);

        std::printf("%4s %10s %32s %8s %18s %10s %12s %8s %8s %12s\n",
                    "iter", "f(x)", "|| x - P(x - grad f(x)) ||_inf", "|d|", "max stepsize", "stepsize",
                    "|working|", "|active L|", "|active U|", "|| x - xp ||_inf");
        constexpr int maxLineWidth = 140;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    auto status = NOMAD::BCQPSolverStatus::MAX_ITER_REACHED;
    double a_max = -1, ak = -1;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        double nQproj = computeFirstOrderError(x, gradQ, lb, ub);

        // Compute the set of active variables
        for (int i = 0; i < n; ++i)
        { // We will correct if needed (see for example SolverTools.jl)
            if (std::fabs(x.get(i, 0) - lb.get(i, 0)) <= 1E-15)
                x.set(i, 0, lb[i]);
            if (std::fabs(x.get(i, 0) - ub.get(i, 0)) <= 1E-15)
                x.set(i, 0, ub[i]);
            activeLowerVars[i] = x.get(i, 0) == lb.get(i, 0);
            activeUpperVars[i] = x.get(i, 0) == ub.get(i, 0);
            activeVars[i] = activeLowerVars[i] || activeUpperVars[i];
        }

        if (verbose)
        {
            std::printf("%-4d", iter);
            std::printf(" %+9e", QPModelUtils::getModelObj(QPModel, x));
            std::printf(" %20e", nQproj);
            std::printf(" %22e", d.norm());
            if (a_max < 0)
            {
                std::printf(" %8s", "-");
            }
            else
            {
                std::printf(" %10e", a_max);
            }
            if (ak < 0)
            {
                std::printf(" %12s    ", "-");
            }
            else
            {
                std::printf(" %12e", ak);
            }
            std::printf(" %7zu", std::count(activeVars.begin(), activeVars.end(), false));
            std::printf(" %8zu", std::count(activeLowerVars.begin(), activeLowerVars.end(), true));
            std::printf(" %10zu", std::count(activeUpperVars.begin(), activeUpperVars.end(), true));
            std::printf(" %20e\n", distXLoop);
        }

        // Stopping criteria
        if (nQproj <= tau * std::max(1.0, ng0))
        {
            status = NOMAD::BCQPSolverStatus::SOLVED;
            break;
        }

        // When the set of free variables is empty, stop the algorithm.
        if (std::count(activeVars.begin(), activeVars.end(), true) == n)
        {
            status = NOMAD::BCQPSolverStatus::SOLVED;
            break;
        }

        // When the algorithm is stagnating, stop
        if (!areActiveVarsSetsChanged && ((distXLoop <= tol_dist_successive_x) || (distQQm <= 1e-9)))
        {
            status = NOMAD::BCQPSolverStatus::STAGNATION_ITERATES;
            break;
        }

        // Save current iterate ...
        xm = x;
        const double fxm = QPModelUtils::getModelObj(QPModel, xm);

        // ... and current set of active variables
        previousActiveLowerVars = activeLowerVars;
        previousActiveUpperVars = activeUpperVars;

        // Generate a sequence of projected gradient steps.
        // Update activeLowerVars and activeUpperVars at the same time.
        const bool verbose_projStep = verbose_level > 1;
        projectedGradientStep(x, gradQ, activeLowerVars, activeUpperVars, QPModel, lb, ub, kappa, verbose_projStep);

        // Update the set of active variables
        for (int i = 0; i < n; ++i)
        {
            activeVars[i] = activeLowerVars[i] || activeUpperVars[i];
        }
        QPModelUtils::getModelObjGrad(gradQ, QPModel, x);

        // When the set of working variables is empty, stop the algorithm.
        if (std::count(activeVars.begin(), activeVars.end(), true) == n)
        {
            status = BCQPSolverStatus::SOLVED;
            break;
        }

        // Compute a reduced quadratic model
        QPModelUtils::getModelObjHessian(H, QPModel, x);

        SGTELIB::Matrix gz = [](const SGTELIB::Matrix& gradQ, const std::vector<bool>& activeVars)
        {
            const int n = static_cast<int>(activeVars.size());
            const int nfree = static_cast<int>(std::count(activeVars.begin(), activeVars.end(), false));
            SGTELIB::Matrix gz("gz", nfree, 1);
            int k = 0;
            for (int i = 0; i < n; ++i)
            {
                if (activeVars[i])
                    continue;

                gz.set(k, 0, gradQ.get(i, 0));
                ++k;
            }
            return gz;
        }(gradQ, activeVars);

        SGTELIB::Matrix ZHZ = [](const SGTELIB::Matrix& H, const std::vector<bool>& activeVars)
        {
            const int n = static_cast<int>(activeVars.size());
            const int nfree = static_cast<int>(std::count(activeVars.begin(), activeVars.end(), false));
            SGTELIB::Matrix ZHZ("ZHZ", nfree, nfree);
            int ki = 0;
            for (int i = 0; i < n; ++i)
            {
                if (activeVars[i])
                    continue;

                int kj = 0;
                for (int j = 0; j < n; ++j)
                {
                    if (activeVars[j])
                        continue;

                    ZHZ.set(ki, kj, H.get(i, j));
                    kj += 1;
                }
                ki += 1;
            }
            return ZHZ;
        }(H, activeVars);

        // Generate conjugate gradient iterates to improve the reduced quadratic model
        int nfree = static_cast<int>(std::count(activeVars.begin(), activeVars.end(), false));
        SGTELIB::Matrix dz("dz", nfree, 1);
        dz.fill(0);
        double sufficient_decrease_cg = 1e-3;
        const bool verbose_cg = verbose_level > 1;
        bool hasNegativeCurvature = conjugateGradientStep(dz, ZHZ, gz, sufficient_decrease_cg, verbose_cg);

        // Project on full space
        d.fill(0);
        int kfree = 0;
        for (int i = 0; i < n; ++i)
        {
            if (activeVars[i])
                continue;

            d[i] = dz[kfree];
            kfree++;
        }

        if (d.has_nan())
        {
            status = BCQPSolverStatus::NUM_ERROR;
            break;
        }

        // Execute a line search along d and update the active set of variables
        auto linesearch_params = lineSearchStep(x, gradQ, xp, d, QPModel, lb, ub,
                                                activeLowerVars, activeUpperVars, hasNegativeCurvature);
        a_max = linesearch_params.first;
        ak = linesearch_params.second;

        // Update set of active vars
        for (int i = 0; i < n; ++i)
        {
            activeVars[i] = activeLowerVars[i] || activeUpperVars[i];
        }

        computeBindingVarsSet(bindingVars, gradQ, activeLowerVars, activeUpperVars);
        if (bindingVars != activeVars)
        {
            // Update stagnation's criteria
            distXLoop = SGTELIB::Matrix::distNorm2(x, xm);
            distQQm = std::fabs(fxm - QPModelUtils::getModelObj(QPModel, x));
            areActiveVarsSetsChanged = (previousActiveLowerVars == activeLowerVars) &&
                                       (previousActiveUpperVars == activeUpperVars);
            continue;
        }

        // Reapply the conjugate gradient method with warm-start from dz with tightening
        // stopping criterion
        sufficient_decrease_cg = 1e-5;
        hasNegativeCurvature = conjugateGradientStep(dz, ZHZ, gz, sufficient_decrease_cg, verbose_cg);

        // Project on full space
        d.fill(0);
        kfree = 0;
        for (int i = 0; i < n; ++i)
        {
            if (activeVars[i])
                continue;

            d[i] = dz[kfree];
            kfree++;
        }

        if (d.has_nan())
        {
            status = BCQPSolverStatus::NUM_ERROR;
            break;
        }

        // Execute a line search along d and update the active set of variables
        linesearch_params = lineSearchStep(x, gradQ, xp, d, QPModel, lb, ub,
                                           activeLowerVars, activeUpperVars, hasNegativeCurvature);
        a_max = linesearch_params.first;
        ak = linesearch_params.second;

        // Update stagnation's criteria
        distXLoop = SGTELIB::Matrix::distNorm2(x, xm);
        distQQm = std::fabs(fxm - QPModelUtils::getModelObj(QPModel, x));
        areActiveVarsSetsChanged = (previousActiveLowerVars != activeLowerVars) ||
                                   (previousActiveUpperVars != activeUpperVars);
    }

    if (verbose)
    {
        std::printf("\nStatus: ");
        std::printf("f(x*) = %e\n", QPModelUtils::getModelObj(QPModel, x));
        if (status == BCQPSolverStatus::SOLVED)
        {
            std::printf("Has reached the minimum tolerance:\n");
            std::printf("|| x - P[x - grad f(x)] ||_inf = %e <= 1e-7 * max(1, grad f(x0)) = %e or\n",
                        computeFirstOrderError(x, gradQ, lb, ub), 1e-7 * std::max(1.0, ng0));
            std::printf("no more working variables: |W| = %zu, |active L| = %zu, |active U| = %zu\n",
                        std::count(activeVars.begin(), activeVars.end(), false),
                        std::count(activeLowerVars.begin(), activeLowerVars.end(), true),
                        std::count(activeUpperVars.begin(), activeUpperVars.end(), true));
        }
        else if (status == BCQPSolverStatus::STAGNATION_ITERATES)
        {
            std::printf("The algorithm has stagnated:\n");
            std::printf("Active sets of lower and upper variables have not changed between successive iterations and\n");
            std::printf("|| x - xp || = %e <= %e or\n", distXLoop, tol_dist_successive_x);
            std::printf("|| f(x) - f(xp) || = %e <= %e\n", distQQm, 1e9);
        }
        else if (status == BCQPSolverStatus::MAX_ITER_REACHED)
        {
            std::printf("The maximum number of iterations has been reached\n");
        }
        else
        {
            std::printf("Unknown stopping criterion\n");
        }
    }

    return status;
}

void NOMAD::BCQPSolver::projectedGradientStep(SGTELIB::Matrix& x,
                                              SGTELIB::Matrix& gradQ,
                                              std::vector<bool>& activeLowerVars,
                                              std::vector<bool>& activeUpperVars,
                                              const SGTELIB::Matrix& QPModel,
                                              const SGTELIB::Matrix& lb,
                                              const SGTELIB::Matrix& ub,
                                              const double kappa,
                                              const bool verbose)
{
    const int n = x.get_nb_rows();

    // Allocate memory for some matrices
    SGTELIB::Matrix armijo_x("armijo_x", n, 1);
    SGTELIB::Matrix armijo_gradQ("armijo_gradQ", n, 1);

    // Compute active sets of variables
    std::vector<bool> activeVars(n);
    for (int i = 0; i < n; ++i)
    {
        activeVars[i] = activeLowerVars[i] || activeUpperVars[i];
    }
    std::vector<bool> previousActiveVars(n);
    previousActiveVars = activeVars;

    if (verbose)
    {
        std::printf("\nProjected gradient step\n");
        std::printf("Number of total variables: %d\n\n", n);

        std::printf("%5s %10s %12s %12s %10s\n", "iter", "f(x)", "|| d ||", "step size", "|Active|");
        constexpr int maxLineWidth = 53;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    // Initialize the value between qm and qmp
    double diff_qmqmp = 0;
    double ak = 1;
    double qm=0, qmp=0;
    for (int iter = 0; iter < n; ++iter)
    {
        qm = QPModelUtils::getModelObj(QPModel, x);
        if (verbose)
        {
            std::printf("%-4d %8e %12e %8e %6zu\n", iter, qm, gradQ.norm(), ak,
                        std::count(activeVars.begin(), activeVars.end(), true));
        }

        // Make a descent direction
        gradQ.multiply(-1.0);

        // Apply an Armijo projected line search along -grad Q(x). NB: at the end of the search, ak is >= 0.
        const double slope = - gradQ.normsquare();
        ak = 1; // min(ak, 1) is the initial value for Armijo line search.
        ak = projectedArmijoLineSearch(armijo_x, armijo_gradQ, x, QPModel, lb, ub, gradQ, qm, slope,
                                       ak); // ak starts from the previous value

        // Update iterate: x = P[x + ak dk], where dk = - grad Q(x) and P[X] is the projection of X on [lb, ub]
        gradQ.multiply(ak); x.add(gradQ);
        projectOnBounds(x, lb, ub);

        // Compute active sets
        for (int i = 0; i < n; ++i)
        {
            const double xi = x.get(i, 0);
            const double li = lb.get(i, 0);
            const double ui = ub.get(i, 0);
            activeLowerVars[i] = std::fabs(xi - li) <= 1e-15;
            activeUpperVars[i] = std::fabs(xi - ui) <= 1e-15;
            activeVars[i] = activeUpperVars[i] || activeLowerVars[i];
        }

        // Compute stopping criteria
        const bool haveActiveVarsSetsChanged = previousActiveVars == activeVars;
        if (haveActiveVarsSetsChanged)
            break;

        qmp = QPModelUtils::getModelObj(QPModel, x);
        QPModelUtils::getModelObjGrad(gradQ, QPModel, x);
        if ((qm - qmp) <= kappa * diff_qmqmp)
            break;

        previousActiveVars = activeVars;
        diff_qmqmp = std::max(qm - qmp, diff_qmqmp);
    }
    if (verbose)
    {
        std::printf("\nProjected gradient step has stopped because:\n");
        std::printf("* the set of active variables (|A| = %zu) has not changed ",
                    std::count(activeVars.begin(), activeVars.end(), true));
        std::printf("between two successive iterations or\n");
        std::printf("* tolerance has been reached:\n");
        std::printf("f(x) - f(xm) = %e <= kappa max_{1 <= r < iter} (f(x^{r+1} - f(x^r) = %e\n",
                    qm - qmp, kappa * diff_qmqmp);
        std::printf("\n");
    }
}

double NOMAD::BCQPSolver::projectedArmijoLineSearch(SGTELIB::Matrix& x,
                                                    SGTELIB::Matrix& gradQ,
                                                    const SGTELIB::Matrix& x_start,
                                                    const SGTELIB::Matrix& QPModel,
                                                    const SGTELIB::Matrix& lb,
                                                    const SGTELIB::Matrix& ub,
                                                    const SGTELIB::Matrix& d,
                                                    const double fvalue_start,
                                                    const double slope,
                                                    const double t_max)
{
    // A classical algorithm to compute a step-length satisfying Armijo/Wolfe conditions can be found at
    //
    // "Line Search Algorithms with Guaranteed Sufficient Decrease",
    // by J.J. More, D.J. Thuente, ACM Transactions on Mathematical Software, 20 (1994), Issue 3, pp 286–307
    //
    // Here, we use a procedure described in SolverTools.jl, which is an improved Armijo line search.

    // Line search parameters
    constexpr double armijo_tol = 1E-4; // > 0
    constexpr double t_small = 1E-15; // small
    constexpr double t_decrease = 2.5; // > 1
    constexpr double wolfe_tol = 0.9999; // < 1
    constexpr int wolfe_iter_max = 5; // non-negative integer
    constexpr double t_increase = 5; // > 1

    // Initialize step length
    double tk = std::min(1.0, t_max);

    // x = P[x_start + tk d], where P[X] is the projection on [lb, ub]
    x = d; x.multiply(tk); x.add(x_start);
    projectOnBounds(x, lb, ub);

    double fkp = QPModelUtils::getModelObj(QPModel, x);
    QPModelUtils::getModelObjGrad(gradQ, QPModel, x);
    double slope_t = SGTELIB::Matrix::dot(d, gradQ);

    // First try to increase tk as long as (strong) Wolfe conditions are satisfied.
    // NB: we do not enter often in this loop.
    // * Wolfe condition: - d' grad Q(x + t d) <= - wolfe_tol d' grad Q(x)
    // * Armijo sufficient decrease conditions: f(x + t d) <= f(x) + armijo_tol * t * d' grad Q(x)
    bool satisfies_wolfe_condition = slope_t < wolfe_tol * slope;
    bool satisfies_armijo_condition = fkp <= fvalue_start - armijo_tol * tk * std::fabs(slope);
    int wolfe_iter = 0;
    while (satisfies_wolfe_condition &&  satisfies_armijo_condition &&
            (wolfe_iter < wolfe_iter_max) && (tk <= t_max))
    {
        tk *= t_increase;

        // x = P[x + tk d], where P[X] is the projection of X on [lb, ub]
        x = d; x.multiply(tk); x.add(x_start);
        projectOnBounds(x, lb, ub);

        // Recompute Wolfe and Armijo conditions.
        QPModelUtils::getModelObjGrad(gradQ, QPModel, x);
        fkp = QPModelUtils::getModelObj(QPModel, x);
        slope_t = SGTELIB::Matrix::dot(d, gradQ);
        satisfies_wolfe_condition = slope_t < wolfe_tol * slope;
        satisfies_armijo_condition = fkp <= fvalue_start - armijo_tol * tk * std::fabs(slope);

        ++wolfe_iter;
    }

    // Then, try to satisfy Armijo's condition
    int armijo_iter = 0;
    satisfies_armijo_condition = fkp <= fvalue_start - armijo_tol * tk * std::fabs(slope);
    // Enrich Armijo condition with these numerical tricks, taken from
    // "A new conjugate gradient method with guaranteed descent and an efficient line search",
    // by W. W. HAGER AND H. ZHANG, SIAM Journal on Optimization, 16 (2005), pp. 170–192.
    // armijo = armijo || ((fkp <= fk + 4E-10 * abs(fk)) && (slope_t <= fact * slope));
    while (!satisfies_armijo_condition && tk > t_small) // TODO: add this stopping condition: && (nbA < bA_max))
    {
        tk /= t_decrease;

        // Update iterate: x = P[x + tk d], where P[X] is the projection of X on [lb, ub]
        x = d; x.multiply(tk); x.add(x_start);
        projectOnBounds(x, lb, ub);

        // Recompute Armijo's condition
        fkp = QPModelUtils::getModelObj(QPModel, x);
        satisfies_armijo_condition = fkp <= fvalue_start - armijo_tol * tk * std::fabs(slope);
        armijo_iter++;
    }

    if (!satisfies_armijo_condition)
    {
        return 0.0;
    }

    return tk;
}


bool NOMAD::BCQPSolver::conjugateGradientStep(SGTELIB::Matrix& x,
                                              const SGTELIB::Matrix& H,
                                              const SGTELIB::Matrix& g,
                                              const double xi,
                                              const bool verbose)
{
    const int n = x.get_nb_rows();

    // Pre-allocation of initial vectors
    SGTELIB::Matrix v("v", n, 1);
    for (int i = 0; i < n; i++)
    {
        v[i] = x[i];
    }

    // Pre-allocation of matrices used in the conjugate gradient algorithm
    SGTELIB::Matrix Ask("Ask", n, 1); // The matrix vector product A sk
    SGTELIB::Matrix rk("rk", n, 1);   // The current residual
    SGTELIB::Matrix sk("sk", n, 1);   // The conjugate gradient direction
    SGTELIB::Matrix Avk("Avk", n, 1); // The matrix vector product A vk
    SGTELIB::Matrix::inplace_product(Avk, H, v);

    // Initialization
    // Compute initial residual vector r0 := gW + HW v0
    rk = g;
    rk.add(Avk);

    // Initialize the conjugate direction s0 := -r0
    sk = rk;
    sk.multiply(-1.0);

    // Initialize gamma := rk^t rk
    double gamma =  SGTELIB::Matrix::dot(rk, rk);
    double sNormSquare = gamma;

    // Define tolerance
    double rNorm = std::sqrt(gamma);
    constexpr double atol = 1e-7;
    constexpr double rtol = 1e-7;
    const double tol = atol + rtol * rNorm;

    // Declare some variables for verbose information
    double sAs = 0;

    // Define stopping criteria
    bool solved = rNorm <= tol;
    bool zeroCurvature = false;

    double diff_qmqmp = 0;
    double qmp = 0.5 * SGTELIB::Matrix::dot(Avk, v) + SGTELIB::Matrix::dot(v, g);
    double qm = qmp;

    // Main iteration
    const int max_iter = 2 * n;
    int iter = 0;
    while (!((iter >= max_iter) || solved || zeroCurvature))
    {
        // Compute sk^t HW sk
        SGTELIB::Matrix::inplace_product(Ask, H, sk);
        sAs = SGTELIB::Matrix::dot(Ask, sk);

        // Negative curvature detection
        if (sAs <= atol * atol * sNormSquare)
        {
            if (std::fabs(sAs) <= atol * sNormSquare)
            {
                zeroCurvature = true;
            }

            // When at iteration 0, set v := gW (the negative gradient)
            // otherwise, the iterate candidate at the previous iteration will be
            // returned
            if (iter == 0)
            {
                v = g;
                v.multiply(-1.0);
            }

            solved = true;
        }

        if (zeroCurvature || solved)
            continue;

        // Compute alpha := rk^t rk / sk^t A sk, where A = HW
        double alpha = gamma / sAs;

        // vk := vk + alpha * sk
        for (int i = 0; i < n; ++i)
        {
            v[i] += alpha * sk[i];
        }

        // rk := rk + alpha * A sk, where A = HW
        for (int i = 0; i < n; ++i)
        {
            rk[i] += alpha * Ask[i];
        }

        double gammap = SGTELIB::Matrix::dot(rk, rk);
        rNorm = std::sqrt(gammap);

        // Compute stopping criterion
        solved = rNorm <= tol;
        if (solved)
        {
            continue;
        }

        // Update parameters
        // beta := rk+1^t rk+1 / rk^t rk
        const double beta = gammap / gamma;

        // sk+1 := -rk+1 + beta * sk
        sk.multiply(beta);
        sk.sub(rk);

        gamma = gammap;

        SGTELIB::Matrix::inplace_product(Avk, H, v);

        // Compute special criterion for sufficient conjugate gradient decrease
        qmp = 0.5 * SGTELIB::Matrix::dot(Avk, v) + SGTELIB::Matrix::dot(v, g);
        solved = solved || ((qm - qmp <= xi * diff_qmqmp));
        diff_qmqmp = std::max(qm - qmp, diff_qmqmp);
        qm = qmp;
        iter += 1;
    }

    verbose && std::cout << "CG tol: " << tol;
    verbose && std::cout << "; CG total niter: " << iter;
    verbose && std::cout << "; CG residual norm: " << rNorm;
    verbose && solved && sAs <= 0 && std::cout << "; non positive curvature detected";
    verbose && std::cout << std::endl;

    for (int i = 0; i < n; i++)
    {
        x[i] = v[i];
    }
    return (solved && sAs <= 0);
}

std::pair<double, double> NOMAD::BCQPSolver::lineSearchStep(SGTELIB::Matrix& x,
                                                            SGTELIB::Matrix& gradQ,
                                                            SGTELIB::Matrix& xp,
                                                            const SGTELIB::Matrix& d,
                                                            const SGTELIB::Matrix& QPModel,
                                                            const SGTELIB::Matrix& lb,
                                                            const SGTELIB::Matrix& ub,
                                                            std::vector<bool>& activeLowerVars,
                                                            std::vector<bool>& activeUpperVars,
                                                            const bool hasNegativeCurvature)
{
    // Compute max step size such that x + a_max d remains in [lb, ub]
    const double a_max = [](const SGTELIB::Matrix& x,
                            const SGTELIB::Matrix& d,
                            const SGTELIB::Matrix& lb,
                            const SGTELIB::Matrix& ub)
    {
        const int n = x.get_nb_rows();

        double t_max = NOMAD::INF;
        for (int i = 0; i < n ; ++i)
        {
            const double di = d.get(i, 0);
            double gamma = INF;
            if (di > 0)
            {
                gamma = (ub.get(i, 0) - x.get(i, 0)) / std::fabs(di);
            }
            else if (di < 0)
            {
                gamma = (x.get(i, 0) - lb.get(i, 0)) / std::fabs(di);
            }
            if (gamma < t_max)
            {
                t_max = gamma;
            }
        }
        return t_max;
    }(x, d, lb, ub);

    double a_k;
    if (hasNegativeCurvature)
    {
        // d is a direction of negative curvature: we find the smallest step gamma such that
        // x + gamma d is at the boundary of Omega^{k,j} = {y : lb <= y <= ub}.
        // We need the objective function to decrease (to compute the predicted reduction in
        // the trust-region test), we then add a test.
        xp = d; xp.multiply(a_max); xp.add(x);
        projectOnBounds(xp, lb, ub);
        if (QPModelUtils::getModelObj(QPModel, xp) <= QPModelUtils::getModelObj(QPModel, x))
        {
            // We accept the step
            a_k = a_max;
        }
        else
        {
            // We perform a backtracking strategy, starting from a_max and decreasing along the line
            const double qm = QPModelUtils::getModelObj(QPModel, x);
            double stepsize = a_max;
            constexpr int max_iter_trials = 10;
            for (int iter_trial = 0; iter_trial < max_iter_trials; ++iter_trial)
            {
                xp = d;
                xp.multiply(stepsize);
                xp.add(x);
                projectOnBounds(xp, lb, ub);
                const double qxp = QPModelUtils::getModelObj(QPModel, xp);
                const double slope = SGTELIB::Matrix::dot(gradQ, xp - x);
                if (qxp <= qm + 1e-4 * slope)
                    break;

                stepsize /= 3.0;
            }
            a_k = QPModelUtils::getModelObj(QPModel, xp) <= QPModelUtils::getModelObj(QPModel, x) ? stepsize : 0;
        }
    }
    else
    {
        // d is not a direction of negative curvature: we execute a projected Armijo line search.
        // 1e-15 is the smallest step allowed.
        const double slope = SGTELIB::Matrix::dot(d, gradQ);
        const double qm = QPModelUtils::getModelObj(QPModel, x);
        a_k = a_max > 1e-15 ? projectedArmijoLineSearch(xp, gradQ, x, QPModel, lb, ub, d, qm, slope, a_max) : 0;
    }

    // Update the incumbent.
    if (a_k > 0)
    {
        const int n = x.get_nb_rows();
        for (int i = 0; i < n; ++i)
        {
            x.set(i, 0, x.get(i, 0) + a_k * d.get(i, 0));
        }
        projectOnBounds(x, lb, ub);

        // Update the set of active constraints
        for (int i = 0; i < n; ++i)
        {
            const double xi = x.get(i, 0);
            const double li = lb.get(i, 0);
            const double ui = ub.get(i, 0);
            activeLowerVars[i] = std::fabs(xi - li) <= 1e-15;
            activeUpperVars[i] = std::fabs(xi - ui) <= 1e-15;
        }
    }
    QPModelUtils::getModelObjGrad(gradQ, QPModel, x);
    return {a_max, a_k};
}

bool NOMAD::BCQPSolver::checkDimensions(const SGTELIB::Matrix& x,
                                        const SGTELIB::Matrix& QPModel,
                                        const SGTELIB::Matrix& lb,
                                        const SGTELIB::Matrix& ub)
{
    const int n = x.get_nb_rows();
    if (n != std::max(x.get_nb_rows(), x.get_nb_cols()) && (x.get_nb_cols() != 1))
    {
        std::string err = "BCQPSolver::solve error: x must be a column vector";
        std::printf("%s\n", err.c_str());
        return false;
    }

    if (n != lb.get_nb_rows() || n != ub.get_nb_rows())
    {
        std::string err = "BCQPSolver::solve error: bound constraints dimensions ";
        err += "nlb = " + std::to_string(lb.get_nb_cols()) + " nub = " + std::to_string(ub.get_nb_cols());
        err += " are not compatible with dimension of x (n = " + std::to_string(n) + ")";
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbParams = QPModel.get_nb_cols();
    if (nbParams != (n + 1) + n * (n + 1) / 2)
    {
        std::string err = "BCQPSolver::solve error: ";
        err += "the number of params of the model nbParams = (n+1) * (n+2) / 2 = " + std::to_string(nbParams);
        err += " is not compatible with the dimension of the solution n = " + std::to_string(n);
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbCons = QPModel.get_nb_rows() - 1;
    if (nbCons > 0)
    {
        std::string err = "L1AugLagSolver::solve error: ";
        err += "the model has constraints";
        std::printf("%s\n", err.c_str());
        return false;
    }

    return true;
}

bool NOMAD::BCQPSolver::checkBoundsCompatibilities(const SGTELIB::Matrix& lb,
                                                   const SGTELIB::Matrix& ub)
{
    const int n = lb.get_nb_cols();
    for (int i = 0; i < n; ++i)
    {
        const bool areBoundsCompatible = lb.get(i, 0) <= ub.get(i, 0);
        if (!areBoundsCompatible)
        {
            std::string err = "BCQPSolver::solve error: ";
            err += "no compatibility between lower bound and upper bound for index " + std::to_string(i);
            std::printf("%s\n", err.c_str());
            return false;
        }
    }
    return true;
}

double NOMAD::BCQPSolver::computeFirstOrderError(const SGTELIB::Matrix& x,
                                                 const SGTELIB::Matrix& gradQ,
                                                 const SGTELIB::Matrix& lb,
                                                 const SGTELIB::Matrix& ub)
{
    const int n = x.get_nb_rows();
    SGTELIB::Matrix dualFeas("dualFeas", n, 1);
    for (int i = 0; i < n; ++i)
    {
        dualFeas.set(i, 0, x.get(i, 0) - gradQ.get(i, 0));
    }
    projectOnBounds(dualFeas, lb, ub);
    dualFeas.sub(x);

    return dualFeas.norm_inf();
}

void NOMAD::BCQPSolver::projectOnBounds(SGTELIB::Matrix& x,
                                        const SGTELIB::Matrix& lb,
                                        const SGTELIB::Matrix& ub)
{
    const int n = x.get_nb_rows();
    for (int i = 0; i < n; ++i)
    {
        const double xi = std::clamp(x.get(i, 0), lb.get(i, 0), ub.get(i, 0));
        x.set(i, 0, xi);
    }
}
