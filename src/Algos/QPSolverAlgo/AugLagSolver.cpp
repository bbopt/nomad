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

#include "AugLagSolver.hpp"

#include "BCQPSolver.hpp"
#include "LevenbergMarquardtSolver.hpp"
#include "QPModelUtils.hpp"
#include "../../Util/utils.hpp"

NOMAD::AugLagSolverStatus NOMAD::AugLagSolver::solve(SGTELIB::Matrix& x,
                                                     SGTELIB::Matrix& QPModel,
                                                     SGTELIB::Matrix& lb,
                                                     SGTELIB::Matrix& ub) const
{
    if (!checkDimensions(x, QPModel, lb, ub))
    {
        return NOMAD::AugLagSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }

    if (!checkBoundsCompatibilities(lb, ub))
    {
        return NOMAD::AugLagSolverStatus::BOUNDS_ERROR;
    }

    const int n = x.get_nb_rows();

    // Particular case : when the gap between lb and ub is too small, we immediately leave the procedure.
    // NB: A scaling could be applied to prevent this behavior. But at extremely low tolerance, this is not
    // really useful. The solver could have difficulty to move in the decision space.
    for (int i = 0; i < n; ++i)
    {
        if (std::fabs(lb.get(i, 0) - ub.get(i, 0)) <= 1e-8)
        {
            return NOMAD::AugLagSolverStatus::TIGHT_VAR_BOUNDS;
        }
    }

    projectOnBounds(x, lb, ub);

    const int ncons = QPModel.get_nb_rows() - 1;
    SGTELIB::Matrix cons("cons", ncons, 1);
    QPModelUtils::getModelCons(cons, QPModel, x);

    // Initialize x + s, where s are the slack variables associated to the constraints
    const int nvar = n + ncons;
    SGTELIB::Matrix XS("XS", nvar, 1);
    for (int i = 0; i < n; ++i)
    {
        XS.set(i, 0, x.get(i, 0));
    }
    // s = std::max(0, - c(x))
    for (int i = 0; i < ncons; ++i)
    {
        XS.set(i + n, 0, std::max(0.0, -cons.get(i, 0)));
    }

    // Initialize the lower and upper bounds associated to x+s
    SGTELIB::Matrix lvar("lvar", nvar, 1);
    lvar.fill(0.0);
    for (int i = 0; i < n; ++i)
    {
        lvar.set(i, 0, lb.get(i, 0));
    }
    SGTELIB::Matrix uvar("uvar", nvar, 1);
    uvar.fill(NOMAD::INF);
    for (int i = 0; i < n; ++i)
    {
        uvar.set(i, 0, ub.get(i, 0));
    }

    // Note: it is not mandatory for the algorithm to start from a feasible point. But it is more efficient
    // in practice. It is a well-known heuristics for augmented Lagrangian algorithms.
    NOMAD::LevenbergMarquardtSolver LMalgo = {mu_init, // feasibility_tol
                                              omega_init, // tolerance
                                              1e-15, // tol_dist_successive_x
                                              30, // Max number of iterations
                                              false, // sol_be_strict
                                              verbose_level > 0 ? verbose_level - 1 : 0}; // verbose level
    LMalgo.solve(x, XS, QPModel, lvar, uvar, cons);

    // Check feasibility of the starting point: if not, stop the procedure.
    for (int i = 0; i < nvar; ++i)
    {
        const bool feasible = (XS.get(i, 0) >= lvar.get(i, 0)) && (XS.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            return NOMAD::AugLagSolverStatus::LM_FAILURE;
        }
    }

    // Stopping criterion
    // a- Compute || x0 - P[x - grad f(x0)] ||_inf, where P[X] is the projection of x on [lb, ub]
    SGTELIB::Matrix gradL("gradL", n, 1);
    QPModelUtils::getModelObjGrad(gradL, QPModel, x);
    const double nprojFx0 = computeFirstOrderError(x, gradL, lb, ub);

    // b- Compute || c(x0) + s0 ||_inf
    SGTELIB::Matrix cslack("c+s", ncons, 1);
    for (int j = 0; j < ncons; ++j)
    {
        cslack.set(j, 0, XS.get(j + n, 0) + cons.get(j, 0));
    }
    const double ncx0s0 = cslack.norm_inf();

    // Initialize Lagrange multipliers
    SGTELIB::Matrix lambda_l("lambda_l", ncons, 1);
    lambda_l.fill(0.0);

    // Allocate memory for other vectors
    SGTELIB::Matrix XSp("XSp", nvar, 1);
    SGTELIB::Matrix xp("xp", n, 1);

    // Outer parameters
    double mu_l = mu_init;
    double eta_l = eta_init;
    double omega_l = omega_init;

    // Tolerances
    constexpr double tol_opt = 1e-7;
    constexpr double tol_feas = 1e-7;

    size_t successive_failure = 0;
    size_t successive_acceptable = 0;
    constexpr size_t successive_failure_before_update = 2;
    constexpr size_t max_successive_failure = 5;

    const bool verbose = verbose_level > 0;
    if (verbose)
    {
        std::printf("\nAugmented lagrangian algorithm\n");
        std::printf("Number of total variables: %d\n", n);
        std::printf("Number of inequality constraints: %d\n", ncons);
        std::printf("Stopping criterion tolerance for optimality: %e\n", tol_opt * std::max(1.0, nprojFx0));
        std::printf("Stopping criterion tolerance for feasibility: %e\n", tol_feas * std::max(1.0, ncx0s0));
        std::printf("Maximum number of iterations allowed for outer loop: %zu\n", max_iter_outer);
        std::printf("Maximum number of iterations allowed for inner loop: %zu\n", max_iter_inner);
        std::printf("Maximum of successive failures: %zu\n", max_successive_failure);
        std::printf("Maximum of successive failures before update: %zu\n\n", successive_failure_before_update);

        std::printf("%12s %15s %22s %28s %16s %15s %8s %12s %12s %20s %31s %12s\n",
                    "iter (outer)", "f(x)", "|| max(c(x), 0) ||_inf", "|| x - P[x-grad L(x)] ||_inf",
                    "|| c(x) + s ||_inf", "|| lambda ||_inf", "mu", "omega", "eta", "|| x - xp ||_inf",
                    "successive slow progress iter", "status inner");
        constexpr int maxLineWidth = 220;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    // NB: for the outer loop, we do not consider the distance between successive
    // iterates. For example, if the current iterate has a stopping criterion
    // value far below the threshold omega, the algorithm will update the
    // lagrangian multipliers, even if the iterate remains identical.
    // Even though the inner sub-problem is not solved at optimality
    // in next iterates, we could hope that the algorithm will keep on.
    // This is only kept for logging.
    double distXOuterLoop = NOMAD::INF;
    double auglagValXS = computeAugLagObj(QPModel, XS, lambda_l, mu_l);

    auto status = NOMAD::AugLagSolverStatus::MAX_ITER_REACHED;
    auto status_inner_loop = BoundAugLagSolverStatus::UNDEFINED;
    for (int iter = 0; iter < max_iter_outer; ++iter)
    {
        // Compute || XS - P[XS -grad L(XS)] ||_inf where XS is the projection
        // of XS in [lvar, uvar]
        QPModelUtils::getModelLagrangianGrad(gradL, QPModel, x, lambda_l);
        QPModelUtils::getModelCons(cons, QPModel, x);
        for (int j = 0; j < ncons; ++j)
        {
            cslack.set(j, 0, XS.get(j + n, 0) + cons.get(j, 0));
        }
        double nprojL = computeFirstOrderError(x, gradL, lb, ub);
        // We have computed || x - P[x - grad L(x)] ||_inf
        // To compute || XS - P[XS -grad L(XS)] ||_inf, we need to compute
        // max(|| x - P[x - grad L(x)] ||_inf,
        //     || s - P[s - grad L(s)] ||_inf = || s - P[s + lambda] ||_inf)
        // where P[S] is the projection on [0, +inf]
        for (int j = 0; j < ncons; ++j)
        {
            const double sj = XS.get(j + n, 0);
            const double lambdaj = lambda_l.get(j, 0);
            const double nprojLs = std::fabs(sj - std::clamp(sj + lambdaj, 0.0, NOMAD::INF));
            nprojL = std::max(nprojLs, nprojL);
        }

        if (verbose)
        {
            const double objVal = QPModelUtils::getModelObj(QPModel, x);
            double normCxMax0 = 0.0;
            for (int i = 0; i < ncons; ++i)
            {
                normCxMax0 = std::max(cons.get(i, 0), normCxMax0);
            }

            const std::string statusLog = [](const BoundAugLagSolverStatus status) -> std::string
            {
                if (status == BoundAugLagSolverStatus::UNDEFINED)
                {
                    return "-";
                }
                if (status == BoundAugLagSolverStatus::MAX_ITER_REACHED)
                {
                    return "Max iteration reached";
                }
                if (status == BoundAugLagSolverStatus::SOLVED)
                {
                    return "Solved";
                }
                if (status == BoundAugLagSolverStatus::STAGNATION_ITERATES)
                {
                    return "Stagnation iterates";
                }
                if (status == BoundAugLagSolverStatus::ONE_STEP_MADE)
                {
                    return "Has progressed";
                }

                return "error";
            }(status_inner_loop);

            std::printf(" %-12d %14e %16e %25e %23e %16e %16e %10e %10e %12e %18zu %27s\n",
                        iter, objVal, normCxMax0, nprojL, cslack.norm_inf(),
                        lambda_l.norm(), mu_l, omega_l, eta_l, distXOuterLoop,
                        successive_failure, statusLog.c_str());
        }

        // Check first-order conditions for optimality
        if (nprojL <= tol_opt * std::max(1.0, nprojFx0) &&
            cslack.norm_inf() <= tol_feas * std::max(1.0, ncx0s0))
        {
            status = NOMAD::AugLagSolverStatus::SOLVED;
            break;
        }

        // Detect stagnation
        if (mu_l <= tol_opt / mu_decrease ||
            successive_failure >= max_successive_failure
        )
        {
            status = NOMAD::AugLagSolverStatus::STAGNATION_ITERATES;
            break;
        }

        // Save current iterate
        XSp = XS;

        // Solve inner problem from XSp
        status_inner_loop = solveBoundAugLag(XSp, QPModel, lvar, uvar, lambda_l, mu_l, omega_l);
        if (status_inner_loop == BoundAugLagSolverStatus::NUM_ERROR)
        {
            status = AugLagSolverStatus::NUM_ERROR;
            break;
        }
        if (status_inner_loop == BoundAugLagSolverStatus::STAGNATION_ITERATES)
        {
            status = AugLagSolverStatus::STAGNATION_ITERATES;
            break;
        }

        for (int i = 0; i < n; ++i)
        {
            xp.set(i, 0, XSp.get(i, 0));
        }
        QPModelUtils::getModelCons(cons, QPModel, xp);
        for (int j = 0; j < ncons; ++j)
        {
            cslack.set(j, 0, XSp.get(j + n, 0) + cons.get(j, 0));
        }
        double cxp = cslack.norm_inf();
        const double auglagValXSp = computeAugLagObj(QPModel, XSp, lambda_l, mu_l);

        distXOuterLoop = SGTELIB::Matrix::distNorm2(x, xp);

        // Update parameters
        if (cxp <= eta_l)
        {
            // Update lambda multipliers
            cslack.multiply(- 1 / mu_l);
            lambda_l.add(cslack);

            // Update eta
            eta_l = std::max(eta_l * std::pow(mu_l, 0.9), tol_opt);

            // Update omega
            if (status_inner_loop == BoundAugLagSolverStatus::SOLVED)
            {
                successive_acceptable = 0;
                omega_l *= mu_l;
            }
            else if (successive_acceptable >= successive_failure_before_update)
            {
                // Decrease the tolerance for the subproblem slower
                successive_acceptable += 1;
                omega_l *= std::sqrt(mu_l);
            }
            else
            {
                successive_acceptable += 1; // Try rerun from new point
            }
            omega_l = std::max(omega_l, tol_opt);
        }
        else // not feasible
        {
            // Update mu
            if (status_inner_loop == BoundAugLagSolverStatus::SOLVED)
            {
                successive_acceptable = 0;
                mu_l /= mu_decrease;
            }
            else if (successive_acceptable >= successive_failure_before_update)
            {
                // The subproblem has not been solved at optimality, but we still want to keep on
                // iterating. If mu decreases too fast, the algorithm can fail quick.
                successive_acceptable += 1;
                mu_l /= sqrt(mu_decrease);
            }
            else
            {
                successive_acceptable += 1; // Try rerun from new point
            }

            // Update eta and omega
            eta_l = std::max(std::pow(mu_l, 0.1), tol_opt);
            omega_l = mu_l;
        }

        if (status_inner_loop == BoundAugLagSolverStatus::SOLVED || (auglagValXSp < 0.99 * auglagValXS))
        {
            // Update the current iterate
            XS = XSp;
            x = xp;
            auglagValXS = computeAugLagObj(QPModel, XS, lambda_l, mu_l);
            successive_failure = 0;
            continue;
        }

        // The algorithm is started to stagnate, and we allow it to continue a bit
        // for next iterations.
        successive_failure += 1;

        // Heuristics: try to decrease the infeasibility one iteration before the end
        if (successive_failure < max_successive_failure - 1)
        {
            LMalgo.tol_dist_successive_x = 1e-15;
            LMalgo.feasibility_tol = omega_l;
            LMalgo.max_iter = 30;
            auto feasibility_status = LMalgo.solve(xp, XSp, QPModel, lvar, uvar, cons);
            if (feasibility_status == NOMAD::LMSolverStatus::SOLVED ||
                feasibility_status == NOMAD::LMSolverStatus::IMPROVED)
            {
                distXOuterLoop = SGTELIB::Matrix::distNorm2(XS, XSp);
                XS = XSp;
                x = xp;

                auglagValXS = computeAugLagObj(QPModel, XS, lambda_l, mu_l);
                // If there is a real change, reset it.
                if (distXOuterLoop > tol_dist_successive_x)
                {
                    successive_failure = 0;
                }
            }
            continue;
        }

        // Keep on going
        XS = XSp;
        x = xp;
        auglagValXS = computeAugLagObj(QPModel, XS, lambda_l, mu_l);
    }

    if (verbose)
    {
        std::printf("\nStatus: ");
        std::printf("f(x*) = %e\n", QPModelUtils::getModelObj(QPModel, x));
        double normCxMax0 = 0;
        for (int i = 0; i < ncons; ++i)
        {
            normCxMax0 = std::max(cons.get(i, 0), normCxMax0);
        }
        std::printf("|| max (c(x), 0) ||_inf = %e\n", normCxMax0);

        if (status == AugLagSolverStatus::SOLVED)
        {
            std::printf("Has reached the minimum tolerance\n");
            std::printf("|| x - P[x - grad L(x)] ||_inf = %e <= max(1.0, || x - P[x - grad f(x)] ||) * tol_opt = %e and\n",
                        computeFirstOrderError(x, gradL, lb, ub), std::max(1.0, nprojFx0) * tol_opt);
            std::printf("|| c(x) + s ||_inf = %e <= max(1.0, || c(x0) + s0 ||_inf) * tol_feas = %e\n",
                         cslack.norm_inf(), std::max(1.0, ncx0s0) * tol_feas);
        }
        else if (status == AugLagSolverStatus::STAGNATION_ITERATES)
        {
            std::printf("Outer steps have stagnated:\n");
            std::printf("|| x - xp || = %e <= %e or\n", distXOuterLoop, tol_dist_successive_x);
            std::printf("mu = %e <= tol_opt / mu_decrease = %e or\n", mu_l, tol_opt / mu_decrease);
            std::printf("The inner solver has not made sufficient progress %zu times, superior to the following threshold = %zu\n",
                        successive_failure, max_successive_failure);
        }
        else if (status == AugLagSolverStatus::MAX_ITER_REACHED)
        {
            std::printf("The maximum number of outer iterations has been reached\n");
        }
        else
        {
            std::printf("Unknown stopping criterion\n");
        }
        printf("\n");
    }

    return status;
}

NOMAD::AugLagSolver::BoundAugLagSolverStatus NOMAD::AugLagSolver::solveBoundAugLag(SGTELIB::Matrix& XS,
                                                                                   const SGTELIB::Matrix& QPModel,
                                                                                   const SGTELIB::Matrix& lvar,
                                                                                   const SGTELIB::Matrix& uvar,
                                                                                   const SGTELIB::Matrix& lambda,
                                                                                   const double mu,
                                                                                   const double omega) const
{
    const int nvar = XS.get_nb_rows();
    const int ncons = QPModel.get_nb_rows() - 1;
    const int n = nvar - ncons;

    // Allocate memory for the following matrices
    const int nbParamsBCQPModel = nvar+1 + nvar * (nvar+1) / 2;
    SGTELIB::Matrix BCQP_model("BCQP_model", 1, nbParamsBCQPModel);
    SGTELIB::Matrix XScan("XScan", nvar, 1);
    SGTELIB::Matrix xcan("xcan", n, 1);
    SGTELIB::Matrix gradLa("gradLa", nvar, 1);
    SGTELIB::Matrix gradLaXScan("gradLaXScan", nvar, 1);
    SGTELIB::Matrix dlvar("dlvar", nvar, 1);
    SGTELIB::Matrix duvar("duvar", nvar, 1);
    SGTELIB::Matrix d("d", nvar, 1);
    d.fill(0);

    // Inner loop parameters
    // Trust region update parameters
    constexpr double epsilon_1 = 0.05;
    constexpr double epsilon_2 = 0.9;
    constexpr double gamma_1 = 0.5; // 0.01; // 0.25;
    constexpr double gamma_2 = 2; // 100; // 2.5

    // Trust region radius parameters
    constexpr double min_tr_radius = 1e-15;
    constexpr double max_tr_radius = 1e15;

    // Tolerance for this algorithm
    const double tol_opt = std::min(1e-8, omega);

    // Others
    constexpr size_t max_unsuccessful_iter = 40;
    constexpr double tol_tr_ratio = 1e-15;

    // Compute initial trust-region radius:
    // Delta0 = 0.1 * || XS - P[XS - grad La(XS, lambda0, mu)] ||_inf, where
    // P[X] is the projection of X on [lvar, uvar]
    computeBoundAugLagQPModel(BCQP_model, QPModel, XS, lambda, mu);
    computeAugLagGrad(gradLa, QPModel, XS, lambda, mu);
    const double ngproj0 = computeFirstOrderError(XS, gradLa, lvar, uvar);
    double delta = 0.1 * ngproj0;

    // Compute tolerance success
    const double tol_success = omega;
    //const double tol_success = omega >= 1 ? omega : omega * (1.0 + ngproj0);

    // Initialize parameters for trust-region non-monotone strategy
    constexpr size_t max_iter_non_monotone_steps = 10;
    size_t iter_non_monotone_steps = 0;
    double sig_ref = 0.0, sig_can = 0.0;
    double f_current = computeAugLagObj(QPModel, XS, lambda, mu);
    double f_min = f_current, f_ref = f_current, f_can = f_current;
    double rho_his = 0;

    const bool verbose = verbose_level >= 2;

    // Initial logging
    if (verbose)
    {
        std::printf("\nAugmented Lagrangian bound constrained algorithm\n");
        std::printf("Number of variables: %d\n", nvar);
        std::printf("Number of inequality constraints: %d\n", ncons);
        std::printf("Maximum number of iterations allowed for inner loop: %zu\n", max_iter_inner);
        std::printf("Maximum number of unsuccessful iterations allowed for inner loop: %zu\n", max_unsuccessful_iter);
        std::printf("Initial trust-region radius: %e\n", delta);
        std::printf("Trust-region increase parameter: %e\n", gamma_2);
        std::printf("Trust-region decrease parameter: %e\n", gamma_1);
        std::printf("Trust-region acceptance tolerance: %e\n", epsilon_1);
        std::printf("Trust-region increasing tolerance: %e\n", epsilon_2);
        std::printf("mu parameter: %e\n", mu);
        std::printf("tolerance success: %e\n\n", std::max(tol_success, tol_opt * ngproj0));

        std::printf("%5s %s %9s %10s %46s %15s %10s %12s %12s %16s %12s\n",
                    "iter (inner)", "nb unsuccessful iter", "non monotone iter", "La(x)",
                    "|| (x,s) - P[(x,s) - grad La(x,s)] ||_inf", "|| XS - XSp ||_inf",  "|| d ||", "step size",
                    "TR ratio", "TR radius", "status BQP");
        constexpr int maxLineWidth = 200;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    // For logging
    double pred = 0;
    double ared = 0;
    double alpha = 1.0;

    // Main loop
    size_t successive_unsuccessful = 0;
    double ngproj = ngproj0;
    double distXInnerLoop = NOMAD::INF;
    auto status = BoundAugLagSolverStatus::MAX_ITER_REACHED;
    auto bcqp_status = BCQPSolverStatus::UNDEFINED;
    for (int iter = 0; iter < max_iter_inner; ++iter)
    {
        if (verbose)
        {
            const std::string statusLog = [](const BCQPSolverStatus status) -> std::string
            {
                if (status == BCQPSolverStatus::UNDEFINED)
                {
                    return "-";
                }
                if (status == BCQPSolverStatus::MAX_ITER_REACHED)
                {
                    return "Max iteration reached";
                }
                if (status == BCQPSolverStatus::SOLVED)
                {
                    return "Solved";
                }
                if (status == BCQPSolverStatus::STAGNATION_ITERATES)
                {
                    return "Stagnation iterates";
                }
                if (status == BCQPSolverStatus::TIGHT_VAR_BOUNDS)
                {
                    return "Bounds too tight";
                }
                return "error";
            }(bcqp_status);

            if (pred != 0)
            {
                std::printf(" %-12d %12zu %14zu %24e %28e %28e %16e %10e %15e %10e %8s\n",
                            iter, successive_unsuccessful, iter_non_monotone_steps,
                            f_current, ngproj, distXInnerLoop, d.norm_inf(),
                            alpha, ared / pred, delta, statusLog.c_str());
            }
            else
            {
                std::printf(" %-12d %12zu %14zu %24e %28e %28e %16e %10e %13e/0 %10e %8s\n",
                            iter, successive_unsuccessful, iter_non_monotone_steps,
                            f_current, ngproj, distXInnerLoop, d.norm_inf(),
                            alpha, ared, delta, statusLog.c_str());
            }
        }

        // Detect stopping criteria
        // Optimality conditions
        if ((ngproj <= tol_opt * ngproj0) || (ngproj <= tol_success))
        {
            status = BoundAugLagSolverStatus::SOLVED;
            break;
        }

        // Stagnation
        if ((distXInnerLoop <= tol_dist_successive_x) || (successive_unsuccessful > max_unsuccessful_iter))
        {
            if (status != BoundAugLagSolverStatus::ONE_STEP_MADE)
            {
                status = BoundAugLagSolverStatus::STAGNATION_ITERATES;
            }
            break;
        }

        // Compute lower and upper bounds for the BCQP
        for (int i = 0; i < nvar; ++i)
        {
            dlvar.set(i, 0, std::max(lvar.get(i, 0) - XS.get(i, 0), -delta));
            duvar.set(i, 0, std::min(uvar.get(i, 0) - XS.get(i, 0), delta));
        }

        bool checkdFeasible = true;
        for (int i = 0; i < nvar; ++i)
        {
            const double dli = dlvar.get(i, 0);
            const double dui = duvar.get(i, 0);
            const double di = d.get(i, 0);
            if ((di < dli) || (di > dui))
            {
                checkdFeasible = false;
                break;
            }
        }

        // When XS has not been updated (meaning the gradient and the hessian of the trust-region model is still
        // the same) and the last direction d is still feasible, we do not need to solve again the BCQP sub problem,
        // as the solution will be the same.
        bcqp_status = BCQPSolverStatus::UNDEFINED;
        if ((successive_unsuccessful == 0) || !checkdFeasible)
        {
            d.fill(0);

            BCQPSolver bcqp_solver = {tol_dist_successive_x, max_iter_bcqp, verbose_level > 1 ? verbose_level - 2 : 0};
            bcqp_status = bcqp_solver.solve(d, BCQP_model, dlvar, duvar);
            if (bcqp_status == NOMAD::BCQPSolverStatus::BOUNDS_ERROR ||
                bcqp_status == NOMAD::BCQPSolverStatus::NUM_ERROR ||
                bcqp_status == NOMAD::BCQPSolverStatus::MATRIX_DIMENSIONS_FAILURE)
            {
                status = BoundAugLagSolverStatus::NUM_ERROR;
                break;
            }
            if (bcqp_status == NOMAD::BCQPSolverStatus::TIGHT_VAR_BOUNDS)
            {
                if (status != BoundAugLagSolverStatus::ONE_STEP_MADE)
                {
                    status = BoundAugLagSolverStatus::STAGNATION_ITERATES;
                }
                break;
            }
        };

        // Trust-region update
        // Step 1: Compute trust-region ratio = ared / pred
        // The formula takes into account numerical errors and is inspired by the technique described
        // in Section 17.4.2 of:
        //
        // Trust-region methods
        // by A.R. Conn, N.I.M. Gould and P. Toint (2000)
        //
        // and the package SolverTools.jl
        XScan = XS;
        XScan.add(d);
        projectOnBounds(XScan, lvar, uvar);
        const double f_trial = computeAugLagObj(QPModel, XScan, lambda, mu);

        // Compute magic step: reset s_can to s_can := max(0, -c(x_can) + lambda mu)
        for (int i = 0; i < n; ++i)
        {
            xcan.set(i, 0, XScan.get(i,  0));
        }
        for (int j = 0; j < ncons; ++j)
        {
            const double cj = QPModelUtils::getModelCons(QPModel, j, xcan);
            const double lambdaj = lambda.get(j, 0);
            XScan.set(j + n, 0, std::max(0.0, -cj + lambdaj * mu));
        }
        const double f_trial_magic = computeAugLagObj(QPModel, XScan, lambda, mu);

        const double qm = QPModelUtils::getModelObj(BCQP_model, d);
        pred = -qm + f_trial - f_trial_magic + tol_tr_ratio * std::max(1.0, std::fabs(f_current));

        // Compute ared := f(x) - f(x+d)
        ared = f_current - f_trial_magic + tol_tr_ratio * std::max(1.0, std::fabs(f_current));

        // Correct for rounding errors
        if ((std::fabs(qm) < 1000 * tol_tr_ratio) || (std::fabs(ared) < 1000 * tol_tr_ratio * std::fabs(f_current)))
        {
            const double slope = SGTELIB::Matrix::dot(gradLa, d);
            computeAugLagGrad(gradLaXScan, QPModel, XScan, lambda, mu);
            const double slope_trial = SGTELIB::Matrix::dot(gradLaXScan, d);
            ared = (slope_trial + slope) / 2.0;
        }

        // Update non-monotone strategy parameters
        rho_his = (f_ref - f_trial_magic) / (sig_ref + pred);
        const double rho = std::max(rho_his, ared / pred);

        // Update trust-region radius and candidate
        alpha = 1.0;
        if (rho >= epsilon_1)
        {
            // Set d according to magic step
            for (int j = 0; j < nvar; ++j)
            {
                d.set(j, 0, XScan.get(j, 0) - XS.get(j, 0));
            }

            // Accept the point
            XS = XScan;

            const double nd = d.norm_inf();
            if (rho >= epsilon_2)
            {
                delta = std::min(gamma_2 * std::max(delta, nd), max_tr_radius);
            }

            // Update non monotonicity parameters (taken from npl.py)
            sig_ref += pred;
            sig_can += pred;
            if (f_trial_magic < f_min)
            {
                f_can = f_trial_magic;
                f_min = f_trial_magic;
                sig_can = 0;
                iter_non_monotone_steps = 0;
            }
            else
            {
                iter_non_monotone_steps += 1;
            }

            if (f_trial_magic > f_can)
            {
                f_can = f_trial_magic;
                sig_can = 0;
            }

            if (iter_non_monotone_steps == max_iter_non_monotone_steps)
            {
                f_ref = f_can;
                sig_ref = sig_can;
            }

            successive_unsuccessful = 0;

            // Update model
            computeBoundAugLagQPModel(BCQP_model, QPModel, XS, lambda, mu);
            computeAugLagGrad(gradLa, QPModel, XS, lambda, mu);
            f_current = computeAugLagObj(QPModel, XS, lambda, mu);

            // Update stopping criteria metrics
            ngproj = computeFirstOrderError(XS, gradLa, lvar, uvar);
            distXInnerLoop = sqrt(alpha) * d.norm();

            status = BoundAugLagSolverStatus::ONE_STEP_MADE;

            continue;
        }

        // Failure
        // First try a backtracking linesearch along the direction d: follow Nocedal and Yuan
        const double slope = SGTELIB::Matrix::dot(d, gradLa);
        constexpr size_t max_iter_backtracking = 5;
        size_t iter_backtracking = 0;
        const double f_start = f_current;
        double f_trial_backtracking = f_start;
        bool satisfies_Armijo_condition = f_trial_backtracking <= f_start + 1e-4 * slope;
        while ((iter_backtracking < max_iter_backtracking) && !satisfies_Armijo_condition)
        {
            // Set XScan := XS + alpha d
            alpha /= 1.2;
            XScan = d;
            XScan.multiply(alpha);
            XScan.add(XS);

            // For the slack variables, apply the magic step to decrease the most efficiently
            // possible the function
            for (int i = 0; i < n; ++i)
            {
                xcan.set(i, 0, XScan.get(i,  0));
            }
            for (int j = 0; j < ncons; ++j)
            {
                const double cj = QPModelUtils::getModelCons(QPModel, j, xcan);
                const double lambdaj = lambda.get(j, 0);
                XScan.set(j + n, 0, std::max(0.0, -cj + lambdaj * mu));
            }
            projectOnBounds(XScan, lvar, uvar);
            f_trial_backtracking = computeAugLagObj(QPModel, XScan, lambda, mu);
            satisfies_Armijo_condition = f_trial_backtracking <= f_start + 1e-4 * slope;
            iter_backtracking += 1;
        }

        if (satisfies_Armijo_condition)
        {
            // Set d according to magic step
            for (int j = 0; j < nvar; ++j)
            {
                d.set(j, 0, XScan.get(j, 0) - XS.get(j, 0));
            }
            const double nd = d.norm_inf();

            // Accept new candidate
            XS = XScan;
            successive_unsuccessful = 0;

            // Update delta
            delta = std::max(std::min(alpha * nd, delta), min_tr_radius);

            // Update model
            computeBoundAugLagQPModel(BCQP_model, QPModel, XS, lambda, mu);
            computeAugLagGrad(gradLa, QPModel, XS, lambda, mu);
            f_current = computeAugLagObj(QPModel, XS, lambda, mu);

            // Update stopping criteria metrics
            ngproj = computeFirstOrderError(XS, gradLa, lvar, uvar);
            distXInnerLoop = d.norm();

            status = BoundAugLagSolverStatus::ONE_STEP_MADE;
            continue;
        }

        // Do not accept the candidate and reduce trust-region radius
        const double nd = d.norm_inf();
        delta = std::max(gamma_1 * std::min(delta, nd), min_tr_radius);
        distXInnerLoop = sqrt(alpha) * d.norm();

        // NB: No need to update the model and gradient as XS remains identical
        successive_unsuccessful += 1;
    }

    if (verbose)
    {
        std::printf("\nStatus: ");
        std::printf("La(x*, s*) = %e\n", f_current);
        if (status == BoundAugLagSolverStatus::SOLVED)
        {
            std::printf("Has reached the minimum tolerance:\n");
            std::printf("|| (x s) - P[(x,s) - grad La(x,s)] ||_inf = %e <= tol_success = %e\n",
                        ngproj,  std::max(tol_opt * ngproj0, tol_success));
        }
        else if (status == BoundAugLagSolverStatus::STAGNATION_ITERATES)
        {
            std::printf("Inner steps have stagnated:\n");
            std::printf("|| (x,s) - (xp, sp) || = %e <= %e or\n", distXInnerLoop, tol_dist_successive_x);
            std::printf("the maximum number of unsuccessful iterations %zu has been reached\n",
                        max_unsuccessful_iter);
        }
        else if (status == BoundAugLagSolverStatus::ONE_STEP_MADE)
        {
            std::printf("At least one step has been made\n");
        }
        else if (status == BoundAugLagSolverStatus::FAILURE)
        {
            std::printf("Has failed\n");
        }
        else
        {
            std::printf("Unknown stopping criterion\n");
        }
        std::printf("\n");
    }

    return status;
}

bool NOMAD::AugLagSolver::checkParameters() const
{
    return true;
}

bool NOMAD::AugLagSolver::checkDimensions(const SGTELIB::Matrix& x,
                                          const SGTELIB::Matrix& QPModel,
                                          const SGTELIB::Matrix& lb,
                                          const SGTELIB::Matrix& ub)
{
    const int n = x.get_nb_rows();
    if (n != std::max(x.get_nb_rows(), x.get_nb_cols()) && (x.get_nb_cols() != 1))
    {
        std::string err = "AugLagSolver::solve error: x must be a column vector";
        std::printf("%s\n", err.c_str());
        return false;
    }

    if (n != lb.get_nb_rows() || n != ub.get_nb_rows())
    {
        std::string err = "AugLagSolver::solve error: bound constraints dimensions ";
        err += "nlb = " + std::to_string(lb.get_nb_cols()) + " nub = " + std::to_string(ub.get_nb_cols());
        err += " are not compatible with dimension of x (n = " + std::to_string(n) + ")";
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbParams = QPModel.get_nb_cols();
    if (nbParams != (n + 1) + n * (n + 1) / 2)
    {
        std::string err = "AugLagSolver::solve error: ";
        err += "the number of params of the model nbParams = (n+1) * (n+2) / 2 = " + std::to_string(nbParams);
        err += " is not compatible with the dimension of the solution n = " + std::to_string(n);
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbCons = QPModel.get_nb_rows() - 1;
    if (nbCons < 1)
    {
        std::string err = "AugLagSolver::solve error: ";
        err += "the model has no constraints";
        std::printf("%s\n", err.c_str());
        return false;
    }

    return true;
}

bool NOMAD::AugLagSolver::checkBoundsCompatibilities(const SGTELIB::Matrix& lb,
                                                     const SGTELIB::Matrix& ub)
{
    const int n = lb.get_nb_cols();
    for (int i = 0; i < n; ++i)
    {
        const bool areBoundsCompatible = lb.get(i, 0) <= ub.get(i, 0);
        if (!areBoundsCompatible)
        {
            std::string err = "AugLagSolver::solve error: ";
            err += "no compatibility between lower bound and upper bound for index " + std::to_string(i);
            std::printf("%s\n", err.c_str());
            return false;
        }
    }
    return true;
}

void NOMAD::AugLagSolver::computeBoundAugLagQPModel(SGTELIB::Matrix& model,
                                                    const SGTELIB::Matrix& QPModel,
                                                    const SGTELIB::Matrix& XS,
                                                    const SGTELIB::Matrix& lambda,
                                                    const double mu)
{
    const int nvar = XS.get_nb_rows();

    // Construct Q(d) = grad La(x, s)' d + 1/2 d' H Lag(x, s) d
    SGTELIB::Matrix gradLa("gradLa", nvar, 1);
    computeAugLagGrad(gradLa, QPModel, XS, lambda, mu);
    SGTELIB::Matrix HLa("HLa", nvar, nvar);
    computeAugLagHessian(HLa, QPModel, XS, lambda, mu);

    // Scalar part
    model.set(0, 0, 0);

    // Linear part of Q(x)
    for (int j = 0; j < nvar; ++j)
    {
        model.set(0, j + 1, gradLa.get(j, 0));
    }

    // Quadratic part of Q(x)
    int k = nvar + 1;
    for (int j = 0; j < nvar; ++j)
    {
       model.set(0, k, HLa.get(j, j));
       ++k;
    }

    for (int i = 0; i < nvar; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            model.set(0, k, HLa.get(i, j));
            ++k;
        }
    }
}

double NOMAD::AugLagSolver::computeAugLagObj(const SGTELIB::Matrix& QPModel,
                                             const SGTELIB::Matrix& XS,
                                             const SGTELIB::Matrix& lambda,
                                             const double mu)
{
    const int nvar = XS.get_nb_rows();
    const int ncons = QPModel.get_nb_rows() - 1;
    const int n = nvar - ncons;

    SGTELIB::Matrix x("x", n, 1);
    for (int i = 0; i < n; ++i)
    {
        x.set(i, 0, XS.get(i, 0));
    }
    SGTELIB::Matrix s("s", ncons, 1);
    for (int j = 0; j < ncons; ++j)
    {
        s.set(j, 0, XS.get(j + n, 0));
    }

    double auglagVal = QPModelUtils::getModelObj(QPModel, x);
    for (int j = 0; j < ncons; ++j)
    {
        const double sj = s.get(j, 0);
        const double lambdaj = lambda.get(j, 0);
        const double cj = QPModelUtils::getModelCons(QPModel, j, x);
        auglagVal -= lambdaj * (cj + sj);
    }

    double nsquarecslack = 0.0;
    for (int j = 0; j < ncons; ++j)
    {
        const double cj = QPModelUtils::getModelCons(QPModel, j, x);
        const double sj = s.get(j, 0);
        nsquarecslack += (cj + sj) * (cj + sj);
    }
    auglagVal += nsquarecslack / (2.0 * mu);

    return auglagVal;
}

void NOMAD::AugLagSolver::computeAugLagGrad(SGTELIB::Matrix& gradL,
                                            const SGTELIB::Matrix& QPModel,
                                            const SGTELIB::Matrix& XS,
                                            const SGTELIB::Matrix& lambda,
                                            const double mu)
{
    const int nvar = XS.get_nb_rows();
    const int ncons = QPModel.get_nb_rows() - 1;
    const int n = nvar - ncons;

    SGTELIB::Matrix x("x", n, 1);
    for (int i = 0; i < n; ++i)
    {
        x.set(i, 0, XS.get(i, 0));
    }
    SGTELIB::Matrix s("s", ncons, 1);
    for (int j = 0; j < ncons; ++j)
    {
        s.set(j, 0, XS.get(j + n, 0));
    }

    // Compute grad La(x) = grad L(x) + (1 / mu) J(x)' (c(x) + s)
    SGTELIB::Matrix gradLx("gradLx", n, 1);
    QPModelUtils::getModelLagrangianGrad(gradLx, QPModel, x, lambda);

    SGTELIB::Matrix cslack("cons", ncons, 1);
    QPModelUtils::getModelCons(cslack, QPModel, x);
    for (int j = 0; j < ncons; ++j)
    {
        const double cj = cslack.get(j, 0);
        cslack.set(j, 0, cj + s.get(j, 0));
    }

    SGTELIB::Matrix Jx("Jx", ncons, n);
    QPModelUtils::getModelJacobianCons(Jx, QPModel, x);

    SGTELIB::Matrix gradLax("gradLax", n, 1);
    SGTELIB::Matrix::inplace_product(gradLax, Jx.transpose(), cslack);
    gradLax.multiply(1.0 / mu);
    gradLax.add(gradLx);

    // Compute grad La(s) = -lambda + (1/mu) (c(x) + s)
    SGTELIB::Matrix gradLas("gradLas", ncons, 1);
    for (int j = 0; j < ncons; ++j)
    {
        const double cslackj = cslack.get(j, 0);
        const double gj = -lambda.get(j, 0) + cslackj / mu;
        gradLas.set(j, 0, gj);
    }

    // Compute grad La(x, s) = [grad La(x) grad La(s)]^T
    for (int i = 0; i < n; ++i)
    {
        gradL.set(i, 0, gradLax.get(i, 0));
    }
    for (int j = 0; j < ncons; ++j)
    {
        gradL.set(j + n, 0, gradLas.get(j, 0));
    }
}

void NOMAD::AugLagSolver::computeAugLagHessian(SGTELIB::Matrix& HLag,
                                               const SGTELIB::Matrix& QPModel,
                                               const SGTELIB::Matrix& XS,
                                               const SGTELIB::Matrix& lambda,
                                               const double mu)
{
    const int nvar = XS.get_nb_rows();
    const int ncons = QPModel.get_nb_rows() - 1;
    const int n = nvar - ncons;

    SGTELIB::Matrix x("x", n, 1);
    for (int i = 0; i < n; ++i)
    {
        x.set(i, 0, XS.get(i, 0));
    }
    SGTELIB::Matrix s("s", ncons, 1);
    for (int j = 0; j < ncons; ++j)
    {
        s.set(j, 0, XS.get(j + n, 0));
    }

    // Compute H La(x, s) = H f(x, s) - sum y_j Hcslackj(x, s) + 1/mu J(x, s)^T J(x, s)
    // where:
    // H f(x, s) = [ H(x)  0 ], Hcslack(x, s) = [ Hcj(x) 0 ], J(x, s) = [ J(x) I ]
    //             [  0    0 ]                  [  0     0 ]
    // and y_j = lambdaj - 1/mu (cj(x) + sj)
    // NB: This formula is taken from:
    // CME338, Large-Scale Numerical Optimization, Note 11 by M. Saunders.
    HLag.fill(0.0);
    auto addHx = [](SGTELIB::Matrix& HLag, const SGTELIB::Matrix& H)
    {
        const int n = H.get_nb_rows();
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                HLag.set(i, j, HLag.get(i, j) + H.get(i, j));
            }
        }
    };

    SGTELIB::Matrix H("H", n, n);
    QPModelUtils::getModelObjHessian(H, QPModel, x);
    addHx(HLag, H);
    for (int j = 0; j < ncons; ++j)
    {
        const double sj = s.get(j, 0);
        const double cj = QPModelUtils::getModelCons(QPModel, j, x);
        const double lambdaj = lambda.get(j, 0);
        const double yj = lambdaj - (cj + sj) / mu;
        QPModelUtils::getModelConsHessian(H, QPModel, j, x);
        H.multiply(-yj);
        addHx(HLag, H);
    }
    SGTELIB::Matrix Jxs("Jxs", ncons, nvar);
    Jxs.fill(0);
    SGTELIB::Matrix Jx("Jx", ncons, n);
    QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
    for (int i = 0; i < ncons; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            Jxs.set(i, j, Jx.get(i, j));
        }
        Jxs.set(i, i + n, 1.0);
    }
    SGTELIB::Matrix JtJ = SGTELIB::Matrix::product(Jxs.transpose(), Jxs);
    JtJ.multiply(1.0 / mu);
    HLag.add(JtJ);
}

double NOMAD::AugLagSolver::computeFirstOrderError(const SGTELIB::Matrix& x,
                                                   const SGTELIB::Matrix& gradL,
                                                   const SGTELIB::Matrix& lb,
                                                   const SGTELIB::Matrix& ub)
{
    const int n = x.get_nb_rows();
    SGTELIB::Matrix dualFeas("dualFeas", n, 1);
    for (int i = 0; i < n; ++i)
    {
        dualFeas.set(i, 0, x.get(i, 0) - gradL.get(i, 0));
    }
    projectOnBounds(dualFeas, lb, ub);
    dualFeas.sub(x);

    return dualFeas.norm_inf();
}

void NOMAD::AugLagSolver::projectOnBounds(SGTELIB::Matrix& x,
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
