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
 \file   TRIPMSolver.cpp
 \brief  Trust-Region Interior Point Method: implementation
 \author Tangi Migot and Ludovic Salomon
 \see    TRIPMSolver.hpp
 */
#include "TRIPMSolver.hpp"

#include "../../Nomad/nomad.hpp"
#include "../../Algos/QPSolverAlgo/DoglegTRSolver.hpp"
#include "../../Algos/QPSolverAlgo/LevenbergMarquardtSolver.hpp"
#include "../../Algos/QPSolverAlgo/ProjectedConjugateGradientSolver.hpp"
#include "../../Algos/QPSolverAlgo/QPModelUtils.hpp"
#include "../../Math/MatrixUtils.hpp"

NOMAD::TRIPMSolverStatus NOMAD::TRIPMSolver::solve(SGTELIB::Matrix& x,
                                                   const SGTELIB::Matrix& QPModel,
                                                   const SGTELIB::Matrix& lb,
                                                   const SGTELIB::Matrix& ub) const
{
    if (!checkParams())
    {
        return NOMAD::TRIPMSolverStatus::PARAM_ERROR;
    }

    if (!checkDimensions(x, QPModel, lb, ub))
    {
        return NOMAD::TRIPMSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }

    if (!checkBoundsCompatibilities(lb, ub))
    {
        return NOMAD::TRIPMSolverStatus::BOUNDS_ERROR;
    }

    // Fix tolerance
    constexpr double atol = 1e-8;

    // Fix verbosity
    const bool verbose = verbose_level > 0;

    const int n = x.get_nb_rows();

    // Collect all fixed variables
    size_t nbFixedVars = 0;
    std::vector<bool> fixedVars(n, false);
    for (int i = 0; i < n; ++i)
    {
        if (std::abs(ub.get(i, 0) - lb.get(i, 0)) <= atol)
        {
            nbFixedVars += 1;
            fixedVars[i] = true;
        }
    }

    // When the difference between lower and upper bounds is too small, stop the algorithm
    if (nbFixedVars == n)
    {
        return NOMAD::TRIPMSolverStatus::TIGHT_VAR_BOUNDS;
    }

    if (verbose)
    {
        std::printf("\nTrust-region interior point method algorithm\n");
        std::printf("Number of total variables: %d\n", n);
        std::printf("Number of fixed variables: %zu\n", nbFixedVars);
    }

    // Solve problem in full dimensions
    if (nbFixedVars == 0)
    {
        const auto status = solveReducedPb(x, QPModel, lb, ub);
        return status;
    }

    // Solve reduced problem with fixed variables
    // NB: due to numerical errors, the values returned by the QPModel and the reduced QPModel can differ.
    const int nfree = n - (int) nbFixedVars;
    const auto QPModelReduced = QPModelUtils::getReducedQPModel(QPModel, x, fixedVars);
    SGTELIB::Matrix xRed("xRed", nfree, 1);
    SGTELIB::Matrix lbRed("lbRed", nfree, 1);
    SGTELIB::Matrix ubRed("ubRed", nfree, 1);
    int curInd = 0;
    for (int i = 0; i < n; ++i)
    {
        if (fixedVars[i])
            continue;

        xRed.set(curInd, 0, x.get(i, 0));
        lbRed.set(curInd, 0, lb.get(i, 0));
        ubRed.set(curInd, 0, ub.get(i, 0));
        curInd++;
    }

    const auto status = solveReducedPb(xRed, QPModelReduced, lbRed, ubRed);

    // NB: whatever the status (even if there is a numerical error), we recompute the full-space solution
    curInd = 0;
    for (int i = 0; i < n; ++i)
    {
        if (fixedVars[i])
            continue;

        x.set(i, 0, xRed.get(curInd, 0));
        curInd++;
    }

    return status;
}

NOMAD::TRIPMSolverStatus NOMAD::TRIPMSolver::solveReducedPb(SGTELIB::Matrix& x,
                                                            const SGTELIB::Matrix& QPModel,
                                                            const SGTELIB::Matrix& lb,
                                                            const SGTELIB::Matrix& ub) const
{
    const int n = x.get_nb_rows();
    const int nbCons = QPModel.get_nb_rows() - 1;
    const int nbVar = n + nbCons;

    // Fix verbosity
    const bool verbose = verbose_level > 0;

    // Fix tolerance
    constexpr double atol_opt = 1e-6;
    constexpr double atol_feas = 1e-6;
    constexpr double smallest_mu_tol = std::min(atol_opt, atol_feas) / 100;

    // Initialize lvar and uvar, respectively lower and upper bounds for (x, s),
    // where s are the slack variables
    SGTELIB::Matrix lvar("lvar", nbVar, 1);
    SGTELIB::Matrix uvar("uvar", nbVar, 1);
    for (int i = 0; i < n; ++i)
    {
        lvar.set(i, 0, lb.get(i, 0));
        uvar.set(i, 0, ub.get(i, 0));
    }
    for (int i = 0; i < nbCons; ++i)
    {
        lvar.set(i + n, 0, 0.0);
        uvar.set(i + n, 0, NOMAD::INF);
    }

    // Compute starting point
    SGTELIB::Matrix cons("cons", nbCons, 1);
    SGTELIB::Matrix XS("XS", nbVar, 1);
    if (!computeStrictlyFeasiblePoint(XS, cons, QPModel, lvar, uvar))
    {
        return NOMAD::TRIPMSolverStatus::STRICT_PT_FAILURE;
    }

    double mu = mu_init;
    double tol_mu = mu_init;
    NOMAD::LevenbergMarquardtSolver LMalgo = {mu, // feasibility_tol
                                              tol_mu, // tolerance
                                              1e-15, // tol_dist_successive_x
                                              30, // Max number of iterations
                                              true, // is_sol_feas
                                              verbose_level > 0 ? verbose_level - 1 : 0}; // verbose level
    LMalgo.solve(x, XS, QPModel, lvar, uvar, cons);

    // Check feasibility of the starting point: if not, stop the procedure.
    for (int i = 0; i < nbVar; ++i)
    {
        const bool feasible = (XS.get(i, 0) >= lvar.get(i, 0)) && (XS.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            return NOMAD::TRIPMSolverStatus::LM_FAILURE;
        }
    }

    // Initialize error metrics
    TRIPMErrorMetric err_metric = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < nbCons; ++i)
    {
        err_metric.cxInitNorm = std::max(err_metric.cxInitNorm, cons.get(i, 0));
    }
    TRIPMErrorMetric err_metric_barrier = {0.0, 0.0, 0.0, 0.0, 0.0};
    err_metric_barrier.cxInitNorm = err_metric.cxInitNorm;

    // Compute cslack variables, i.e., cslack = c(x) + s.
    SGTELIB::Matrix cslack("c+s", nbCons, 1);
    for (int j = 0; j < nbCons; ++j)
    {
        cslack.set(j, 0, XS.get(j + n, 0) + cons.get(j, 0));
    }
    double cx = cslack.norm();
    double cxp = cx;

    // TRIPM parameters initialization
    const size_t maxSuccessiveFailures = 3;
    size_t successiveFailure = 0;
    size_t successiveAcceptable = 0;
    size_t successiveBeforeUpdate = 2;
    double distXOuterLoop = NOMAD::INF;

    // Compute Lagrange multipliers estimates
    SGTELIB::Matrix Gk("Gk", n, 1);
    QPModelUtils::getModelObjGrad(Gk, QPModel, x);
    SGTELIB::Matrix Jx("Jx", nbCons, n);
    QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
    SGTELIB::Matrix lambda("lambda", nbCons, 1);
    computeSlackMultipliers(lambda, XS, Jx, Gk, 0.0);

    // Allocate memory for the following vectors
    SGTELIB::Matrix XSp("XSp", nbVar, 1); // Next iterate
    SGTELIB::Matrix p("p", nbVar, 1); // The primal-dual step
    SGTELIB::Matrix xp(x); xp.set_name("xp");

    // Compute outer stopping criteria
    computeErrorFunctionMetric(err_metric, XS, QPModel, lvar, uvar, lambda, 0.0, false);

    if (verbose)
    {
        std::printf("Number of inequality constraints: %d\n", nbCons);
        std::printf("Stopping criterion tolerance for optimality: %f\n", atol_opt);
        std::printf("Stopping criterion tolerance for feasibility: %f\n", std::max(1.0, err_metric.cxInitNorm) * atol_feas);
        std::printf("Maximum number of iterations allowed for outer loop: %zu\n", max_iter_outer);
        std::printf("Maximum number of iterations allowed for inner loop: %zu\n\n", max_iter_inner);

        std::printf("%10s %10s %26s %16s %22s %18s %9s %15s %20s %18s %15s\n",
                    "iter (outer)", "f(x)", "|| max(c(x), 0) ||_inf", "|| c(x) + s ||",
                    "|| Proj grad L(x) ||_inf", "|| S lambda ||_inf", "mu", "tol_mu", "|| lambda ||_inf", "|| x - xp ||_inf", "status inner");
        constexpr int maxLineWidth = 195;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    TRIPMSolverStatus status = TRIPMSolverStatus::MAX_ITER_REACHED;
    BarrierSolverStatus innerStatus = BarrierSolverStatus::UNDEFINED;
    for (int iter = 0; iter < (int)max_iter_outer; ++iter)
    {
        if (verbose)
        {
            const double objVal = QPModelUtils::getModelObj(QPModel, x);
            double normCxMax0 = 0.0;
            for (int i = 0; i < nbCons; ++i)
            {
                normCxMax0 = std::max(cons.get(i, 0), normCxMax0);
            }

            const std::string statusLog = [](const BarrierSolverStatus status) -> std::string
            {
                if (status == BarrierSolverStatus::UNDEFINED)
                {
                    return "-";
                }
                if (status == BarrierSolverStatus::SOLVED)
                {
                    return "Solved";
                }
                if (status == BarrierSolverStatus::ONE_STEP_MADE)
                {
                    return "Has progressed";
                }
                if (status == BarrierSolverStatus::STAGNATION_ITERATES)
                {
                    return "Stagnation iterates";
                }
                return "error";
            }(innerStatus);

            std::printf(" %-9d %+15e %15e %23e %18e %21e %18e %13e %14e %18e %18s\n",
                        iter, objVal, normCxMax0, cx, err_metric.projlagGradNorm, err_metric.slackLambdaMuNorm, mu, tol_mu,
                        lambda.norm_inf(), distXOuterLoop, statusLog.c_str());
        }

        // This stopping criterion is adapted from:
        //
        //  "An interior algorithm for nonlinear optimization that combines line search and trust region steps"
        //  by R. Waltz, J. Morales, J. Nocedal, and D. Orban
        //  Math. Program. 107, 391–408 (2006). https://doi.org/10.1007/s10107-004-0560-5
        //
        const bool outerSuccess = err_metric.projlagGradNorm <= std::max(1.0, err_metric.projObjGrad) * atol_opt &&
                                  err_metric.slackLambdaMuNorm <= std::max(1.0, err_metric.projObjGrad) * atol_opt &&
                                  err_metric.cxNorm <= std::max(1.0, err_metric.cxInitNorm) * atol_feas;
        if (outerSuccess)
        {
            status = TRIPMSolverStatus::SOLVED;
            break;
        }

        const bool outerFailure = (distXOuterLoop <= tol_dist_successive_x) ||
                                  (mu <= std::min(atol_opt, atol_feas) / mu_decrease) ||
                                  (successiveFailure >= maxSuccessiveFailures);
        if (outerFailure)
        {
            status = TRIPMSolverStatus::STAGNATION_ITERATES;
            break;
        }

        // Run inner barrier solver
        innerStatus = solveBarrierSubproblem(xp, XSp, cslack, lambda, Gk,
                                             cons, Jx, err_metric_barrier, QPModel, XS, lvar, uvar,
                                             mu, atol_opt, atol_feas,
                                             (int) verbose_level - 1);
        if (innerStatus == BarrierSolverStatus::NUM_ERROR)
        {
            verbose && std::printf("TRIPMSolver::solve numerical error: stop the procedure\n");
            return TRIPMSolverStatus::NUM_ERROR;
        }

        const bool innerSuccess = innerStatus == BarrierSolverStatus::SOLVED ||
                                  innerStatus == BarrierSolverStatus::ONE_STEP_MADE ||
                                  innerStatus == BarrierSolverStatus::STAGNATION_ITERATES;

        for (int i = 0; i < n; ++i)
        {
            xp[i] = XSp.get(i, 0);
        }
        QPModelUtils::getModelCons(cons, QPModel, xp);
        for (int j = 0; j < nbCons; ++j)
        {
            cslack.set(j, 0, XSp.get(j + n, 0) + cons.get(j, 0));
        }
        cxp = cslack.norm();

        // Update
        if (innerSuccess)
        {
            distXOuterLoop = SGTELIB::Matrix::distNorm2(XS, XSp);
            XS = XSp;
            x = xp;
            cx = cxp;
            successiveFailure = 0;
            if (innerStatus == BarrierSolverStatus::SOLVED)
            {
                successiveAcceptable = 0;
                mu /= mu_decrease;
                tol_mu /= mu_decrease;
            }
            else if (successiveAcceptable >= successiveBeforeUpdate)
            { // slower decrease
                successiveAcceptable = 0;
                mu /= sqrt(mu_decrease);
                tol_mu /= sqrt(mu_decrease);
            }
            else
            {
                // Try to re-run from new point.
                successiveAcceptable += 1;
            }
            // Compute stopping criterion
            QPModelUtils::getModelObjGrad(Gk, QPModel, x);
            // ng = Gk.norm_inf();
            // tol = std::max(ng, 1.0) * 1e-6;
            computeErrorFunctionMetric(err_metric, XS, QPModel, lvar, uvar, lambda, 0.0, false);

            // Recompute Jacobian of the constraints
            QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
        }
        else
        {
            successiveAcceptable = 0;
            successiveFailure += 1;

            if (cxp > tol_mu)
            {
                LMalgo.feasibility_tol = mu;
                LMalgo.tol = tol_mu;
                const auto LMStatus = LMalgo.solve(xp, XSp, QPModel, lvar, uvar, cons);

                if (LMStatus == LMSolverStatus::SOLVED)
                {
                    XS = XSp;
                    x = xp;
                    for (int j = 0; j < nbCons; ++j)
                    {
                        cslack.set(j, 0, XS.get(j + n, 0) + cons.get(j, 0));
                    }
                    cxp = cslack.norm();
                    cx = cxp;

                    // Recompute Lagrange multipliers
                    QPModelUtils::getModelObjGrad(Gk, QPModel, x);
                    QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
                    computeSlackMultipliers(lambda, XS, Jx, Gk, 0.0);

                    // Compute stopping criterion
                    computeErrorFunctionMetric(err_metric, XS, QPModel, lvar, uvar, lambda, 0.0, false);

                    successiveFailure = 0;
                }
                else
                {
                    mu /= mu_decrease;
                    tol_mu /= mu_decrease;
                }
            }
            else
            {
                mu /= mu_decrease;
                tol_mu /= mu_decrease;
            }
        }
        tol_mu = std::max(smallest_mu_tol, tol_mu);
    }

    if (verbose)
    {
        std::printf("\nStatus: ");
        std::printf("f(x*) = %e\n", QPModelUtils::getModelObj(QPModel, x));
        std::printf("|| c(x*) + s* || = %e\n", cxp);
        double normCxMax0 = 0;
        for (int i = 0; i < nbCons; ++i)
        {
            normCxMax0 = std::max(cons.get(i, 0), normCxMax0);
        }
        std::printf("|| max (c(x), 0) ||_inf = %e\n", normCxMax0);
        if (status == TRIPMSolverStatus::SOLVED)
        {
            std::printf("Has reached the minimum tolerance\n");
            std::printf("|| x - P[x - grad L(x)] ||_inf = %e <= max(1.0, || x - P[x - grad f(x)] ||) * atol_opt = %e and\n",
                        err_metric.projlagGradNorm, std::max(1.0, err_metric.projObjGrad) * atol_opt);
            std::printf("|| S lambda ||_inf = %e <= max(1.0, || x - P[x - grad f(x)] ||_inf) * atol_opt = %e and\n",
                        err_metric.slackLambdaMuNorm, std::max(1.0, err_metric.projObjGrad) * atol_opt);
            std::printf("|| max(0, c(x)) ||_inf = %e <= max(1.0, || max(0, c(x0)) ||_inf) * atol_feas = %e\n",
                        err_metric.cxNorm, std::max(1.0, err_metric.cxInitNorm) * atol_feas);
        }
        else if (status == TRIPMSolverStatus::STAGNATION_ITERATES)
        {
            std::printf("Outer steps have stagnated:\n");
            std::printf("|| x - xp || = %e <= %e or\n", distXOuterLoop, tol_dist_successive_x);
            std::printf("barrier parameter mu = %e is below the following tolerance %e or\n", mu,
                        std::min(atol_opt, atol_feas) / mu_decrease);
            std::printf("the number of successive failure %zu has reached its maximum value %zu",
                        successiveFailure, maxSuccessiveFailures);
        }
        else if (status == TRIPMSolverStatus::MAX_ITER_REACHED)
        {
            std::printf("The maximum number of outer iterations has been reached\n");
        }
        else
        {
            std::printf("Unknown stopping criterion\n");
        }
    }
    return status;
}

bool NOMAD::TRIPMSolver::checkDimensions(const SGTELIB::Matrix& x,
                                         const SGTELIB::Matrix& QPModel,
                                         const SGTELIB::Matrix& lb,
                                         const SGTELIB::Matrix& ub)
{
    const int n = x.get_nb_rows();
    if (n != std::max(x.get_nb_rows(), x.get_nb_cols()) && (x.get_nb_cols() != 1))
    {
        std::string err = "TRIPMSolver::solve error: x must be a column vector";
        std::printf("%s\n", err.c_str());
        return false;
    }

    if (n != lb.get_nb_rows() || n != ub.get_nb_rows())
    {
        std::string err = "TRIPMSolver::solve error: bound constraints dimensions ";
        err += "nlb = " + std::to_string(lb.get_nb_cols()) + " nub = " + std::to_string(ub.get_nb_cols());
        err += " are not compatible with dimension of x (n = " + std::to_string(n) + ")";
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbParams = QPModel.get_nb_cols();
    if (nbParams != (n+1) + n * (n+1) / 2)
    {
        std::string err = "TRIPMSolver::solve error: ";
        err += "the number of params of the model nbParams = (n+1) * (n+2) / 2 = " + std::to_string(nbParams);
        err += " is not compatible with the dimension of the solution n = " + std::to_string(n);
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbCons = QPModel.get_nb_rows() - 1;
    if (nbCons < 1)
    {
        std::string err = "TRIPMSolver::solve error: ";
        err += "the model has no constraints";
        std::printf("%s\n", err.c_str());
        return false;
    }

    return true;
}

bool NOMAD::TRIPMSolver::checkBoundsCompatibilities(const SGTELIB::Matrix& lb,
                                                    const SGTELIB::Matrix& ub)
{
    const int n = lb.get_nb_cols();
    for (int i = 0; i < n; ++i)
    {
        const bool areBoundsCompatible = lb.get(i, 0) <= ub.get(i, 0);
        if (!areBoundsCompatible)
        {
            std::string err = "TRIPMSolver::solve error: ";
            err += "no compatibility between lower bound and upper bound for index " + std::to_string(i);
            std::printf("%s\n", err.c_str());
            return false;
        }
    }
    return true;
}

bool NOMAD::TRIPMSolver::checkParams() const
{
    if (mu_decrease <= 1)
    {
        std::string err = "TRIPMSolver::solve error: ";
        err += "mu_decrease parameter value must be superior to 1";
        std::printf("%s\n", err.c_str());
        return false;
    }

    return true;
}

bool NOMAD::TRIPMSolver::computeStrictlyFeasiblePoint(SGTELIB::Matrix& x,
                                                      SGTELIB::Matrix& cons,
                                                      const SGTELIB::Matrix& QPModel,
                                                      const SGTELIB::Matrix& lvar,
                                                      const SGTELIB::Matrix& uvar)
{
    const int nbCons = cons.get_nb_rows();
    const int n = uvar.get_nb_rows() - nbCons;

    SGTELIB::Matrix xtmp("xtmp", n, 1);
    for (int i = 0; i < n; ++i)
    {
        double xi = x.get(i, 0);
        const double lb = lvar.get(i, 0);
        const double ub = uvar.get(i, 0);
        if ((xi <= lb) || (xi >= ub))
        {
            if ((lb > NOMAD::M_INF) && (ub == NOMAD::INF))
            {
                xi = lb + 0.5;
            }
            else if ((lb == NOMAD::M_INF) && (ub < NOMAD::INF))
            {
                xi = ub - 0.5;
            }
            else if ((lb > NOMAD::M_INF) && (ub < NOMAD::INF))
            {
                const double mid = uvar.get(i, 0) - lvar.get(i, 0);
                xi = lb + mid / 2;
            }
            else
            {
                xi = 0.0;
            }
        }
        x.set(i, 0, xi);
        xtmp.set(i, 0, xi);
    }

    QPModelUtils::getModelCons(cons, QPModel, xtmp);
    for (int j = 0; j < nbCons; ++j)
    {
        x.set(j + n, 0, std::max(-cons.get(j, 0), 0.5)); // s = -c(x) or s > 0
    }

    for (int i = 0; i < n; i++)
    {
        const double xi = x.get(i, 0);
        const double ui = uvar.get(i, 0);
        const double li = lvar.get(i, 0);
        const bool isStrictlyFeasible = (xi > li) && (xi < ui);
        if (!isStrictlyFeasible)
        {
            std::string err = "TRIPMSolver::solve warning: ";
            err += "x is not strictly feasible for index variable i = " + std::to_string(i);
            std::printf("%s: lb[i] = %f, ub[i] = %f, x[i] = %f\n", err.c_str(), li, ui, xi);
            return false;
        }
    }

    return true;
}

void NOMAD::TRIPMSolver::computeErrorFunctionMetric(NOMAD::TRIPMSolver::TRIPMErrorMetric& errMetric,
                                                    const SGTELIB::Matrix& XS,
                                                    const SGTELIB::Matrix& QPModel,
                                                    const SGTELIB::Matrix& lvar,
                                                    const SGTELIB::Matrix& uvar,
                                                    const SGTELIB::Matrix& lambda,
                                                    const double mu,
                                                    const bool isBarrierProblem)
{
    const int nbVar = XS.get_nb_rows();
    const int nbCons = QPModel.get_nb_rows() - 1;
    const int n = nbVar - nbCons;

    SGTELIB::Matrix x("x", n, 1);
    for (int i = 0; i < n; ++i)
    {
        x.set(i, 0, XS.get(i, 0));
    }

    // Compute the gradient of the lagrangian with respect to x
    SGTELIB::Matrix gradLx("gradLx", n, 1);
    QPModelUtils::getModelLagrangianGrad(gradLx, QPModel, x, lambda);

    // Compute || x - P[x - grad L(x)] ||_inf, where P[X] is the projection of x on [lvar, uvar]
    errMetric.projlagGradNorm = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double dualFeas = x.get(i, 0) - gradLx.get(i, 0);
        dualFeas = std::clamp(dualFeas, lvar.get(i, 0), uvar.get(i, 0));
        dualFeas = x.get(i, 0) - dualFeas;
        errMetric.projlagGradNorm = std::max(errMetric.projlagGradNorm, std::abs(dualFeas));
    }

    // Compute || -S y - mu ||_inf
    errMetric.slackLambdaMuNorm = 0.0;
    for (int i = 0; i < nbCons; i++)
    {
        errMetric.slackLambdaMuNorm = std::max(std::abs(-XS.get(i + n, 0) * lambda.get(i, 0) - mu),
                                               errMetric.slackLambdaMuNorm);
    }

    // Compute || x - P[x - grad f(x)] ||_inf, where P[X] is the projection of x on [lvar, uvar]
    SGTELIB::Matrix gradObj("gradObj", n, 1);
    QPModelUtils::getModelObjGrad(gradObj, QPModel, x);
    for (int i = 0; i < n; ++i)
    {
        double dualFeas = x.get(i, 0) - gradObj.get(i, 0);
        dualFeas = std::clamp(dualFeas, lvar.get(i, 0), uvar.get(i, 0));
        dualFeas = x.get(i, 0) - dualFeas;
        errMetric.projObjGrad = std::max(errMetric.projObjGrad, std::abs(dualFeas));
    }

    // Compute || c(x) + s ||_inf (for the barrier subproblem) or || max(c(x), 0) ||_inf (for the outer loop)
    SGTELIB::Matrix cons("cons", nbCons, 1);
    QPModelUtils::getModelCons(cons, QPModel, x);
    errMetric.cxNorm = 0.0;
    if (isBarrierProblem)
    {
        for (int i = 0; i < nbCons; ++i)
        {
            const double si = XS.get(i + n, 0);
            errMetric.cxNorm = std::max(std::abs(cons.get(i, 0) + si), errMetric.cxNorm);
        }
    }
    else
    {
        for (int i = 0; i < nbCons; ++i)
        {
            errMetric.cxNorm = std::max(cons.get(i, 0), errMetric.cxNorm);
        }
    }
}

bool NOMAD::TRIPMSolver::computeSlackMultipliers(SGTELIB::Matrix& slackMultipliers,
                                                 const SGTELIB::Matrix& XS,
                                                 const SGTELIB::Matrix& Jx,
                                                 const SGTELIB::Matrix& Gx,
                                                 const double mu)
{
    const int nbVar = XS.get_nb_rows();
    const int nbCons = Jx.get_nb_rows();
    const int n = Jx.get_nb_cols();

    // Compute matrices W and bls defined as:
    // 1- W = [ grad c1(x) ... grad cm(x) S ], where S = diag(s),
    //    with s the slack variables of the problem, and m the number of constraints.
    //    NB: grad c1(x) ... grad cm(x) are accessible with Jx = [ grad c1(x) ... grad cm(x) ]'
    // 2- bls =  [ grad f(x) ]
    //           [ -mu e     ]
    //    with e = [1, ..., 1] m times.
    SGTELIB::Matrix W("W", nbVar, nbCons);
    SGTELIB::Matrix bls("bls", nbVar, 1);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nbCons; j++)
        {
            W.set(i, j, Jx.get(j, i));
        }
        bls.set(i, 0, Gx.get(i, 0));
    }
    for (int i = 0; i < nbCons; i++)
    {
        for (int j = 0; j < nbCons; j++)
        {
            if (i == j)
            {
                W.set(i + n, j, XS.get(i + n));
            }
            else
            {
                W.set(i + n, j, 0);
            }
        }
        bls.set(i + n, 0, -mu);
    }

    // slackMultipliers is the least square solution of || W y - bls ||
    slackMultipliers = SGTELIB::Matrix::solve_least_squares_SVD(W, bls);

    // Enforce sign of slackMultipliers, which must remain negative.
    for (int i = 0; i < nbCons; i++)
    {
        if (slackMultipliers.get(i, 0) >= 0)
        {
            const double si = XS.get(i + n, 0);
            slackMultipliers.set(i, 0, -std::min(std::abs(1E-3), std::abs(mu / si)));
        }
    }
    return true;
}

NOMAD::TRIPMSolver::BarrierSolverStatus NOMAD::TRIPMSolver::solveBarrierSubproblem(SGTELIB::Matrix& x,
                                                                                   SGTELIB::Matrix& XSp,
                                                                                   SGTELIB::Matrix& cslack,
                                                                                   SGTELIB::Matrix& lambda,
                                                                                   SGTELIB::Matrix& Gx,
                                                                                   SGTELIB::Matrix& cons,
                                                                                   SGTELIB::Matrix& Jx,
                                                                                   NOMAD::TRIPMSolver::TRIPMErrorMetric& errMetric,
                                                                                   const SGTELIB::Matrix& QPModel,
                                                                                   const SGTELIB::Matrix& XS,
                                                                                   const SGTELIB::Matrix& lvar,
                                                                                   const SGTELIB::Matrix& uvar,
                                                                                   const double mu,
                                                                                   const double atolOpt,
                                                                                   const double atolFeas,
                                                                                   const int verbose_degree) const
{
    auto isPointStrictlyFeasible = [](const SGTELIB::Matrix& xs,
                                      const SGTELIB::Matrix& lvar,
                                      const SGTELIB::Matrix& uvar,
                                      const int n) -> bool
    {
        for (int i = 0; i < n; ++i)
        {
            if ((xs.get(i, 0) <= lvar.get(i, 0)) || (xs.get(i, 0) >= uvar.get(i, 0)))
            {
                std::string err = "TRIPMSolver::solveBarrierSubproblem error: ";
                err += "current iterate x is not strictly feasible for index i = " + std::to_string(n);
                std::printf("%s: lb[i] = %f, ub[i] = %f, x[i] = %f\n", err.c_str(),
                            lvar.get(i, 0), uvar.get(i, 0), xs.get(i, 0));
                return false;
            }
        }
        return true;
    };

    auto computeQuadFctValue = [](const double g0,
                                  const SGTELIB::Matrix& g,
                                  const SGTELIB::Matrix& H,
                                  const SGTELIB::Matrix& x) -> double
    {
        const int n = x.get_nb_rows();

        // NB: Do not use SGTELIB matrix operations for faster computation
        double quadVal = g0;
        for (int i = 0; i < n; ++i)
        {
            quadVal += g.get(i,0) * x.get(i,0);
            double phxi = 0;
            for (int j = 0 ; j < n ; j++)
            {
                phxi += H.get(i,j) * x.get(j,0);
            }
            quadVal += 0.5 * x.get(i,0) * phxi;
        }
        return quadVal;
    };

    const bool verbose = verbose_degree > 0;

    // Slack variables are used
    const int n = x.get_nb_rows();
    const int nbCons = QPModel.get_nb_rows() - 1;
    const int nbVar = n + nbCons;

    // Initialize X
    XSp = XS;
    if (!isPointStrictlyFeasible(XSp, lvar, uvar, n))
    {
        return BarrierSolverStatus::NUM_ERROR;
    }

    // Allocation of matrices and vectors for solver_barrier
    SGTELIB::Matrix Xcan("Xcan", nbVar, 1);
    SGTELIB::Matrix xp("xp", n, 1);
    SGTELIB::Matrix r("r", nbCons, 1); // The residual
    double nr = 0;
    SGTELIB::Matrix checkCons("checkCons", nbCons, 1);
    SGTELIB::Matrix checkSlack("checkSlack", nbCons, 1);

    // Allocation of matrices for the normal step
    SGTELIB::Matrix W("W", nbCons, nbCons + n);
    SGTELIB::Matrix wq("wq", nbCons + n, 1);
    SGTELIB::Matrix vxs("vxs", nbCons + n, 1); // The normal step
    double nv = 0;
    SGTELIB::Matrix zer("zer", nbCons + n, 1);
    zer.fill(0.0);

    // Compute p
    SGTELIB::Matrix p("p", nbVar, 1);
    p.fill(1e15); // Prevent the activation of the stopping criterion based on || p ||
    double np = p.norm();
    SGTELIB::Matrix Q("Q", nbCons + n, nbCons + n);
    Q.fill(0);
    SGTELIB::Matrix qc("qc", nbCons + n, 1);
    SGTELIB::Matrix HLag("HLag", n, n);

    // Allocation of matrix for the second order correction (soc) step
    SGTELIB::Matrix soc_step("soc_step", nbVar, 1);
    SGTELIB::Matrix consSOC("consSOC", nbCons, 1);
    SGTELIB::Matrix cslackSOC("cslackSOC", nbCons, 1);
    SGTELIB::Matrix vxsSOC("vxsSOC", nbVar, 1);

    // Barrier solver parameters:
    constexpr double epsilon_1 = 1e-8; // trust-region successful ratio
    constexpr double epsilon_2 = 0.9; // trust-region very successful ratio
    constexpr double gamma_1 = 0.5; // trust-region decrease factor
    constexpr double gamma_2 = 2; // trust-region increase factor
    constexpr double tol_TR_radius = 1e-8; // Below this value, declare success.

    // tolerances on mu
    const double tol_mu_opt = std::max(mu, atolOpt - mu);
    const double tol_mu_feas = std::max(mu, atolFeas);

    // Tolerance bounds
    constexpr double tol_bounds = 1e-8;

    // Merit function parameters
    double nu = 1;

    double Delta = 1; // trust-region initial radius
    constexpr double smallestDelta = 1e-15;
    constexpr double largestDelta = 1e15;

    constexpr double tau = 0.995;
    // const double rho = 0.5;
    constexpr double rho = 0.1;

    constexpr double Delta_normal_step_factor = 0.8; // Factor of Delta used in normal step
    constexpr double small_p = 1e-8; // Below this value `p` is considered 0 (declare success).

    // Compute cons, Jx and Gx
    QPModelUtils::getModelObjGrad(Gx, QPModel, x);
    QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
    QPModelUtils::getModelCons(cons, QPModel, x);
    for (int i = 0; i < nbCons; ++i)
    {
        cslack.set(i, 0, cons.get(i, 0) + XSp.get(i + n, 0));
    }

    // Stopping criteria for the barrier solver
    // 1- Inner success conditions
    computeSlackMultipliers(lambda, XSp, Jx, Gx, mu);
    computeErrorFunctionMetric(errMetric, XSp, QPModel, lvar, uvar, lambda, mu, true);

    // Initial logging
    if (verbose)
    {
        std::printf("\nBarrier subproblem algorithm\n");
        std::printf("Number of variables: %d\n", n);
        std::printf("Number of inequality constraints: %d\n", nbCons);
        std::printf("Maximum number of iterations allowed for inner loop: %zu\n", max_iter_inner);
        std::printf("Initial trust-region radius: %e\n", Delta);
        std::printf("mu parameter: %e\n", mu);
        std::printf("tol_mu_opt parameter: %e\n", tol_mu_opt);
        std::printf("tol_mu_feas parameter: %e\n\n", tol_mu_feas);

        std::printf("%5s %s %9s %27s %15s %13s %13s %10s %19s %17s %17s %10s %13s %13s %16s\n",
                    "iter (inner)", "nb successive unsuccessful iter", "f(x)", "|| max(0, c(x)) ||_inf",
                    "|| c(x) + s ||", "|| x - P[x - grad L(x)] ||", "|| S lambda + mu ||",
                    "nu", "norm residual",
                    "|| lambda ||_inf", "|| x - xp ||_inf", "|| p ||", "|| v ||", "TR ratio", "TR radius");
        constexpr int maxLineWidth = 270;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    double ared =0 , pred = 0;
    size_t successiveUnsuccessful = 0;
    double distXInnerLoop = NOMAD::INF;

    auto status = BarrierSolverStatus::FAILURE;
    for (int iter = 0; iter < (int)max_iter_inner; ++iter)
    {
        if (verbose)
        {
            const double objVal = QPModelUtils::getModelObj(QPModel, x);
            const double cx = cslack.norm();
            double normMaxCx0 = 0;
            for (int i = 0; i < n; ++i)
            {
                normMaxCx0 = std::max(cons.get(i, 0), normMaxCx0);
            }

            if (pred != 0)
            {
                std::printf(" %-12d %-10zu %+33e %16e %21e %21e %22e %18e %14e %14e %18e %15e %13e %13e %15e\n",
                            iter, successiveUnsuccessful, objVal, normMaxCx0, cx,
                            errMetric.projlagGradNorm, errMetric.slackLambdaMuNorm, nu, nr,
                            lambda.norm_inf(), distXInnerLoop, np, nv, ared/pred, Delta);
            }
            else
            {
                std::printf(" %-12d %-10zu %+33e %16e %21e %21e %22e %18e %14e %14e %18e %15e %13e %13e/0 %13e\n",
                            iter, successiveUnsuccessful, objVal, normMaxCx0, cx,
                            errMetric.projlagGradNorm, errMetric.slackLambdaMuNorm, nu, nr,
                            lambda.norm_inf(), distXInnerLoop, np, nv, ared, Delta);
            }
        }

        // This stopping criterion is adapted from:
        //
        //  "An interior algorithm for nonlinear optimization that combines line search and trust region steps"
        //  by R. Waltz, J. Morales, J. Nocedal, and D. Orban
        //  Math. Program. 107, 391–408 (2006). https://doi.org/10.1007/s10107-004-0560-5
        //
        const bool innerSuccess = (errMetric.projlagGradNorm <= std::max(1.0, errMetric.projObjGrad) * tol_mu_opt &&
                                   errMetric.slackLambdaMuNorm <= std::max(1.0, errMetric.projObjGrad) * tol_mu_opt &&
                                   errMetric.cxNorm <= std::max(1.0, errMetric.cxInitNorm) * tol_mu_feas)
                                  || np <= small_p || Delta < tol_TR_radius;
        if (innerSuccess)
        {
            status = BarrierSolverStatus::SOLVED;
            break;
        }

        const bool innerFailure = (distXInnerLoop <= tol_dist_successive_x);
        if (innerFailure)
        {
            status = BarrierSolverStatus::STAGNATION_ITERATES;
            break;
        }

        // Compute the normal step vxs := (vx, vs), i.e., solve:
        //  min (|| W vxs + clack ||_2)^2
        //  vxs
        //  s.t. || vxs || <= delta_normal_step
        //       vs        >= -tau/2
        //       tau/2 * l <= vx + tau/2 * x <= tau/2 * u
        //
        // W = [ Jx  S ] with S = diag(s), where s are the slack variables of the inequality constraints.
        //
        // There is no need to recompute the matrix W when the iteration has failed, as XS remains the same.
        //
        // NB: the quadratic function is convex, i.e. W^t W is at least semi-definite positive
        if (successiveUnsuccessful == 0)
        {
            for (int i = 0; i < nbCons; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    W.set(i, j, Jx.get(i, j));
                }
                for (int j = 0; j < nbCons; j++)
                {
                    if (i == j)
                    {
                        W.set(i, j + n, XSp.get(i + n, 0));
                    }
                    else
                    {
                        W.set(i, j + n, 0.0);
                    }
                }
            }
        }
        vxs.fill(0.0);
        NOMAD::DoglegTRSolver::solve(vxs, W, cslack, Delta_normal_step_factor * Delta);

        for (int i = 0; i < n; i++)
        {
            Xcan.set(i, 0, XSp.get(i, 0) + vxs.get(i, 0));
            xp[i] = x[i] + vxs.get(i, 0);
        }
        for (int i = 0; i < nbCons; i++)
        {
            Xcan.set(i + n, 0, XSp.get(i + n, 0) + vxs.get(i + n, 0));
        }
        QPModelUtils::getModelCons(checkCons, QPModel, xp);
        for (int j = 0; j < nbCons; ++j)
        {
            checkSlack.set(j, 0, XSp.get(j + n, 0) + checkCons.get(j, 0));
        }

        // Backtrack to satisfy the inequality vs >= -tau/2
        double backtrackStepsize = 1;
        for (int i = 0; i < nbCons; i++)
        {
            const double vsi = vxs.get(i + n, 0);
            if (vsi < -tau / 2.0)
            {
                backtrackStepsize = std::min(backtrackStepsize, - tau / (2.0 * vsi));
            }
        }
        // Backtrack to satisfy vx + tau/2 * x >= tau/2 * l and tau/2 * u >= vx + tau/2 * x
        for (int i = 0; i < n; i++)
        {
            const double vxi = vxs.get(i, 0);
            const double xi = XSp.get(i, 0);
            const double li = lvar.get(i, 0);
            const double ui = uvar.get(i, 0);
            if (vxi == 0)
            {
                if (xi < li)
                {
                    std::printf("Warning: xi = %e < li = %e for index %d\n", xi, li, i);
                }
                if (xi > ui)
                {
                    std::printf("Warning: xi = %e > ui = %e for index %d\n", xi, ui, i);
                }
            }
            else
            {
                if (vxi < tau * (li - xi))
                {
                    backtrackStepsize = std::min(backtrackStepsize, tau * (li - xi) / (2.0 * vxi));
                }
                if (vxi > tau * (ui - xi))
                {
                    backtrackStepsize = std::min(backtrackStepsize, tau * (ui - xi) / (2.0 * vxi));
                }
            }
        }
        vxs.multiply(backtrackStepsize);
        if (verbose)
        {
            nv = vxs.norm();
        }

        // Compute p:= (px, ptildes) with projected CG, i.e. solve:
        // min (1/2) p' Q p + qc' p
        //  p
        // s.t. W p = W v
        //      || p || <= Delta
        //      ptildes + tau >= 0,
        //      tau * l <= px + tau * x <= tau * u
        //
        // where:
        //
        // Q = [ HLag + mu * (diag(1/(x-l)^2) + diag(1/(u-x)^2) ]
        //     [           - S Sigma S                          ]
        //
        // qc = [ grad f - mu (diag(1/(x-l)) - diag(1/(u-x))) ]
        //      [         -mu e                               ]
        // where Sigma = S^-1 Y with Y the estimated Lagrange multipliers (always negative)
        //
        // As before, there is no need to compute the Q and qc matrices when the iteration has failed,
        // as XSp does not change.
        if (successiveUnsuccessful == 0)
        {
            QPModelUtils::getModelLagrangianHessian(HLag, QPModel, x, lambda);
            for (int i = 0; i < n; i++)
            {
                const double xi = XSp.get(i, 0);
                const double li = lvar.get(i, 0);
                const double ui = uvar.get(i, 0);
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        Q.set(i, j, HLag.get(i, j) + mu * (1 / std::pow(xi - li, 2) + 1 / std::pow(ui - xi, 2)));
                    }
                    else
                    {
                        Q.set(i, j, HLag.get(i, j));
                    }
                }
                qc.set(i, 0, Gx.get(i, 0) - mu / (xi - li) + mu / (ui - xi));
            }
            for (int i = 0; i < nbCons; i++)
            {
                const double si = XSp.get(i + n, 0);
                const double lambda_i = lambda.get(i, 0);
                Q.set(i + n, i + n, - lambda_i * si); // NB: S Sigma S = S (S^-1 lambda) S = lambda S
                qc.set(i + n, 0, -mu);
            }
        }

        // Compute the residual r:= W v
        SGTELIB::Matrix::inplace_product(r, W, vxs);
        p = vxs; // This starting point is supposed to be feasible and below the threshold
        p.set_name("p");
        const auto solverStatus = NOMAD::ProjectedConjugateGradientSolver::solve(p, Q, qc, W, r, Delta, verbose_degree - 1 > 0);

        // Special case: when there is a NO_INIT_SOLUTION, it means that the initial point solution of:
        // min || p ||
        // s.t W p = W v
        // is above Delta.
        // In this case, we only try a vertical step, i.e.
        // p = vxs; still the backtracking needs to be done to be sure the point satisfies the bounds.
        if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::NO_INIT_SOLUTION)
        {
            p = vxs;
            p.set_name("p");
        }

        if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::QUAD_ROOTS_ERROR ||
            solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::TR_PARAM_ERROR ||
            solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::MATRIX_DIMENSIONS_FAILURE)
        {

            std::printf("TRIPMSolver::solveBarrierSubproblem error: Projected Conjugate Gradient has failed.\n");
            std::printf("Status:\n");
            if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::QUAD_ROOTS_ERROR)
            {
                std::printf("Quadratic resolution root error\n");
            }
            else
            {
                std::printf("Wrong parameters\n");
            }
            return BarrierSolverStatus::NUM_ERROR;
        }

        // Backtrack to satisfy ptildes >= -tau e
        backtrackStepsize = p.norm() > Delta ? Delta / p.norm() : 1.0;
        for (int i = 0; i < nbCons; i++)
        {
            const double psi = p.get(i + n, 0);
            if (psi < -tau)
            {
                backtrackStepsize = std::min(backtrackStepsize, - tau / psi);
            }
        }
        // Backtrack to satisfy px + tau * x >= tau * l and tau * u >= px + tau * x
        for (int i = 0; i < n; i++)
        {
            const double pxi = p.get(i, 0);
            const double xi = XSp.get(i, 0);
            const double li = lvar.get(i, 0);
            const double ui = uvar.get(i, 0);
            if (pxi == 0)
            {
                if (xi < li)
                {
                    std::printf("Warning: xi = %e < li = %e for index %d\n", xi, li, i);
                }
                if (xi > ui)
                {
                    std::printf("Warning: xi = %e > ui = %e for index %d\n", xi, ui, i);
                }
            }
            else
            {
                if (pxi < tau * (li - xi))
                {
                    backtrackStepsize = std::min(backtrackStepsize, tau * (li - xi) / pxi);
                }
                if (pxi > tau * (ui - xi))
                {
                    backtrackStepsize = std::min(backtrackStepsize, tau * (ui - xi) / pxi);
                }
            }
        }
        p.multiply(backtrackStepsize);

        // Compute the residual: r = J(x) * px + ps + (c(x) + s)
        // (since ps = S * pstildes, W p = Jx * px + ps).
        SGTELIB::Matrix::inplace_product(r, W, p);
        r.add(cslack);
        nr = r.norm();

        // Update nu_l
        double feasibilityGap = cslack.norm() - nr + tol_mu_opt;
        if (feasibilityGap < 0)
        {
            // We try a feasibility step: i.e. p = vxs
            nr = checkSlack.norm();
            feasibilityGap = cslack.norm() - nr;
            if (feasibilityGap <= 0)
            {
                nu = NOMAD::INF;
            }
            else
            { // feasibilityGap > 0
                p = vxs;
                nu = 1e15;
            }
        }
        else if (feasibilityGap > 0)
        {
            const double nu_trial = computeQuadFctValue(0, qc, Q, p) / ((1 - rho) * feasibilityGap);
            nu = nu >= nu_trial ? nu : nu_trial + 1.0;
            // nu = std::max(nu, -nu); // std::max(nu, smallest_nu);
        }

        // Compute pred (> 0):
        pred = nu * feasibilityGap - computeQuadFctValue(0, qc, Q, p);

        // Compute ared (> 0)
        // a- Perform the following change of variables: ps := S ptildes:
        for (int i = 0; i < nbCons; i++)
        {
            const double si = XSp.get(i + n, 0);
            p.set(i + n, 0, p.get(i + n, 0) * si);
        }
        np = p.norm();

        // b- Update candidates
        for (int i = 0; i < n; i++)
        {
            Xcan.set(i, 0, XSp.get(i, 0) + p.get(i, 0));
            xp[i] = x[i] + p.get(i, 0);
        }

        // c- Make sure the problem is strictly feasible, relatively to the bound constraints
        for (int i = 0; i < n; ++i)
        {
            const double li = lvar.get(i, 0);
            const double ui = uvar.get(i, 0);
            double xican = std::max(Xcan.get(i, 0), li + tol_bounds);
            xican = std::min(xican, ui - tol_bounds);
            Xcan.set(i, 0, xican);

            double xi = std::max(xp.get(i, 0), li + tol_bounds);
            xi = std::min(xi, ui - tol_bounds);
            xp.set(i, 0, xi);
        }

        // Before the trial point is tested for acceptance by the merit function,
        // set the slack variables for which -ci(Xp) > sip to: sip := - ci(Xp),
        // where sp = s + ps.
        // It enables to decrease further the barrier term.
        // This slack reset strategy is taken from:
        //
        // "An interior algorithm for nonlinear optimization that
        //  combines line search and trust region steps"
        // by R.A. Waltz, J.L. Morales, J. Nocedal and D. Orban
        // Mathematical Programming 107, 391–408 (2006)
        // https://doi.org/10.1007/s10107-004-0560-5
        //
        QPModelUtils::getModelCons(cons, QPModel, xp);
        for (int i = 0; i < nbCons; i++)
        {
            const double sip = XSp.get(i + n, 0) + p.get(i + n, 0);
            Xcan.set(i + n, 0, std::max(-cons.get(i, 0), sip));
        }

        // Compute ared (> 0)
        ared = computeMeritFctBarrier(QPModel, lvar, uvar, x, XSp, mu, nu) -
                computeMeritFctBarrier(QPModel, lvar, uvar, xp, Xcan, mu, nu);

        // Trust-region update
        if ((ared >= pred * epsilon_1) && (pred > 0))
        { // r >= epsilon_1
            // Accept x and s steps
            XSp = Xcan;
            distXInnerLoop = SGTELIB::Matrix::distNorm2(x, xp);
            x = xp;

            // Increase Delta
            if (ared >= pred * epsilon_2) // r >= epsilon_2
            {
                Delta = std::min(gamma_2 * Delta, std::max(1 / mu, largestDelta)); // std::max(d.norm(), delta);
            }

            status = BarrierSolverStatus::ONE_STEP_MADE;
            successiveUnsuccessful = 0;

            // Check optimality
            QPModelUtils::getModelCons(cons, QPModel, x);
            for (int j = 0; j < nbCons; ++j)
            {
                cslack.set(j, 0, XSp.get(j + n, 0) + cons.get(j, 0));
            }
            QPModelUtils::getModelObjGrad(Gx, QPModel, x);
            QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
            computeSlackMultipliers(lambda, XSp, Jx, Gx, mu);
            computeErrorFunctionMetric(errMetric, XSp, QPModel, lvar, uvar, lambda, mu, true);
            continue;
        }

        // Try a second order correction step
        QPModelUtils::getModelCons(consSOC, QPModel, xp);

        // 1- Compute c(x + px) + s + ps
        for (int j = 0; j < nbCons; ++j)
        {
            cslackSOC.set(j, 0, Xcan.get(j + n, 0) + consSOC.get(j, 0));
        }

        // 2- Get vxs: we need to perform the following change of variables, i.e. vs := S vs.
        // (p has already been updated)
        for (int i = 0; i < n; ++i)
        {
            vxsSOC.set(i, 0, vxs.get(i, 0));
        }
        for (int i = 0; i < nbCons; ++i)
        {
            const double si = XSp.get(i + n, 0);
            const double vsi = vxs.get(i + n, 0);
            vxsSOC.set(i + n, 0, si * vsi);
        }

        // Compute the second order correction step
        const bool isSecondOrderStepComputed = computeSecondOrderCorrectionStep(soc_step, XSp, cslackSOC, Jx, p, vxsSOC);
        if (isSecondOrderStepComputed)
        {
            // Check that ps + soc_step_s >= -tau s
            bool isInFeasibleRegion = true;
            for (int j = 0; j < nbCons; ++j)
            {
                const double psi = p.get(j + n, 0);
                const double si = XSp.get(j + n, 0);
                const double ysi = soc_step.get(j + n, 0);
                if ((psi + ysi) < (- tau * si))
                {
                    isInFeasibleRegion = false;
                    break;
                }
            }

            // Check that px + soc_step_x + tau * x >= tau * l and tau * u >= px + soc_step_x + tau * x
            if (isInFeasibleRegion)
            {
                for (int i = 0; i < n; ++i)
                {
                    const double pxi = p.get(i, 0);
                    const double li = lvar.get(i, 0);
                    const double ui = uvar.get(i, 0);
                    const double xi = XSp.get(i, 0);
                    const double yxi = soc_step.get(i, 0);
                    if ((pxi + yxi + tau * xi < tau * li) || (tau * ui < pxi + yxi + tau * xi))
                    {
                        isInFeasibleRegion = false;
                        break;
                    }
                }
            }

            if (isInFeasibleRegion)
            {
                // Try the trust-region test with new candidate Xcan := Xcan + soc_step
                for (int i = 0; i < n; i++)
                {
                    Xcan.set(i, 0, XSp.get(i, 0) + soc_step.get(i, 0));
                    xp[i] = x[i] + soc_step.get(i, 0);
                }

                // Make sure the problem is strictly feasible, relatively to the bound constraints
                for (int i = 0; i < n; ++i)
                {
                    const double li = lvar.get(i, 0);
                    const double ui = uvar.get(i, 0);
                    double xican = std::max(Xcan.get(i, 0), li + tol_bounds);
                    xican = std::min(xican, ui - tol_bounds);
                    Xcan.set(i, 0, xican);

                    double xi = std::max(xp.get(i, 0), li + tol_bounds);
                    xi = std::min(xi, ui - tol_bounds);
                    xp.set(i, 0, xi);
                }

                // Apply the magic step
                QPModelUtils::getModelCons(cons, QPModel, xp);
                for (int i = 0; i < nbCons; i++)
                {
                    const double sip = XSp.get(i + n, 0) + soc_step.get(i + n, 0) + p.get(i + n, 0);
                    Xcan.set(i + n, 0, std::max(-cons.get(i, 0), sip));
                }

                // Compute ared (> 0)
                ared = computeMeritFctBarrier(QPModel, lvar, uvar, x, XSp, mu, nu) -
                       computeMeritFctBarrier(QPModel, lvar, uvar, xp, Xcan, mu, nu);

                if (ared >= pred * epsilon_1)
                {
                    // Accept x and s steps
                    XSp = Xcan;
                    distXInnerLoop = SGTELIB::Matrix::distNorm2(x, xp);
                    x = xp;

                    status = BarrierSolverStatus::ONE_STEP_MADE;
                    successiveUnsuccessful = 0;

                    // Check optimality
                    QPModelUtils::getModelCons(cons, QPModel, x);
                    for (int j = 0; j < nbCons; ++j) {
                        cslack.set(j, 0, XSp.get(j + n, 0) + cons.get(j, 0));
                    }
                    QPModelUtils::getModelObjGrad(Gx, QPModel, x);
                    QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
                    computeSlackMultipliers(lambda, XSp, Jx, Gx, mu);
                    computeErrorFunctionMetric(errMetric, XSp, QPModel, lvar, uvar, lambda, mu, true);
                    continue;
                }
            }
        }

        // Decrease trust-region radius
        if (pred < 0)
        {
            //verbose && std::cout << "quad = " << getModelObj(p, Q, qc) << " <= " << getModelObj(zer, Q, qc);
            //verbose && std::cout << " |Wscalp| = " << SGTELIB::Matrix::product(Wscal, p).norm() << " |Wp| = " << SGTELIB::Matrix::product(W, p).norm() << " |Wvxs| = " << SGTELIB::Matrix::product(W, vxs).norm();
            //verbose && std::cout << " |r| = " << nr << " |rv| = " << SGTELIB::Matrix::add(SGTELIB::Matrix::product(W, vxs), cslack).norm();
            //verbose && std::cout << " |c+s| = " << cslack.norm() << " nu = " << nu << std::endl;
            //verbose && std::cout << "pred = " << nu * cslack.norm();
            //verbose && std::cout << " + " << - getModelObj(p, Q, qc);
            //verbose && std::cout << " + " << - nu * nr << std::endl;
        }
        // Decrease Delta
        Delta = std::max(gamma_1 * std::min(Delta, np), smallestDelta);
        successiveUnsuccessful += 1;
    }
    if (verbose)
    {
        std::printf("\nStatus: ");
        std::printf("f(x*) = %e\n", QPModelUtils::getModelObj(QPModel, x));
        std::printf("|| c(x*) + s* || = %e\n", cslack.norm());
        double normMaxCx0 = 0;
        for (int i = 0; i < n; ++i)
        {
            normMaxCx0 = std::max(cons.get(i, 0), normMaxCx0);
        }
        std::printf("|| max(0, c(x*)) ||_inf = %e\n", normMaxCx0);
        if (status == BarrierSolverStatus::SOLVED)
        {
            std::printf("Has reached the minimum tolerance:\n");
            //std::printf("E(x,s,y;mu) = %e <= tol = %e or\n", errorInnerVal, tol_mu);
            std::printf("Trust-region radius %e below tol = %e or\n", Delta, tol_TR_radius);
            std::printf("|| p || = %e below tol = %e\n", np, small_p);
        }
        else if (status == BarrierSolverStatus::STAGNATION_ITERATES)
        {
            std::printf("Inner steps have stagnated:\n");
            std::printf("|| x - xp || = %e <= %e\n", distXInnerLoop, tol_dist_successive_x);
        }
        else if (status == BarrierSolverStatus::ONE_STEP_MADE)
        {
            std::printf("At least one step has been made\n");
        }
        else if (status == BarrierSolverStatus::FAILURE)
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

double NOMAD::TRIPMSolver::computeMeritFctBarrier(const SGTELIB::Matrix& QPModel,
                                                  const SGTELIB::Matrix& lvar,
                                                  const SGTELIB::Matrix& uvar,
                                                  const SGTELIB::Matrix& x,
                                                  const SGTELIB::Matrix& XS,
                                                  const double mu,
                                                  const double nu)
{
    const int n = x.get_nb_rows();
    const int nbVars = XS.get_nb_rows();
    const int nbCons = nbVars - n;

    const double fx = QPModelUtils::getModelObj(QPModel, x);

    double logPartVal = 0;
    for (int i = 0; i < nbCons; i++)
    {
        const double si = XS.get(i + n, 0);
        logPartVal -= mu * std::log(si);
    }

    for (int i = 0; i < n; i++)
    {
        const double xi = XS.get(i, 0);
        const double ui = uvar.get(i, 0);
        const double li = lvar.get(i, 0);
        logPartVal -= mu * std::log(xi - li);
        logPartVal -= mu * std::log(ui - xi);
    }

    SGTELIB::Matrix cons("cons", nbCons, 1);
    QPModelUtils::getModelCons(cons, QPModel, x);
    SGTELIB::Matrix cslack("cslack", nbCons, 1);
    for (int j = 0; j < nbCons; ++j)
    {
        cslack.set(j, 0, XS.get(j + n, 0) + cons.get(j, 0));
    }

    return fx + logPartVal + nu * cslack.norm();
}


bool NOMAD::TRIPMSolver::computeSecondOrderCorrectionStep(SGTELIB::Matrix& y,
                                                          const SGTELIB::Matrix& XS,
                                                          const SGTELIB::Matrix& cslackXcan,
                                                          const SGTELIB::Matrix& Jx,
                                                          const SGTELIB::Matrix& pxs,
                                                          const SGTELIB::Matrix& vxs)
{
    // In procedure SOC, the second order correction step is computed only when
    // || tangent_step || <= 0.1 || vxs ||
    // || pxs || = || vxs + tangent_step || <= || vxs || + || tangent_step ||
    // so we apply this step, when || pxs || <= 1.1 || vxs ||
    constexpr double ratio_p_v_trigger = 1.1;
    if (pxs.norm() > ratio_p_v_trigger * vxs.norm())
    {
        return false;
    }

    const int nbVar = XS.get_nb_rows();
    const int nbCons = Jx.get_nb_rows();
    const int n = Jx.get_nb_cols();

    // Compute matrix W defined as:
    // W = [Jx S]
    SGTELIB::Matrix W("W", nbCons, nbVar);
    for (int i = 0; i < nbCons; i++)
    {
        for (int j = 0; j < n; j++)
        {
            W.set(i, j, Jx.get(i, j));
        }
        for (int j = 0; j < nbCons; j++)
        {
            if (i == j)
            {
                W.set(i, j + n, XS.get(i + n, 0));
            }
            else
            {
                W.set(i, j + n, 0.0);
            }
        }
    }

    // y := W' (W W')^(-1) ( cons(x + px) + s + ps ) = W' (W W')^(-1) cslackXcan
    // Minimize || W y - cslackXcan ||^2 with QR decomposition.
    // 1- Compute the QR factorization of W if W has more rows than columns, otherwise
    // the computation of W^T.
    const int nrows = W.get_nb_rows();
    const int ncols = W.get_nb_cols();
    auto Q = new double*[std::max(nrows, ncols)];
    auto R = new double*[std::max(nrows, ncols)];
    auto M = new double*[std::max(nrows, ncols)];
    for (int i = 0; i < std::max(nrows, ncols); ++i)
    {
        Q[i] = new double[std::max(nrows, ncols)];
        R[i] = new double[std::min(nrows, ncols)];
        M[i] = new double[std::min(nrows, ncols)];
    }
    for (int i = 0; i < std::max(nrows, ncols); ++i)
    {
        for (int j = 0; j < std::min(nrows, ncols); ++j)
        {
            if (nrows >= ncols)
            {
                // M := W
                M[i][j] = W.get(i, j);
            }
            else
            {
                // M := W^T
                M[i][j] = W.get(j, i);
            }
        }
    }
    std::string error_str;
    bool factorization_success = qr_factorization(error_str, M, Q, R, std::max(nrows, ncols), std::min(nrows, ncols));

    if (!factorization_success)
    {
        for (int i = 0; i < std::max(nrows, ncols); ++i)
        {
            delete [] Q[i];
            delete [] R[i];
            delete [] M[i];
        }
        delete [] Q;
        delete [] R;
        delete [] M;
        return false;
    }

    auto c = new double[nrows];
    if (nrows >= ncols)
    {
        // The solution is the least-square solution of || W y - cslackXcan ||^2
        // First, compute c := Q^T cslackXcan
        for (int i = 0; i < nrows; ++i)
        {
            c[i] = 0;
            for (int j = 0; j < nrows; ++j)
            {
                const double bi = cslackXcan.get(i, 0);
                c[i] += Q[j][i] * bi;
            }
        }

        // Then compute y by solving R[1:n] y = c with back-substitution
        y.set(ncols - 1, 0, c[ncols - 1] / R[ncols - 1][ncols - 1]);
        for (int i = ncols - 2; i > -1; --i)
        {
            double yi = c[i];
            for (int j = i + 1; j < ncols; ++j)
            {
                yi -= y.get(j, 0) * R[i][j];
            }
            yi /= R[i][i];
            y.set(i, 0, yi);
        }
    }
    else
    {
        // The solution is the least-norm solution of || W y - cslackXcan ||^2
        // First, compute c := R[1:m, 1:m]^-T cslackXcan with forward substitution
        c[0] = cslackXcan.get(0, 0) / R[0][0];
        for (int i = 1; i < nrows; ++i)
        {
            c[i] = cslackXcan.get(i, 0);
            for (int j = 0; j < i; ++j)
            {
                c[i] -= c[j] * R[j][i];
            }
            c[i] /= R[i][i];
        }

        // Then compute y := Q[:, 1:nrows] c
        for (int i = 0; i < ncols; ++i)
        {
            double yi = 0;
            for (int j = 0; j < nrows; ++j)
            {
                yi += Q[i][j] * c[j];
            }
            y.set(i, 0, yi);
        }
    }

    delete[] c;
    for (int i = 0; i < std::max(nrows, ncols); ++i)
    {
        delete [] Q[i];
        delete [] R[i];
        delete [] M[i];
    }
    delete [] Q;
    delete [] R;
    delete [] M;

    return true;
}
