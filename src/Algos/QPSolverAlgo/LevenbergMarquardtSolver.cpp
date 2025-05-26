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
 \file   LevenbergMarquardtSolver.cpp
 \brief  Levenberg-Marquardt algorithm: implementation
 \author Tangi Migot and Ludovic Salomon
 \see    LevenbergMarquardtSolver.hpp
 */
#include "LevenbergMarquardtSolver.hpp"

#include "../../Algos/QPSolverAlgo/QPModelUtils.hpp"
#include "../../Algos/QPSolverAlgo/DoglegTRSolver.hpp"

NOMAD::LMSolverStatus NOMAD::LevenbergMarquardtSolver::solve(SGTELIB::Matrix& x,
                                                             SGTELIB::Matrix& XS,
                                                             const SGTELIB::Matrix& QPModel,
                                                             const SGTELIB::Matrix& lvar,
                                                             const SGTELIB::Matrix& uvar,
                                                             SGTELIB::Matrix& cX) const
{
    if (!checkDimensions(x, XS, QPModel, lvar, uvar, cX))
    {
        return LMSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }

    const int n = x.get_nb_rows();
    if (!checkBoundsCompatibilities(lvar, uvar, n))
    {
        return LMSolverStatus::BOUNDS_ERROR;
    }

    if (!checkBoundsTightness(lvar, uvar, n))
    {
        return LMSolverStatus::TIGHT_VAR_BOUNDS;
    }

    if (!checkStartingPointInBounds(XS, lvar, uvar, n))
    {
        return LMSolverStatus::STRICT_PT_FAILURE;
    }

    // Initialize x
    const int nbCons = QPModel.get_nb_rows() - 1;
    const int nbVar = n + nbCons;
    SGTELIB::Matrix xp("xp", n, 1);
    SGTELIB::Matrix XSp("XSp", nbVar, 1);
    for (int i = 0; i < n; ++i)
    {
        const double xi = XS.get(i, 0);
        x.set(i, 0, xi); // Reload x in case
        xp.set(i, 0, xi);
        XSp.set(i, 0, xi);
    }
    for (int i = n; i < nbVar; ++i)
    {
        const double si = XS.get(i, 0);
        XSp.set(i, 0, si);
    }
    SGTELIB::Matrix Xcan("Xcan", nbVar, 1);

    // Initialize c(x) + s
    SGTELIB::Matrix cons("cx", nbCons, 1);
    QPModelUtils::getModelCons(cons, QPModel, x);
    SGTELIB::Matrix cslack("cx+s", nbCons, 1);
    for (int j = 0; j < nbCons; ++j)
    {
        cslack.set(j, 0, XSp.get(j + n, 0) + cons.get(j, 0));
    }

    // Trust-region parameters
    constexpr double epsilon_1 = 1E-8; // trust-region successful ratio
    constexpr double epsilon_2 = 0.9; // trust-region very successful ratio
    constexpr double gamma_1 = 0.5; // trust-region decrease factor
    constexpr double gamma_2 = 2; // trust-region increase factor

    double Delta = 10000; // trust-region initial radius
    constexpr double smallestDelta = 1E-15;
    constexpr double largestDelta = 1E15;

    // Linesearch threshold
    constexpr double tau = 0.5;

    // Tolerance threshold for v
    constexpr double small_v = 1e-10;

    const bool verbose = verbose_level > 0;
    if (verbose)
    {
        std::printf("\nLevenberg-Marquardt algorithm\n");
        std::printf("Number of variables: %d\n", n);
        std::printf("Number of inequality constraints: %d\n", nbCons);
        std::printf("Stopping criterion tolerance: %f\n", feasibility_tol);
        std::printf("Maximum number of iterations allowed: %zu\n", max_iter);
        std::printf("Initial trust-region radius: %e\n\n", Delta);
    }

    // Compute first stopping criterion and check it is satisfied as soon as possible
    double errVal = cslack.norm();
    if (errVal <= feasibility_tol)
    {
        if (verbose)
        {
            std::printf("\nStatus:\n");
            std::printf("|| c(x*) + s* || = %e\n", errVal);
            std::printf("Has reached the minimum tolerance || c(x*) + s* || = %e <= tol = %e\n", errVal, feasibility_tol);
            return LMSolverStatus::SOLVED;
        }
    }
    double errValInit = errVal;

    // Initialize W: it is the Jacobian of the residual, i.e.
    // W = [ Jx Iq ], where q is the number of constraints.
    SGTELIB::Matrix Jx("Jx", nbCons, n);
    QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
    SGTELIB::Matrix W("W", nbCons, nbCons + n);
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
                W.set(i, j + n, 1.0);
            }
            else
            {
                W.set(i, j + n, 0.0);
            }
        }
    }

    // Second stopping criterion: stop as soon as possible
    SGTELIB::Matrix wq("vxs", nbCons + n, 1);
    SGTELIB::Matrix::inplace_product(wq, W.transpose(), cslack);
    const double wqNorm = wq.norm();
    if (wqNorm <= tol)
    {
        if (verbose)
        {
            std::printf("\nStatus:\n");
            std::printf("|| c(x*) + s* || = %e\n", errVal);
            std::printf("Has reached the minimum tolerance || [ Jx Iq ]' (c(x) + s) || = %e <= tol = %e\n",
                        wqNorm, tol);
            return LMSolverStatus::SOLVED;
        }
    }

    // Allocate memory for other matrices
    SGTELIB::Matrix vxs("vxs", nbVar, 1);
    vxs.fill(1e15);
    SGTELIB::Matrix r("r", nbCons, 1); // Residual
    SGTELIB::Matrix cxp("cxp", nbCons, 1);
    SGTELIB::Matrix checkslack("c(xp) + sp", nbCons, 1);
    SGTELIB::Matrix WtWr("W'W r", nbCons + n, 1);
    WtWr.fill(1e15);

    // Initial logging
    if (verbose)
    {
        std::printf("%5s %14s %16s %14s %17s %10s %14s %13s %24s\n",
                    "iter", "|| c(x) + s ||", "norm residual", "|| W ' W r ||" , "|| x - xp ||_inf", "|| v ||",
                    "TR radius", "TR ratio", "Backtrack step size");
        constexpr int maxLineWidth = 135;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    size_t successiveUnsuccessful = 0;
    double distXLoop = 1e15;
    double backtrackStepSize = 1.0;
    double ared =0, pred = 0;
    auto status = LMSolverStatus::MAX_ITER_REACHED;
    for (int iter = 0; iter < (int) max_iter; ++iter)
    {
        if (verbose)
        {
            if (pred != 0)
            {
                std::printf(" %-4d %+10e %16e %+15e %15e %15e %13e %13e %18e\n",
                            iter, errVal, r.norm(), WtWr.norm(), distXLoop, vxs.norm(), Delta, ared / pred, backtrackStepSize);
            }
            else
            {
                std::printf(" %-4d %+10e %16e %+15e %15e %15e %13e %13e/0 %16e\n",
                            iter, errVal, r.norm(), WtWr.norm(), distXLoop, vxs.norm(), Delta, ared, backtrackStepSize);
            }
        }

        const bool success = (errVal <= feasibility_tol) || (vxs.norm() <= small_v) || (WtWr.norm() <= tol);
        if (success)
        {
            status = LMSolverStatus::SOLVED;
            break;
        }

        const bool failure = distXLoop <= tol_dist_successive_x || backtrackStepSize == 0;
        if (failure)
        {
            status = LMSolverStatus::STAGNATION_ITERATES;
            break;
        }

        // Compute the normal step (vx, vs):
        NOMAD::DoglegTRSolver::solve(vxs, W, cslack, Delta);

        // Backtrack to satisfy vs >= -tau
        backtrackStepSize = 1;
        for (int i = 0; i < nbCons; i++)
        {
            const double vsi = vxs.get(i + n, 0);
            if (vsi < -tau)
            {
                backtrackStepSize = std::min(backtrackStepSize, - tau / vsi);
            }
        }

        // Backtrack to satisfy vx + tau * x >= tau * l and tau * u >= vx + tau * x
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
                    backtrackStepSize = std::min(backtrackStepSize, tau * (li - xi) / vxi);
                }
                if (vxi > tau * (ui - xi))
                {
                    backtrackStepSize = std::min(backtrackStepSize, tau * (ui - xi) / vxi);
                }
            }
        }
        vxs.multiply(backtrackStepSize);

        // Update the candidate; make sure it remains inside [lb, ub]
        const double bounds_tol = sol_be_strict ? 1e-13 : 0.0;
        for (int i = 0; i < n; i++)
        {
            double xcan = XSp.get(i, 0) + vxs.get(i, 0);
            xcan = std::min(xcan, uvar.get(i, 0) - bounds_tol);
            xcan = std::max(xcan, lvar.get(i, 0) + bounds_tol);
            Xcan.set(i, 0, xcan);

            double xpi = x.get(i, 0) + vxs.get(i, 0);
            xpi = std::min(xpi, uvar.get(i, 0) - bounds_tol);
            xpi = std::max(xpi, lvar.get(i, 0) + bounds_tol);
            xp[i] = xpi;
        }

        // Magic step: for all coordinates j for which cj(x) < 0,
        // set sj := -cj(x). Such steps allow to decrease further the
        // objective function. See
        //
        // Trust region methods (2000)
        // A.R. Conn and N.I. Gould and P.L. Toint
        // Society for Industrial and Applied Mathematics.
        //
        QPModelUtils::getModelCons(cxp, QPModel, xp);
        for (int i = 0; i < nbCons; ++i)
        {
            const double ci = cxp.get(i, 0);
            if (ci < 0)
            {
                Xcan.set(i + n, 0, -ci);
            }
            else
            {
                Xcan.set(i + n, 0, XSp.get(i + n, 0) + vxs.get(i + n, 0));
            }

            // Make sure s remains in [lvar, uvar] (not really useful for uvar)
            double si = Xcan.get(i + n);
            si = std::min(si, uvar.get(i + n, 0) - bounds_tol);
            si = std::max(si, lvar.get(i + n, 0) + bounds_tol);
            Xcan.set(i + n, 0, si);
        }

        // Re-adjust the step vxs
        for (int i = 0; i < nbVar; ++i)
        {
            vxs.set(i, 0, Xcan.get(i, 0) - XSp.get(i, 0));
        }

        // Compute trust-region ratio
        for (int j = 0; j < nbCons; ++j)
        {
            checkslack.set(j, 0, Xcan.get(j + n, 0) + cxp.get(j, 0));
        }

        // Compute the residual: r = J(x) * vx + vs + (c(x) + s)
        SGTELIB::Matrix::inplace_product(r, W, vxs);
        r.add(cslack);

        // Update W' W r (for stopping criteria)
        SGTELIB::Matrix::inplace_product(WtWr, W.transpose(), r);

        ared = cslack.norm() - checkslack.norm(); // > 0 if that works
        pred = cslack.norm() - r.norm();

        if ((ared >= pred * epsilon_1) && (pred > 0))
        { // TR ratio >= epsilon_1
            // Accept x and s steps
            XSp = Xcan;
            distXLoop = SGTELIB::Matrix::distNorm2(x, xp);
            x = xp;

            // Increase Delta
            if (ared >= pred * epsilon_2) // TR ratio >= epsilon_2
            {
                Delta = std::min(gamma_2 * Delta, largestDelta); // std::max(d.norm(), delta);
            }

            successiveUnsuccessful = 0;

            // Update parameters
            cons = cxp;
            cslack = checkslack;
            errVal = cslack.norm();
            QPModelUtils::getModelJacobianCons(Jx, QPModel, x);
            for (int i = 0; i < nbCons; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    W.set(i, j, Jx.get(i, j));
                }
            }
        }
        else
        {
            if (pred < 0)
            {
                // Oh...
            }
            // Decrease Delta
            Delta = std::max(gamma_1 * std::min(Delta, vxs.norm()), smallestDelta);
            successiveUnsuccessful += 1;
        }
    }

    if (verbose)
    {
        std::printf("\nStatus:\n");
        std::printf("|| c(x*) + s* || = %e\n", errVal);
        if (status == LMSolverStatus::SOLVED)
        {
            std::printf("Has reached the feasibility tolerance || (c(x*) + s*) || = %e <= tol = %e or\n",
                        errVal, feasibility_tol);
            std::printf("The step norm || v || = %e is below its tolerance tol = %e or\n", vxs.norm(), small_v);
            std::printf("A first order condition has been reached: || W' W r || = %e <= tol = %e\n",
                        WtWr.norm(), tol);
        }
        else if (errVal < errValInit)
        {
            std::printf("Has improved the solution\n");
            status = LMSolverStatus::IMPROVED;
        }
        else if (status == LMSolverStatus::MAX_ITER_REACHED)
        {
            std::printf("Has reached the maximum number of iterations\n");
        }
        else if (status == LMSolverStatus::STAGNATION_ITERATES)
        {
            std::printf("The algorithm has failed: \n");
            std::printf("|| x - xp || = %e below tolerance tol = %e or\n", distXLoop, tol_dist_successive_x);
            std::printf("backtrack step size value = %e is equal 0\n", backtrackStepSize);
        }
        else
        {
            std::printf("Unknown stopping criterion\n");
        }
        std::printf("\n");
    }

    if (status == LMSolverStatus::SOLVED)
    {
        XS = XSp;
        cX = cons;
    }
    else if (errVal < errValInit)
    {
        status = LMSolverStatus::IMPROVED;
        XS = XSp;
        cX = cons;
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            x[i] = XS.get(i, 0);
        }
    }
    return status;
}

bool NOMAD::LevenbergMarquardtSolver::checkDimensions(const SGTELIB::Matrix& x,
                                                      const SGTELIB::Matrix& XS,
                                                      const SGTELIB::Matrix& QPModel,
                                                      const SGTELIB::Matrix& lvar,
                                                      const SGTELIB::Matrix& uvar,
                                                      const SGTELIB::Matrix& cons)
{
    const int n = x.get_nb_rows();
    if (n != std::max(x.get_nb_rows(), x.get_nb_cols()) && (x.get_nb_cols() != 1))
    {
        std::string err = "LevenbergMarquardtSolver::solve error: x must be a column vector";
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbParams = QPModel.get_nb_cols();
    if (nbParams != (n+1) + n * (n+1) / 2)
    {
        std::string err = "LevenbergMarquardtSolver::solve error: ";
        err += "the number of params of the model nbParams = (n+1) * (n+2) / 2 = " + std::to_string(nbParams);
        err += " is not compatible with the dimension of the solution n = " + std::to_string(n);
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbCons = QPModel.get_nb_rows() - 1;
    if (nbCons < 1)
    {
        std::string err = "LevenbergMarquardtSolver::solve error: ";
        err += "the model has no constraints";
        std::printf("%s\n", err.c_str());
        return false;
    }

    if (nbCons != cons.get_nb_rows())
    {
        std::string err = "LevenbergMarquardtSolver::solve error: ";
        err += "the number of constraints of the model nbCons = " + std::to_string(nbCons);
        err += " is not compatible with the dimension of the cons vector q = " + std::to_string(cons.get_nb_rows());
        std::printf("%s\n", err.c_str());
        return false;
    }

    const int nbVar = XS.get_nb_rows();
    if (nbVar != std::max(XS.get_nb_rows(), XS.get_nb_cols()) && (XS.get_nb_cols() != 1))
    {
        std::string err = "LevenbergMarquardtSolver::solve error: XS must be a column vector";
        std::printf("%s\n", err.c_str());
        return false;
    }

    if (nbVar != n + nbCons)
    {
        std::string err = "LevenbergMarquardtSolver::solve error: ";
        err += "the dimension of vector XS (nbVar = n + nbCons = " + std::to_string(nbVar);
        err += " ) is not compatible with the dimension of the vector x (n = " + std::to_string(n);
        err += " ) and the number of constraints (nbCons = " + std::to_string(nbCons) + ")";
        std::printf("%s\n", err.c_str());
        return false;
    }

    if (nbVar != lvar.get_nb_rows() || nbVar != uvar.get_nb_rows())
    {
        std::string err = "LevenbergMarquardtSolver::solve error: bound constraints dimensions ";
        err += " (nlb = " + std::to_string(lvar.get_nb_cols()) + " and nub = " + std::to_string(uvar.get_nb_cols());
        err += " ) are not compatible with dimension of XS (n = " + std::to_string(nbVar) + ")";
        std::printf("%s\n", err.c_str());
        return false;
    }

    return true;
}

bool NOMAD::LevenbergMarquardtSolver::checkBoundsCompatibilities(const SGTELIB::Matrix& lvar,
                                                                 const SGTELIB::Matrix& uvar,
                                                                 const int n)
{
    for (int i = 0; i < n; ++i)
    {
        const bool areBoundsCompatible = lvar.get(i, 0) <= uvar.get(i, 0);
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

bool NOMAD::LevenbergMarquardtSolver::checkBoundsTightness(const SGTELIB::Matrix& lvar,
                                                           const SGTELIB::Matrix& uvar,
                                                           const int n)
{
    constexpr double atol = 1e-8;
    for (int i = 0; i < n; ++i)
    {
        if (std::fabs(uvar.get(i, 0) - lvar.get(i, 0)) <= atol)
        {
            return false;
        }
    }
    return true;
}

bool NOMAD::LevenbergMarquardtSolver::checkStartingPointInBounds(const SGTELIB::Matrix& XS,
                                                                 const SGTELIB::Matrix& lvar,
                                                                 const SGTELIB::Matrix& uvar,
                                                                 const int n) const
{
    const int nbVars = XS.get_nb_rows();
    for (int i = 0; i < n; ++i)
    {
        const double li = lvar.get(i, 0);
        const double ui = uvar.get(i, 0);
        const double xi = XS.get(i, 0);
        if ((sol_be_strict) && ((li >= xi) || ui <= xi))
        {
            std::printf("LevenbergMarquardtSolver::solve error: x is not strictly inside [lvar, uvar]:");
            std::printf(" for index i = %d, lb[i] = %e, ub[i] =  %e, x[i] = %e\n", i, li, ui, xi);
            return false;
        }
        if (li > xi || ui < xi)
        {
            std::printf("LevenbergMarquardtSolver::solve error: x is not inside [lvar, uvar]:");
            std::printf(" for index i = %d, lb[i] = %e, ub[i] =  %e, x[i] = %e\n", i, li, ui, xi);
            return false;
        }
    }
    for (int i = n; i < nbVars; ++i)
    {
        const double li = lvar.get(i, 0);
        const double ui = uvar.get(i, 0);
        const double xi = XS.get(i, 0);
        if (li > xi || ui < xi)
        {
            std::printf("LevenbergMarquardtSolver::solve error: x is not inside [lvar, uvar]:");
            std::printf(" for index i = %d, lb[i] = %e, ub[i] =  %e, x[i] = %e\n", i, li, ui, xi);
            return false;
        }
    }

    return true;
}

