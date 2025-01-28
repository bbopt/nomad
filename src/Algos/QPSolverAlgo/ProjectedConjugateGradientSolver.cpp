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
 \file   ProjectedConjugateGradientSolver.cpp
 \brief  Projected Conjugate Gradient algorithm: implementation
 \author Tangi Migot and Ludovic Salomon
 \see    ProjectedConjugateGradientSolver.hpp
 */
#include "ProjectedConjugateGradientSolver.hpp"

#include <cstdio>

#include "../../Math/MathUtils.hpp"
#include "../../Math/MatrixUtils.hpp"

NOMAD::ProjectedConjugateGradientSolverStatus NOMAD::ProjectedConjugateGradientSolver::solve(SGTELIB::Matrix& x,
                                                                                             const SGTELIB::Matrix& G,
                                                                                             const SGTELIB::Matrix& c,
                                                                                             const SGTELIB::Matrix& A,
                                                                                             const SGTELIB::Matrix& b,
                                                                                             const double delta,
                                                                                             const bool verbose)
{
    // NB: m should be less than n
    const int n = x.get_nb_rows();
    const int m = b.get_nb_rows();

    constexpr double tolTR = 1e-8;
    if (delta < tolTR)
    {
        return NOMAD::ProjectedConjugateGradientSolverStatus::TR_PARAM_ERROR;
    }

    if (n <= 0 || m <= 0)
    {
        return NOMAD::ProjectedConjugateGradientSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }
    if (A.get_nb_rows() != m || A.get_nb_cols() != n)
    {
        return NOMAD::ProjectedConjugateGradientSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }
    if (G.get_nb_rows() != n || G.get_nb_cols() != n)
    {
        return NOMAD::ProjectedConjugateGradientSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }
    if (c.get_nb_rows() != n)
    {
        return NOMAD::ProjectedConjugateGradientSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }

    // NB : Other choices are possible for H
    SGTELIB::Matrix H = SGTELIB::Matrix::identity(n);

    // Define the augmented system matrix
    // M = [ I  A^T ]
    //     [ A  0   ]
    const int npm = n + m;
    auto M = new double*[npm];
    for (int i = 0; i < npm; ++i)
    {
        M[i] = new double[npm];
        for (int j = 0; j < npm; ++j)
        {
            if ((i < n) & (j < n))
            {
                // Top-left block of size nxn
                M[i][j] = H.get(i, j);
            }
            else if (i < n)
            {
                // Top-right block of size nxm
                M[i][j] = A.get(j - n, i);
            }
            else if (j < n)
            {
                // Bottom-left block of size mxn
                M[i][j] = A.get(i - n, j);
            }
            else
            {
                // Bottom-right block of size mxm
                M[i][j] = 0;
            }
        }
    }

    // Compute LDL^T factorization
    std::string error_msg;
    auto pivots = new int[npm];
    NOMAD::LDLt_factorization(error_msg, M, pivots, npm);

    // Algorithm PCG may assume that an initial feasible point x0 satisfying A x0 = b
    // is provided. If it is not the case, we try to find one, by solving
    // [ I  A^T ]  [ x0 ]  = [ 0 ]
    // [ A  0   ]  [ y  ]  = [ b ]
    constexpr double tolAxbSolved = 1e-15;
    auto sol = new double[npm];
    auto rhs = new double[npm];

    SGTELIB::Matrix Ax = SGTELIB::Matrix::product(A, x);
    SGTELIB::Matrix init_ro = b;
    init_ro.sub(Ax);

    // Temporary solution
    SGTELIB::Matrix xtmp("xtmp", n, 1);
    xtmp = x;

    // Allocations of matrices for the iterative refinement procedure
    SGTELIB::Matrix rhog("rhog", n, 1);
    SGTELIB::Matrix rhow("rhow", m, 1);
    SGTELIB::Matrix wp("wp", m, 1);

    // Compute minimal norm solution:
    // min || x ||
    // s.t. Ax = b
    // via the resolution of the augmented system
    for (int i = 0; i < n; ++i)
    {
        rhs[i] = 0;
        sol[i] = 0;
    }
    for (int i = n; i < npm; ++i)
    {
        rhs[i] = b.get(i - n, 0);
        sol[i] = 0;
    }
    NOMAD::ldl_solve(error_msg, M, rhs, sol, pivots, npm);
    for (int i = 0; i < n; ++i)
    {
        x.set(i, 0, sol[i]);
    }
    for (int i = n; i < npm; ++i)
    {
        wp.set(i - n, 0, sol[i]);
    }

    // For some ill-conditioned matrices, refining the solution can help.
    // Apply the iterative refinement procedure three times.
    constexpr int nb_iterative_init = 3;
    for (int iter = 0; iter < nb_iterative_init; ++iter)
    {
        if (x.norm() <= delta)
            break;

        // Compute rhog = - x - A'w+
        SGTELIB::Matrix::inplace_product(rhog, A.transpose(), wp);
        rhog.multiply(-1.0);
        rhog.sub(x);

        // Compute rhow = b - A x
        SGTELIB::Matrix::inplace_product(rhow, A, x);
        rhow.multiply(-1.0);
        rhow.add(b);

        // Solve
        // [ I A' ] [ Delta_x  ] = [ rhog ]
        // [ A 0  ] [ Delta_w+ ] = [ rhow ]
        for (int i = 0; i < n; ++i)
        {
            rhs[i] = rhog.get(i, 0);
            sol[i] = 0;
        }
        for (int i = n; i < npm; ++i)
        {
            rhs[i] = rhow.get(i - n, 0);
            sol[i] = 0;
        }
        NOMAD::ldl_solve(error_msg, M, rhs, sol, pivots, npm);

        // x = x + Delta_x
        for (int i = 0; i < n; ++i)
        {
            x.set(i, 0, x.get(i, 0) + sol[i]);
        }
        // wp = w + Delta_wp
        for (int i = n; i < npm; ++i)
        {
            wp.set(i - n, 0, wp.get(i - n, 0) + sol[i]);
        }
    }

    // Fix initial solution
    if (x.norm() <= delta)
    {
        xtmp = x;
    }
    else
    {
        // No initial solution has been found.
        if ((init_ro.norm() > tolAxbSolved) || xtmp.norm() > delta)
        {
            for (int i = 0; i < npm; ++i) {
                delete[] M[i];
            }
            delete[] M;
            delete[] sol;
            delete[] rhs;

            return NOMAD::ProjectedConjugateGradientSolverStatus::NO_INIT_SOLUTION;
        }
        x = xtmp;
    }

    // r = Gx + c
    SGTELIB::Matrix r = SGTELIB::Matrix::product(G, x);
    r.add(c);

    // r = Pr
    // See Algorithm 6.2 described in:
    // On the Solution of Equality Constrained Quadratic Programming Problems Arising in Optimization
    // N.I. M. Gould, M. E. Hribar, and J. Nocedal,
    // SIAM Journal on Scientific Computing, 23, Issue 4 (2001)
    // https://doi.org/10.1137/S1064827598345667
    //
    // NB: this algorithm only works in the case where H = I.
    for (int i = 0 ; i < n ; ++i )
    {
        rhs[i] = r.get(i, 0);
        sol[i] = 0;
    }
    for (int i = n ; i < npm ; ++i )
    {
        rhs[i] = 0;
        sol[i] = 0;
    }
    NOMAD::ldl_solve(error_msg, M, rhs, sol, pivots, npm);
    for (int i = 0 ; i < n ; ++i)
    {
        r.set(i, 0, sol[i]);
    }

    // Init r+
    SGTELIB::Matrix rp = r;

    // g = Pr
    SGTELIB::Matrix g("g", n, 1);
    for (int i = 0 ; i < n ; ++i )
    {
        rhs[i] = r.get(i, 0);
        sol[i] = 0;
    }
    for (int i = n ; i < npm ; ++i )
    {
        rhs[i] = 0;
        sol[i] = 0;
    }
    NOMAD::ldl_solve(error_msg, M, rhs, sol, pivots, npm);
    for (int i = 0 ; i < n ; ++i)
    {
        g.set(i, 0, sol[i]);
    }

    // Init g+
    SGTELIB::Matrix gp = g;

    // d = -g
    SGTELIB::Matrix d = g;
    d.multiply(-1);

    // Gd = G * d
    SGTELIB::Matrix Gd = SGTELIB::Matrix::product(G, d);

    // Stopping criteria
    const size_t max_iter = npm * 2;
    double rg = SGTELIB::Matrix::dot(r, g);
    const double tol_CG = 0.01 * std::sqrt(rg);
    double dtGd = SGTELIB::Matrix::dot(d, Gd);
    const double tol_itRefinement = 1e-12;

    auto computeCosTheta = [](const SGTELIB::Matrix& A, const SGTELIB::Matrix& g) -> double
    {
        const int m = A.get_nb_rows();
        const int n = A.get_nb_cols();

        double cosTheta = NOMAD::M_INF;
        for (int i = 0; i < m; ++i)
        {
            double cosThetaTmp = 0.0;
            double normAi = 0;
            for (int j = 0; j < n; ++j)
            {
                const double aij = A.get(i, j);
                cosThetaTmp += aij * g.get(j, 0);
                normAi += aij * aij;
            }
            normAi = std::sqrt(normAi);
            cosThetaTmp /= (normAi * g.norm());
            cosTheta = std::max(cosThetaTmp, cosTheta);
        }
        return cosTheta;
    };

    double cosTheta = computeCosTheta(A, g);

    // Given an iterate x, a direction d and a radius delta, compute theta such
    // that || x + theta d || = delta
    // Requires: || d || != 0, and || x || <= delta.
    auto toBoundary = [](const SGTELIB::Matrix& x,
                         const SGTELIB::Matrix& d,
                         const double delta) -> std::pair<double, double>
    {
        const double xtd = SGTELIB::Matrix::dot(x, d);
        const double nd2 = d.normsquare();
        const double nx2 = x.normsquare();
        const double delta2 = delta * delta;

        // Find the quadratic roots of the following problem
        // q2 theta^2 + q1 theta + q0 = 0,
        // where q2 = d'd, q1 = 2 x'd and q0 = x'x - delta^2.
        // We adopt a numerically stable algorithm (taken from Krylov.jl) but other algorithms could be used.
        const double q2 = nd2;
        const double q1 = 2 * xtd;
        const double q0 = nx2 - delta2;
        double root1, root2;
        const bool hasRealRoots = NOMAD::roots_quadratic(q2, q1, q0, root1, root2);
        if (!hasRealRoots)
            return {NOMAD::INF, NOMAD::INF};

        return {root1, root2};
    };

    if (verbose)
    {
        std::printf("\nProjected conjugate gradient algorithm\n");
        std::printf("Number of variables: %d\n", n);
        std::printf("Number of equality constraints: %d\n", m);
        std::printf("Trust-region radius: %f\n\n", delta);
        std::printf("Stopping criterion tolerance: %f\n", tol_CG);
        std::printf("Maximum number of iterations allowed: %zu\n\n", max_iter);

        std::printf("%5s %12s %25s %13s %14s %13s %18s\n", "iter", "f(x)", "Projected residual norm", "|| x ||", "d'Gd", "cos theta", "|| Ax - b ||_inf");
        constexpr int maxLineWidth = 107;
        for (int i = 0; i < maxLineWidth; ++i)
        {
            std::printf("-");
        }
        std::printf("\n");
    }

    auto solverStatus = NOMAD::ProjectedConjugateGradientSolverStatus::MAX_ITER_REACHED;
    for (int iter = 0; iter < (int)max_iter; ++iter)
    {
        if (verbose)
        {
            // Compute objective value
            SGTELIB::Matrix Gx = SGTELIB::Matrix::product(G, x);
            const double objVal = 0.5 * SGTELIB::Matrix::dot(Gx, x) + SGTELIB::Matrix::dot(c,x);

            // Compute || Ax - b ||
            SGTELIB::Matrix::inplace_product(Ax, A, x);
            init_ro = b;
            init_ro.sub(Ax);
            std::printf(" %-4d %+9e %23e %15e %13e %13e %18e\n", iter, objVal, rg, x.norm(), dtGd, cosTheta, init_ro.norm_inf());
        }

        // Detection of a negative curvature: in this case, return a solution x
        // such that ||x||_2 = delta.
        if (dtGd <= 0)
        {
            const auto roots = toBoundary(x, d, delta);
            const double theta = std::max(roots.first, roots.second);

            if (theta == NOMAD::INF)
            {
                // Failure to find || x ||_2 = delta
                solverStatus = NOMAD::ProjectedConjugateGradientSolverStatus::QUAD_ROOTS_ERROR;
            }
            else
            {
                // x := x + theta d
                x.add(SGTELIB::Matrix::product(d, theta));
                solverStatus = NOMAD::ProjectedConjugateGradientSolverStatus::NEGATIVE_CURVATURE;
            }
            break;
        }

        const double alpha = rg / dtGd;

        // x = x + alpha * d
        x.add(SGTELIB::Matrix::product(d, alpha));
        const double xnorm = x.norm();

        // Leave trust-region boundary: return a solution x such that ||x||_2 = delta
        // starting from the temporary solution xtmp
        if (xnorm > delta)
        {
            const auto roots = toBoundary(xtmp, d, delta);
            const double theta = std::max(roots.first, roots.second);
            x = xtmp;

            if (theta == NOMAD::INF)
            {
                // Failure to find || x ||_2 = delta
                solverStatus = NOMAD::ProjectedConjugateGradientSolverStatus::QUAD_ROOTS_ERROR;
            }
            else
            {
                // x := x + theta d
                x.add(SGTELIB::Matrix::product(d, theta));
                solverStatus = NOMAD::ProjectedConjugateGradientSolverStatus::BOUNDARY_REACHED;
            }
            break;
        }

        // r+ = r + alpha * G * d
        rp = r;
        Gd.multiply(alpha);
        rp.add(Gd);

        // g+ = Pr+
        for (int i = 0 ; i < n ; ++i)
        {
            rhs[i] = rp.get(i, 0);
            sol[i] = 0;
        }
        for (int i = n ; i < npm ; ++i )
        {
            rhs[i] = 0;
            sol[i] = 0;
        }
        NOMAD::ldl_solve(error_msg, M, rhs, sol, pivots, npm);
        for (int i = 0 ; i < n ; ++i)
        {
            gp.set(i, 0, sol[i]);
        }
        for (int i = n; i < npm; ++i)
        {
            wp.set(i - n, 0, sol[i]);
        }

        // Iterative refinement: see Algorithm 6.2 described in:
        //
        // On the Solution of Equality Constrained Quadratic Programming Problems Arising in Optimization
        // N.I. M. Gould, M. E. Hribar, and J. Nocedal,
        // SIAM Journal on Scientific Computing, 23, Issue 4 (2001)
        // https://doi.org/10.1137/S1064827598345667
        //
        // NB: this algorithm only works in the case where H = I.
        constexpr int nbIterRefinements = 3;
        for (int iterRefine = 0; iterRefine < nbIterRefinements; ++iterRefine)
        {
            cosTheta = computeCosTheta(A, gp);
            if (std::abs(cosTheta) <= tol_itRefinement)
            {
                break;
            }

            // Compute rhog = r+ - g+ - A'w+
            SGTELIB::Matrix::inplace_product(rhog, A.transpose(), wp);
            rhog.multiply(-1.0);
            rhog.sub(gp);
            rhog.add(rp);

            // Compute rhow = - A g+
            SGTELIB::Matrix::inplace_product(rhow, A, gp);
            rhow.multiply(-1.0);

            // Solve
            // [ I A' ] [ Delta_g+ ] = [ rhog ]
            // [ A 0  ] [ Delta_w+ ] = [ rhow ]
            for (int i = 0 ; i < n ; ++i)
            {
                rhs[i] = rhog.get(i, 0);
                sol[i] = 0;
            }
            for (int i = n ; i < npm ; ++i )
            {
                rhs[i] = rhow.get(i - n, 0);
                sol[i] = 0;
            }
            NOMAD::ldl_solve(error_msg, M, rhs, sol, pivots, npm);

            // g+ := g+ + Delta_g+
            for (int i = 0 ; i < n ; ++i)
            {
                gp.set(i, 0, gp.get(i, 0) + sol[i]);
            }
            // w+ := w+ + Delta_w+
            for (int i = n; i < npm; ++i)
            {
                wp.set(i - n, 0, wp.get(i - n, 0) + sol[i]);
            }
        }

        // Has reached the maximum tolerance
        double rgp = SGTELIB::Matrix::dot(rp, gp);
        if (rgp < tol_CG)
        {
            rg = rgp;
            solverStatus = NOMAD::ProjectedConjugateGradientSolverStatus::SOLVED;
            break;
        }

        // beta = r+'g+ / r'g
        const double beta = rgp / rg;

        // d = -(g+) + beta * d
        d.multiply(beta);
        d.sub(gp);

        g = gp;
        r = rp;
        rg = rgp;
        xtmp = x;
        SGTELIB::Matrix::inplace_product(Gd, G, d);
        dtGd = SGTELIB::Matrix::dot(d, Gd);
    }

    if (verbose)
    {
        std::printf("\nStatus: ");
        if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::SOLVED)
        {
            std::printf("Has reached the minimum tolerance (r'g = %e <= tol_CG = %e, ", rg, tol_CG);
            std::printf("|| x || = %e): solved\n", x.norm());
        }
        else if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::NEGATIVE_CURVATURE)
        {
           std::printf("Detection of a negative curvature ");
            std::printf("(|| x || = %e): stop\n", x.norm());
        }
        else if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::BOUNDARY_REACHED)
        {
            std::printf("Has reached the trust-region boundary (|| x || = %e ~ radius = %f): stop\n", x.norm(), delta);
        }
        else if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::QUAD_ROOTS_ERROR)
        {
            std::printf("Quadratic roots resolution error: stop\n");
        }
        else if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::MAX_ITER_REACHED)
        {
            std::printf("Maximal number of iterations reached ");
            std::printf("(|| x || = %e): stop\n", x.norm());
        }
    }

    for (int i = 0; i < npm; ++i)
    {
        delete [] M[i];
    }
    delete [] M;
    delete [] sol;
    delete [] rhs;

    return solverStatus;
}
