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
#include "../../Math/MathUtils.hpp"
#include "../../Math/MatrixUtils.hpp"

#include "DoglegTRSolver.hpp"

NOMAD::DoglegTRSolverStatus NOMAD::DoglegTRSolver::solve(SGTELIB::Matrix& x,
                                                         const SGTELIB::Matrix& A,
                                                         const SGTELIB::Matrix& b,
                                                         const double delta)
{
    // Check delta parameter
    if (delta <= 1e-8)
    {
        return NOMAD::DoglegTRSolverStatus::TR_PARAM_ERROR;
    }

    // Check dimensions
    const int m = A.get_nb_rows();
    const int n = A.get_nb_cols();
    if (m == 0 || n == 0)
    {
        return NOMAD::DoglegTRSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }
    if (x.get_nb_rows() != n || x.get_nb_cols() != 1)
    {
        return NOMAD::DoglegTRSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }
    if (b.get_nb_rows() != m || b.get_nb_cols() != 1)
    {
        return NOMAD::DoglegTRSolverStatus::MATRIX_DIMENSIONS_FAILURE;
    }

    if (b.norm_inf() <= 1e-13)
    {
        x.fill(0);
        return NOMAD::DoglegTRSolverStatus::SOLVED;
    }

    SGTELIB::Matrix At = A.transpose();
    SGTELIB::Matrix AtA = SGTELIB::Matrix::product(At, A);

    // Compute the Cauchy point, which is obtained by minimizing the quadratic 1/2 || Ax + b ||^2
    // along the steepest descent direction, starting from x = 0, i.e.
    // xcp = - alpha A^T b
    // where alpha = || A^T b ||^2 / (b^T (A A^T)^2 b)

    // g0 is the gradient of x -> 1/2 || A x + b ||^2, at x = 0, i.e. g0 = A^T b.
    SGTELIB::Matrix g0("g0", n, 1);
    SGTELIB::Matrix::inplace_product(g0, At, b);

    // Compute alpha: alpha = || A^T b ||^2 / (b^T (A A^T)^2 b) = || A^T b ||^2 / || A g0 ||^2
    SGTELIB::Matrix Ag0 = SGTELIB::Matrix::product(A, g0);
    const double alpha = g0.normsquare() / Ag0.normsquare();

    // Compute the Cauchy point
    SGTELIB::Matrix xCp = g0;
    xCp.multiply(-alpha);

    // Compute the Newton point, which minimizes || A x + b ||^2 with QR decomposition.
    // 1- Compute the QR factorization of A if A has more rows than columns, otherwise
    // the computation of A^T.
    auto Q = new double*[std::max(m, n)];
    auto R = new double*[std::max(m, n)];
    auto M = new double*[std::max(m, n)];
    for (int i = 0; i < std::max(m, n); ++i)
    {
        Q[i] = new double[std::max(m, n)];
        R[i] = new double[std::min(m, n)];
        M[i] = new double[std::min(m, n)];
    }
    for (int i = 0; i < std::max(m, n); ++i)
    {
        for (int j = 0; j < std::min(m, n); ++j)
        {
            if (m >= n)
            {
                // M := A
                M[i][j] = A.get(i, j);
            }
            else
            {
                // M := A^T
                M[i][j] = A.get(j, i);
            }
        }
    }
    std::string error_str;
    bool factorization_success = qr_factorization(error_str, M, Q, R, std::max(m, n), std::min(m, n));

    if (!factorization_success)
    {
        for (int i = 0; i < std::max(m, n); ++i)
        {
            delete [] Q[i];
            delete [] R[i];
            delete [] M[i];
        }
        delete [] Q;
        delete [] R;
        delete [] M;
        return NOMAD::DoglegTRSolverStatus::QR_FACTORIZATION_FAILURE;
    }

    SGTELIB::Matrix xN("xN", n, 1);
    auto c = new double[m];
    if (m >= n)
    {
        // The solution is the least-square solution of || A x + b ||^2
        // First, compute c := -Q^T b
        for (int i = 0; i < m; ++i)
        {
            c[i] = 0;
            for (int j = 0; j < m; ++j)
            {
                const double bi = b.get(i, 0);
                c[i] += Q[j][i] * bi;
            }
            c[i] = -c[i];
        }

        // Then compute xN by solving R[1:n] x = c with back-substitution
        xN.set(n - 1, 0, c[n - 1] / R[n - 1][n - 1]);
        for (int i = n - 2; i > -1; --i)
        {
            double xi = c[i];
            for (int j = i + 1; j < n; ++j)
            {
                xi -= xN[j] * R[i][j];
            }
            xi /= R[i][i];
            xN.set(i, 0, xi);
        }
    }
    else
    {
        // The solution is the least-norm solution of || A x + b ||^2
        // First, compute c := - R[1:m, 1:m]^-T b with forward substitution
        c[0] = -b[0] / R[0][0];
        for (int i = 1; i < m; ++i)
        {
            c[i] = -b[i];
            for (int j = 0; j < i; ++j)
            {
                c[i] -= c[j] * R[j][i];
            }
            c[i] /= R[i][i];
        }

        // Then compute xN := Q[:, 1:m] c
        for (int i = 0; i < n; ++i)
        {
            double xi = 0;
            for (int j = 0; j < m; ++j)
            {
                xi += Q[i][j] * c[j];
            }
            xN.set(i, 0, xi);
        }
    }

    delete[] c;
    for (int i = 0; i < std::max(m, n); ++i)
    {
        delete [] Q[i];
        delete [] R[i];
        delete [] M[i];
    }
    delete [] Q;
    delete [] R;
    delete [] M;

    // Compute x according to the dogleg method.
    if (xN.norm() <= delta)
    {
        // Newton direction inside the trust-region: x := xN
        x = xN;
    }
    else if (xCp.norm() > delta)
    {
        // Both Newton direction and Cauchy point directions are outside the trust-region.
        // x := (delta / || xCp ||) xCp
        x = xCp;
        x.multiply(delta / xCp.norm());
    }
    else
    {
        // It is at the point of intersection of the dogleg and the trust-region
        // boundary.
        // x is the solution of the scalar quadratic equation:
        // || xCp + (tau - 1) (xN - xCp) ||^2 = delta^2.
        SGTELIB::Matrix xNmxCp = xN - xCp;
        const double q2 = xNmxCp.normsquare();
        const double q1 = 2 * SGTELIB::Matrix::dot(xCp, xNmxCp);
        const double q0 = xCp.normsquare() - delta * delta;

        double r1, r2;
        const bool has_roots = NOMAD::roots_quadratic(q2, q1, q0, r1, r2);
        if (!has_roots)
        {
            return NOMAD::DoglegTRSolverStatus::TR_NUM_ERROR;
        }
        // r := tau -1, with 1 <= tau <= 2.
        const double r = (r1 >= 0) && (r1 <= 1) ? r1 : r2;
        if ((r2 < 0) || (r2 > 1))
        {
            return NOMAD::DoglegTRSolverStatus::TR_NUM_ERROR;
        }

        x = xNmxCp;
        x.multiply(r);
        x = x + xCp;
    }

    return NOMAD::DoglegTRSolverStatus::SOLVED;
}
