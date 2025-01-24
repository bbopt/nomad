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
#ifndef __NOMAD_4_5_MATRIXUTILS__
#define __NOMAD_4_5_MATRIXUTILS__

#include <string>
#include "../Util/defines.hpp"

// Singular Value Decomposition (SVD) constants:
const double SVD_EPS      = 1e-13;      ///< Epsilon for SVD
const int    SVD_MAX_MPN  = 1500;       ///< Matrix maximal size (\c m+n )
const double SVD_MAX_COND = NOMAD::INF; ///< Max. acceptable cond. number

#include "../nomad_platform.hpp"
#include "../nomad_nsbegin.hpp"

/// Utilities for matrix.
/**
 Utilities for Matrix. Not everything that would be expected for a Matrix class is here. The methods are added when they are necessary for the code.
*/

/// SVD decomposition.
/**
 - The \c mxn \c M matrix is decomposed into \c M=U.W.V'.
 \param error_msg Error message when the function returns \c false    -- \b OUT.
 \param M         The input \c mxn matrix; Will be replaced by \c U   -- \b IN/OUT.
 \param W         The output \c nxn diagonal matrix                   -- \b OUT.
 \param V         The output \c nxn matrix                            -- \b OUT.
 \param m         Number of rows in M                                 -- \b IN.
 \param n         Number of columns in M                              -- \b IN.
 \param max_mpn   Maximum allowed value for \c m+n; ignored if \c <=0 -- \b IN
 (Opt) (default = \c 1500).
 \return          \c true if the decomposition worked.
 */
DLL_UTIL_API bool SVD_decomposition(std::string& error_msg,
                       double ** M,
                       double  * W,
                       double ** V,
                       int       m,
                       int       n,
                       int       max_mpn = 1500);


/// LU decomposition.
/**
 - The \c nxn \c M matrix is decomposed into \c M=L.U.
 \param error_msg Error message when the function returns \c false    -- \b OUT.
 \param M         The input \c nxn matrix; Will be replaced by \c LU   -- \b IN/OUT.
 \param n         Number of columns and rows in M                              -- \b IN.
 \param d         Used by determinant                          -- \b OUT.
 \param max_mpn   Maximum allowed value for \c n; ignored if \c <=0 -- \b IN.
 \return          \c true if the decomposition worked.
 */
DLL_UTIL_API bool LU_decomposition(std::string & error_msg,
                       double ** M,
                       int       n,
                       double&   d,
                       int       max_mpn = 1500);

/// QR decomposition.
/**
 * QR factorization for a rectangular m x n matrix.

 - The \c mxn \c M matrix is decomposed into \c M = Q x R
 \param error_msg Error message when the function returns \c false    -- \b OUT.
 \param M         The input \c nxn matrix;  -- \b IN.
 \param Q         The input \c mxm matrix;  -- \b OUT.
 \param R         The input \c mxn matrix;  -- \b OUT.
 \param m         Number of rows in M -- \b IN.
 \param n         Number of columns in M -- \b IN.
 \param max_mpn   Maximum allowed value for \c m or \c n; ignored if \c <=0 -- \b IN.
 \return          \c true if the decomposition worked.

 The QR factorization is done using a variant of the method of Givens rotations, proposed
 by Alexis Montoison.

 Also it is slower than the QR factorization using Householder reflections for dense matrices
 (this is the default implementation in LAPACK library), it is also more stable numerically,
 and easier to parallelize than the Householder approach.

 This implementation could be improved, especially in the storing of the reflection matrices
 coefficients. For the small matrix sizes involved in this context, the memory cost is not
 a problem.
 */
DLL_UTIL_API bool qr_factorization(std::string & error_msg,
                                   double ** M,
                                   double ** Q,
                                   double ** R,
                                   const int m,
                                   const int n,
                                   const int max_mpn = 1500);

/// LDLt decomposition.
/**
 * Block LDL' factorization for a symmetric indefinite matrix

 This code is adapted from LAPACK linear library and bunchkaufmann native implementation from
 the Julia project.

 Given a symmetric matrix M of order n, this procedure computes the decomposition
 P' M P = L D L'
 where L is a unit lower triangular matrix. D is a block diagonal matrix with blocks
 of order 1 or 2 and L[i+1,i] = 0 when D[i+1,i] is non-zero, and P a permutation matrix.

 This procedure uses the `partial pivoting` strategy of Bunch-Kaufman to form the decomposition.
 M is assumed to be stored only in its lower triangular part. L and D are written over M, and the
 diagonal of M will be destroyed. On output, the real array change of length n contains a record
 of the interchanges performed, i.e. the permutation matrix P. This permutation matrix contains the
 ``block indices'', following the convention of LAPACK:
         | k + 1 if ind corresponds to block of size 1x1,
 P[ind] =| else
         | -(k + 1) if ind corresponds to block of size 2x2
 where k is the corresponding permutation index.

 It is then possible to explicitly recompute L and D from the overwritten matrix M.

 \param error_msg Error message when the function returns \c false -- \b OUT.
 \param M         The input-output \c nxn matrix; -- \b IN OUT.
 \param pivots    Permutation vector of size \c n ; -- \b OUT.
 \param n         Number of columns and rows in M   -- \b IN.
 \param max_mpn   Maximum allowed value for \c n; ignored if \c <= 0 -- \b IN.
 \return          \c true if the decomposition worked.
 */
DLL_UTIL_API bool LDLt_factorization(std::string& error_msg,
                                     double** M,
                                     int* pivots,
                                     const int n,
                                     const int max_mpn = 1500);

/**
 * Extract the components L and D and permutation vector p from the overwritten matrix M,
   computed by the LDLt_factorization procedure.

 \param error_msg      Error message when the function returns \c false -- \b OUT.
 \param M              The input-output \c nxn matrix; -- \b IN.
 \param perm_lapack    Permutation vector of size \c n (in LAPACK format); -- \b IN.
 \param L              The lower triangular matrix component \c nxn matrix; -- \b OUT.
 \param D              The block diagonal (1x1 or 2x2) matrix component \c nxn matrix; -- \b OUT.
 \param p              Permutation vector of size \c n; -- \b OUT.
 \param n              Number of columns and rows in M   -- \b IN.
 \param max_mpn        Maximum allowed value for \c n; ignored if \c <= 0 -- \b IN.
 \return               \c true if the extraction worked worked.
 */
DLL_UTIL_API bool LDLt_extract(std::string& error_msg,
                               double** M,
                               const int* perm_lapack,
                               double** L,
                               double** D,
                               int* p,
                               const int n,
                               const int max_mpn = 1500);

/**
 * Block LDL^T factorization for a symmetric indefinite matrix.

 - The \c nxn \c M matrix is decomposed into \c M=L.D.L^t
 \param error_msg Error message when the function returns \c false    -- \b OUT.
 \param M         The input \c nxn matrix;  -- \b IN.
 \param L         The input \c nxn matrix; Triangle matrix with unit diagonal  -- \b OUT.
 \param D         The input \c nxn matrix; Block diagonal matrix with block of size \c 1x1 or \c 2x2 -- \b OUT.
 \param pp        Permutation vector of size \c n ; -- \b OUT.
 \param n         Number of columns and rows in M                              -- \b IN.
 \param max_mpn   Maximum allowed value for \c n; ignored if \c <=0 -- \b IN.
 \return          \c true if the decomposition worked.

 We use the original `partial pivoting` strategy from Bunch and Kaufman.
 */
DLL_UTIL_API bool LDLt_decomposition(std::string & error_msg,
                       double ** M,
                       double ** L,
                       double ** D,
                       int    * pp,
                       int       n,
                       int       max_mpn = 1500);

/// Block-diagonal regularization.
/**
 * Block LDL^T factorization for a symmetric indefinite matrx.

 \param D         The input \c nxn matrix; Block diagonal matrix with block of size \c 1x1 or \c 2x2 -- \b OUT.
 \param n         Number of columns and rows in M                              -- \b IN.
 \return          \c true if the decomposition worked.

Modify the block-diagonal matrix such that it becomes positive definite.
 */
DLL_UTIL_API bool DiagRegularization(double ** D, int n);

/// Block-diagonal regularization.
/**
 * Block LDL^T factorization for a symmetric indefinite matrx.

 \param D         The input \c nxn matrix; Block diagonal matrix with block of size \c 1x1 or \c 2x2 -- \b OUT.
 \param n         Number of columns and rows in M                              -- \b IN.
 \return          \c true if the decomposition worked.

Return the smallest eigenvalue of D.
 */
DLL_UTIL_API double FindSmallestEigenvalue(double ** D, int n);

/// LDL' solve
/**
 Solve linear system
 M x = rhs
 where M has been factored by the LDLt_factorization procedure.
 Note that only the lower triangular part of M is defined (i.e., M[i][j] with j <= i).

 \param error_msg Error message when the function returns \c false    -- \b OUT.
 \param M         The input factored \c nxn matrix; -- \b IN.
 \param rhs       Right-hand side vector of size \c n of the linear system; -- \b IN.
 \param sol       Solution vector of size \c n; -- \b OUT.
 \param pivots    Permutation vector of size \c n; -- \b IN.
 \param n         Number of columns and rows in M; -- \b IN.
 \return          \c true if the resolution has worked.
 */
DLL_UTIL_API bool ldl_solve(std::string& error_msg,
                            double** M,
                            const double* rhs,
                            double* sol,
                            const int* pivots,
                            const int n);

/// LDLt solve.
/**
 * Solve linear system given a block LDL^T factorization for a symmetric indefinite matrix.

 - The \c nxn \c M matrix is decomposed into \c M=L.D.L^t
 \param error_msg Error message when the function returns \c false    -- \b OUT.
 \param L         The input \c nxn matrix; Triangle matrix with unit diagonal  -- \b IN.
 \param D         The input \c nxn matrix; Block diagonal matrix with block of size \c 1x1 or \c 2x2 -- \b IN.
 \param rhs       Right-hand side vector of size \c n of the linear system ; -- \b IN.
 \param sol       Solution vector of size \c n ; -- \b OUT.
 \param pp        Permutation vector of size \c n ; -- \b IN.
 \param n         Number of columns and rows in M                              -- \b IN.
 \return          \c true if the solve worked.

This function successively uses ldl_dsolve, ldl_ltsolve and ldl_lsolve.
 */
DLL_UTIL_API bool ldl_solve(std::string & error_msg,
    double       ** D,
    double       ** L,
    const double * rhs,
    double       * sol,
    const int    * pp,
    int          n);

DLL_UTIL_API bool ldl_dsolve( double ** D, const double * rhs, double * Ly, int n);

DLL_UTIL_API bool ldl_ltsolve( double ** L, const double * rhs, double * Ly, int n);

DLL_UTIL_API bool ldl_lsolve( double ** L, const double * rhs, double * Ly, int n);

// Get rank of a matrix  using SVD decomposition
/**
 - The \c mxn \c M matrix is decomposed into \c M=U.W.V'. The rank equals the size of W
 \param M         The input \c mxn matrix                               -- \b IN.
 \param m         Number of rows in M                                   -- \b IN.
 \param n         Number of columns in M                                -- \b IN.
 \param eps       Precision to detect rank from W                       -- \b IN.
 \return          The rank>0 if the decomposition worked else 0.
 */
DLL_UTIL_API int getRank(double ** M,
            size_t    m,
            size_t    n,
            double    eps = SVD_EPS);


// Get determinant of a matrix
/**
 \param M         The input \c nxn matrix         -- \b IN.
 \param n         Number of row and columns in M  -- \b IN.
 \param det       The determinant of the matrix M -- \b OUT.
 \return True/False if success or not.
 */
DLL_UTIL_API bool getDeterminant(double **M,
                    double & det,
                    size_t n);



#include "../nomad_nsend.hpp"
#endif // __NOMAD_4_5_MATRIXUTILS__
