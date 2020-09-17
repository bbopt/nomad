#ifndef __NOMAD400_MATRIXUTILS__
#define __NOMAD400_MATRIXUTILS__

#include <string>
#include "../Util/defines.hpp"

// Singular Value Decomposition (SVD) constants:
const double SVD_EPS      = 1e-13;      ///< Epsilon for SVD
const int    SVD_MAX_MPN  = 1500;       ///< Matrix maximal size (\c m+n )
const double SVD_MAX_COND = NOMAD::INF; ///< Max. acceptable cond. number

#include "../nomad_nsbegin.hpp"

/// Utilities for matrix.
/**
 Utilities for Matrix. Not everything that would be expected for a Matrix class is here. The methods are added when they are necessary for the code.
 \todo Implement a Matrix class.
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
bool SVD_decomposition(std::string& error_msg,
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
bool LU_decomposition(std::string & error_msg,
                       double ** M,
                       int       n,
                       double&   d,
                       int       max_mpn = 1500);

// Get rank of a matrix  using SVD decomposition
/**
 - The \c mxn \c M matrix is decomposed into \c M=U.W.V'. The rank equals the size of W
 \param M         The input \c mxn matrix                               -- \b IN.
 \param m         Number of rows in M                                   -- \b IN.
 \param n         Number of columns in M                                -- \b IN.
 \param eps       Precision to detect rank from W                       -- \b IN.
 \return          The rank>0 if the decomposition worked else 0.
 */
int getRank(double ** M,
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
bool getDeterminant(double **M,
                    double & det,
                    size_t n);



#include "../nomad_nsend.hpp"
#endif // __NOMAD400_MATRIXUTILS__
