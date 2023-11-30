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
#include "../Math/MatrixUtils.hpp"
#include "../Util/utils.hpp"

#include <math.h>
#include <algorithm>

bool NOMAD::getDeterminant(double ** M,
                            double  & det,
                            size_t n)
{
    std::string error_msg;

    double d = 1 ;

    double ** LU = new double *[n];
    for (size_t i = 0 ; i < n ; i++ )
    {
        LU[i]=new double [n];
        for ( size_t j = 0; j < n ; j++ )
            LU[i][j]=M[i][j];
    }

    NOMAD::LU_decomposition ( error_msg, LU,  static_cast<int>(n), d );


    if ( error_msg.empty() )
    {
        for (size_t i=0;i < n;i++)
            d*=LU[i][i];
    }


    for (size_t i = 0 ; i < n ; i++ )
    {
        delete [] LU[i];
    }
    delete [] LU;

    det = d;

    if ( error_msg.empty() )
        return true;
    else
        return false;
}


int NOMAD::getRank(double ** M,
                    size_t m,
                    size_t n,
                    double    eps)
{
    double  * W = new double  [n];
    double ** V = new double *[n];
    for (size_t i = 0 ; i < n ; ++i )
    {
        V[i]=new double [n];
    }

    std::string error_msg;
    NOMAD::SVD_decomposition ( error_msg, M, W, V, static_cast<int>(m), static_cast<int>(n) );

    for (size_t i=0;i<n;++i)
        delete [] V[i];
    delete [] V;


    if (! error_msg.empty())
    {
        delete [] W;
        return -1;
    }

    int rank=0;
    for (size_t i=0;i<n;i++)
    {
        if (fabs(W[i]) > eps)
            rank++;
    }

    delete [] W;
    return rank;

}


/*--------------------------------------------------------------*/
/*                        SVD decomposition                     */
/*  inspired and recoded from an old numerical recipes version  */
/*--------------------------------------------------------------*/
/*                                                              */
/*           M = U . W . V'                                     */
/*                                                              */
/*           M ( m x n )                                        */
/*           U ( m x n )                                        */
/*           W ( n x n )                                        */
/*           V ( n x n )                                        */
/*                                                              */
/*           U.U' = U'.U = I if m = n                           */
/*           U'.U = I        if m > n                           */
/*                                                              */
/*           V.V' = V'.V = I                                    */
/*                                                              */
/*           W diagonal                                         */
/*                                                              */
/*           M is given as first argument and becomes U         */
/*           W is given as a size-n vector                      */
/*           V is given, not V'                                 */
/*                                                              */
/*--------------------------------------------------------------*/
bool NOMAD::SVD_decomposition ( std::string & error_msg,
                               double     ** M,
                               double      * W,
                               double     ** V,
                               int           m,
                               int           n,
                               int           max_mpn     ) // default = 1500
{
    error_msg.clear();

    if ( max_mpn > 0 && m+n > max_mpn )
    {
        error_msg = "SVD_decomposition() error: m+n > " + NOMAD::itos ( max_mpn );
        return false;
    }

    double * rv1   = new double[n];
    double   scale = 0.0;
    double   g     = 0.0;
    double   norm  = 0.0;

    int      nm1   = n - 1;

    bool   flag;
    int    i, j, k, l = 0, its, jj, nm = 0;
    double s, f, h, c, x, y, z, absf, absg, absh;

    const int NITER = 30;

    // Initialization W and V
    for (i=0; i < n; ++i)
    {
        W[i]=0.0;
        for (j=0; j < n ; ++j)
            V[i][j]=0.0;
    }

    // Householder reduction to bidiagonal form:
    for ( i = 0 ; i < n ; ++i )
    {
        l      = i + 1;
        rv1[i] = scale * g;
        g      = s = scale = 0.0;
        if ( i < m )
        {
            for ( k = i ; k < m ; ++k )
                scale += fabs ( M[k][i] );
            if ( scale )
            {
                for ( k = i ; k < m ; ++k )
                {
                    M[k][i] /= scale;
                    s += M[k][i] * M[k][i];
                }
                f       = M[i][i];
                g       = ( f >= 0.0 ) ? -fabs(sqrt(s)) : fabs(sqrt(s));
                h       = f * g - s;
                M[i][i] = f - g;
                for ( j = l ; j < n ; ++j )
                {
                    for ( s = 0.0, k = i ; k < m ; ++k )
                        s += M[k][i] * M[k][j];
                    f = s / h;
                    for ( k = i ; k < m ; ++k )
                        M[k][j] += f * M[k][i];
                }
                for ( k = i ; k < m ; ++k )
                    M[k][i] *= scale;
            }
        }
        W[i] = scale * g;
        g    = s = scale = 0.0;
        if ( i < m && i != nm1 )
        {
            for ( k = l ; k < n ; ++k )
                scale += fabs ( M[i][k] );
            if ( scale )
            {
                for ( k = l ; k < n ; ++k )
                {
                    M[i][k] /= scale;
                    s       += M[i][k] * M[i][k];
                }
                f       = M[i][l];
                g       = ( f >= 0.0 ) ? -fabs(sqrt(s)) : fabs(sqrt(s));
                h       = f * g - s;
                M[i][l] = f - g;
                for ( k = l ; k < n ; ++k )
                    rv1[k] = M[i][k] / h;
                for ( j = l ; j < m ; ++j )
                {
                    for ( s=0.0, k = l ; k < n ; ++k )
                        s += M[j][k] * M[i][k];
                    for ( k = l ; k < n ; ++k )
                        M[j][k] += s * rv1[k];
                }
                for ( k = l ; k < n ; ++k )
                    M[i][k] *= scale;
            }
        }
        double tmp  = fabs ( W[i] ) + fabs ( rv1[i] );
        norm = ( norm > tmp ) ? norm : tmp;
    }

    // accumulation of right-hand transformations:
    for ( i = nm1 ; i >= 0 ; --i )
    {
        if ( i < nm1 )
        {
            if ( g )
            {
                for ( j = l ; j < n ; ++j )
                    V[j][i] = ( M[i][j] / M[i][l] ) / g;
                for ( j = l ; j < n ; ++j )
                {
                    for ( s = 0.0, k = l ; k < n ; ++k )
                        s += M[i][k] * V[k][j];
                    for ( k = l ; k < n ; ++k )
                        V[k][j] += s * V[k][i];
                }
            }
            for ( j = l ; j < n ; ++j )
                V[i][j] = V[j][i] = 0.0;
        }
        V[i][i] = 1.0;
        g       = rv1[i];
        l       = i;
    }

    // accumulation of left-hand transformations:
    for ( i = ( ( m < n ) ? m : n ) - 1 ; i >= 0 ; --i )
    {
        l = i + 1;
        g = W[i];
        for ( j = l ; j < n ; ++j )
            M[i][j] = 0.0;
        if ( g )
        {
            g = 1.0 / g;
            for ( j = l ; j < n ; ++j )
            {
                for ( s = 0.0, k = l ; k < m ; ++k )
                    s += M[k][i] * M[k][j];
                f = ( s / M[i][i] ) * g;
                for ( k = i ; k < m ; ++k )
                    M[k][j] += f * M[k][i];
            }
            for ( j = i ; j < m ; ++j )
                M[j][i] *= g;
        }
        else
            for ( j = i ; j < m ; ++j )
                M[j][i] = 0.0;
        ++M[i][i];
    }

    // diagonalization of the bidiagonal form:
    for ( k = nm1 ; k >= 0 ; --k )
    {
        for ( its = 1 ; its <= NITER ; its++ )
        {
            flag = true;
            for ( l = k ; l >= 0 ; l-- )
            {
                nm = l - 1;
                if ( nm < 0 || fabs ( rv1[l]) + norm == norm )
                {
                    flag = false;
                    break;
                }
                if ( fabs ( W[nm] ) + norm == norm )
                    break;
            }
            if ( flag )
            {
                c = 0.0;
                s = 1.0;
                for ( i = l ; i <= k ; i++ )
                {
                    f      = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ( fabs(f) + norm == norm )
                        break;
                    g = W[i];

                    absf = fabs(f);
                    absg = fabs(g);
                    h    = ( absf > absg ) ?
                    absf * sqrt ( 1.0 + pow ( absg/absf, 2.0 ) ) :
                    ( ( absg==0 ) ? 0.0 : absg * sqrt ( 1.0 + pow ( absf/absg, 2.0 ) ) );

                    W[i] =  h;
                    h    =  1.0 / h;
                    c    =  g * h;
                    s    = -f * h;
                    for ( j = 0 ; j < m ; ++j )
                    {
                        y = M[j][nm];
                        z = M[j][ i];
                        M[j][nm] = y * c + z * s;
                        M[j][ i] = z * c - y * s;
                    }
                }
            }
            z = W[k];
            if ( l == k)
            {
                if ( z < 0.0 )
                {
                    W[k] = -z;
                    for ( j = 0 ; j < n ; j++ )
                        V[j][k] = -V[j][k];
                }
                break;  // this 'break' is always active if k==0
            }
            if ( its == NITER )
            {
                error_msg = "SVD_decomposition() error: no convergence in " +
                NOMAD::itos ( NITER ) + " iterations";
                delete [] rv1;
                return false;
            }
            x  = W[l];
            nm = k - 1;
            y  = W[nm];
            g  = rv1[nm];
            h  = rv1[k];
            f  = ( (y-z) * (y+z) + (g-h) * (g+h) ) / ( 2.0 * h * y );

            absf = fabs(f);
            g    = ( absf > 1.0 ) ?
            absf * sqrt ( 1.0 + pow ( 1.0/absf, 2.0 ) ) :
            sqrt ( 1.0 + pow ( absf, 2.0 ) );

            f = ( (x-z) * (x+z) +
                 h * ( ( y / ( f + ( (f >= 0)? fabs(g) : -fabs(g) ) ) ) - h ) ) / x;
            c = s = 1.0;

            for ( j = l ; j <= nm ; ++j )
            {
                i = j + 1;
                g = rv1[i];
                y = W[i];
                h = s * g;
                g = c * g;

                absf = fabs(f);
                absh = fabs(h);
                z    = ( absf > absh ) ?
                absf * sqrt ( 1.0 + pow ( absh/absf, 2.0 ) ) :
                ( ( absh==0 ) ? 0.0 : absh * sqrt ( 1.0 + pow ( absf/absh, 2.0 ) ) );

                rv1[j] = z;
                c      = f / z;
                s      = h / z;
                f      = x * c + g * s;
                g      = g * c - x * s;
                h      = y * s;
                y     *= c;
                for ( jj = 0 ; jj < n ; ++jj )
                {
                    x = V[jj][j];
                    z = V[jj][i];
                    V[jj][j] = x * c + z * s;
                    V[jj][i] = z * c - x * s;
                }

                absf = fabs(f);
                absh = fabs(h);
                z    = ( absf > absh ) ?
                absf * sqrt ( 1.0 + pow ( absh/absf, 2.0 ) ) :
                ( ( absh==0 ) ? 0.0 : absh * sqrt ( 1.0 + pow ( absf/absh, 2.0 ) ) );

                W[j] = z;

                if ( z )
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for ( jj = 0 ; jj < m ; ++jj )
                {
                    y = M[jj][j];
                    z = M[jj][i];
                    M[jj][j] = y * c + z * s;
                    M[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            W  [k] = x;
        }
    }

    delete [] rv1;
    return true;
}


/*--------------------------------------------------------------*/
/*                        LU decomposition                      */
/*  Recoded from numerical recipes 3rd ed. 2007                 */
/*--------------------------------------------------------------*/
/*                                                              */
/*           M = L.U                                            */
/*                                                              */
/*           M ( n x n )                                        */
/*                                                              */
/*           M is given as first argument and becomes LU        */
/*           LU = [[b00 b01 ... b0n]                            */
/*                 [a10 b11 ... b1n]                            */
/*                 [a20 a21 ... b2n]                            */
/*                    ...........                               */
/*                 [an0 an1 ... bnn]                            */
/*                                                              */
/*--------------------------------------------------------------*/
bool NOMAD::LU_decomposition ( std::string & error_msg,
                              double     ** LU,
                              int           n,
                              double      & d,
                              int           max_n     ) // default=50
{
    error_msg.clear();

    if ( max_n > 0 && n > max_n )
    {
        error_msg = "LU_decomposition() error: n > " + NOMAD::itos ( max_n );
        return false;
    }

    const double TINY = 1E-40;
    int i, j, k;

    double big, temp;

    double * vv = new double[n]; // stores the implicit scaling of each row
    int *indx = new int[n]; // stores the permutations

    d = 1; // No row interchange yet

    for (i = 0; i < n ; i++ )  // Loop over row to get implicit scaling information
    {
        big = 0.0;
        for ( j = 0 ; j < n ; j++ )
        {
            if ( (temp = fabs(LU[i][j])) > big )
                big = temp;
        }
        if ( big == 0 )
        {
            error_msg = "LU_decomposition() error: no nonzero largest element";
            delete [] vv;
            delete [] indx;

            return false;
        }
        vv[i]= 1.0/big; // Saves the scaling
    }
    for ( k = 0; k < n ; k++) // This is the outermost kij loop
    {
        big = 0.0;            // Initialized the search for largest pivot element
        int imax = k;
        for ( i = k ; i < n ; i++ )
        {
            temp = vv[i]*fabs(LU[i][k]);
            if ( temp > big )
            {
                big = temp;
                imax = i;
            }
        }
        if ( k != imax )
        {
            for ( j = 0; j < n ; j++ )
            {
                temp = LU[imax][j];
                LU[imax][j]=LU[k][j];
                LU[k][j]=temp;
            }
            d = -d;
            vv[imax] = vv[k];
        }
        indx[k] = imax;
        if ( LU[k][k] == 0.0) // TINY for zero
            LU[k][k] = TINY;
        for ( i = k+1 ; i < n ; i++ )
        {
            temp = LU[i][k] /= LU[k][k]; // Divide by pivot element
            for ( j = k+1 ; j < n ; j++ ) // Innermos matrix: reduce remaining submatrix
                LU[i][j] -= temp*LU[k][j];
        }
    }

    delete [] vv;
    delete [] indx;
    return true;

}

/*--------------------------------------------------------------*/
/*                        QR factorization                      */
/*                                                              */
/*--------------------------------------------------------------*/
/*                                                              */
/*           M = Q.R                                            */
/*                                                              */
/*           M ( m x n )                                        */
/*           Q ( m x m ): ....                                  */
/*           R ( m x n ): lower triangular                      */
/*                                                              */
/*--------------------------------------------------------------*/
bool NOMAD::qr_factorization ( std::string & error_msg,
                              double     ** M,
                              double     ** Q,
                              double     ** R,
                              int           m,
                              int           n,
                              int           max_n     ) // default = 1500
{
    error_msg.clear();

    if ( max_n > 0 && (n > max_n || m > max_n) )
    {
        error_msg = "qr_factorization() error: min(m,n) > " + NOMAD::itos ( max_n );
        return false;
    }

    int  i, j, k;
    bool chk = true;

    double ** v = new double *[m]; // size m x n
    for ( i = 0 ; i < m ; ++i )
        v[i] = new double [n];
    double * normv = new double [n];

    // Gram-Schmidt (GS): First step: v1=u1 (first column of A):
    normv[0] = 0.0;
    for ( i = 0 ; i < m ; ++i )
    {
        v[i][0] = M[i][0];
        normv[0] += pow(v[i][0], 2.0);
    }
    if ( normv[0] == 0.0 )
        chk = false;

    // GS: Compute vj=uj-projections:
    for ( j = 1 ; j < n ; ++j )
    {

        if ( !chk )
            break;

        // Init: vj=uj:
        for ( i = 0 ; i < m ; ++i )
            v[i][j] = M[i][j];

        // Compute projections:
        for ( k = 0 ; k < j ; ++k )
        {
            // pk: projection of uj to vk:
            double tmp = 0.0;
            for ( i = 0 ; i < m ; ++i )
                tmp += v[i][k] * M[i][j];
            for ( i = 0 ; i < m ; ++i )
                v[i][j] -= tmp * v[i][k] / normv[k];
        }

        // Compute norm of vj:
        normv[j] = 0.0;
        for ( i = 0 ; i < m ; ++i )
            normv[j] += pow(v[i][j], 2.0);
        if ( normv[j] == 0.0 )
            chk = false;
    }

    // Final step of GS: Normalize the v vectors into the Q matrix:
    if ( chk )
        for ( j = 0 ; j < n ; ++j )
            for ( i = 0 ; i < m ; ++i )
                Q[i][j] = v[i][j] / pow(normv[j],0.5);

    // Delete v and normv:
    for ( i = 0 ; i < m ; ++i )
        delete [] v[i];
    delete [] v;
    delete [] normv;

    // Final step of the QR decomposition: Compute R:
    for ( i = 0 ; i < n ; ++i )
    {
        for ( j = 0 ; j < i ; ++j )
            R[i][j] = 0.0;
        for ( j = i ; j < n ; ++j )
        {
            R[i][j] = 0.0;
            for ( k = 0 ; k < m ; ++k )
                R[i][j] += Q[k][i] * M[k][j];
        }
    }

    // Reset Q and R in case of failure:
    if ( !chk )
    {
        for ( j = 0 ; j < n ; ++j )
        {
            for ( i = 0 ; i < m ; ++i )
                Q[i][j] = 0.0;
            for ( i = 0 ; i < n ; ++i )
                R[i][j] = 0.0;
        }
        return false;
    }

    return true;
}

/*--------------------------------------------------------------*/
/*                        LDLt decomposition                    */
/*                                                              */
/*--------------------------------------------------------------*/
/*                                                              */
/*           M = L.D.L^T                                        */
/*                                                              */
/*           M ( n x n )                                        */
/*           D ( n x n ): block-diagonal of size 1 x 1 or 2 x 2 */
/*           L ( n x n ): lower triangular with unit diagonal   */
/*                                                              */
/*--------------------------------------------------------------*/
bool NOMAD::LDLt_decomposition ( std::string & error_msg,
                              double     ** M,
                              double     ** L,
                              double     ** D,
                              int        * pp,
                              int           n,
                              int           max_n     ) // default = 1500
{
    error_msg.clear();

    int piv = 1; // We use Bunch and Kaufman partial pivoting, but rook pivoting could also be implemented.

    if ( max_n > 0 && n > max_n )
    {
        error_msg = "LDLt_decomposition() error: n > " + NOMAD::itos ( max_n );
        return false;
    }

    // Initialization of L and D as identity matrices
    for (int i = 0 ; i < n ; i++ )
    {
        for ( int j = 0; j < n ; j++ )
            if (j == i)
            {
                D[i][j]=1;
            }
            else{
                D[i][j]=0;
            }
    }

    for (int i = 0 ; i < n ; i++ )
    {
        for ( int j = 0; j < n ; j++ )
            if (j == i)
            {
                L[i][j]=1;
            }
            else{
                L[i][j]=0;
            }
    }

    double * Mik = new double[n];
    double * Mki = new double[n];

    double detE;
    double tmp;

    for (int i = 0 ; i < n ; i++ )
    {
        pp[i] = i;
    }

    double rho = 1; // rho = norm(A, Inf)

    int s = 1; // size of the current block

    int m1=0;
    int m2=0;

    double alpha = (1 + sqrt(27)) / 8;
    double sigma;

    // Temporary variables:
    int vr=0;
    double lambda;
    int r=0;
    bool swap;

    int k = 0;
    while (k < n - 1)
    {
        // lambda, vr = findmax(abs.(A[k+1:n, k]))
        vr = k + 1;
        lambda = abs(M[vr][k]);
        for (int j = k + 2; j < n; j++)
        {
            if (abs(M[j][k]) > lambda)
            {
                vr = int(j);
                lambda = abs(M[vr][k]);
            }
        }

        r = vr;
        if (lambda > 0)
        {

            swap = false;
            if (abs(M[k][k]) >= alpha * lambda)
            {
                s = 1;
            }
            else
            {

                if (piv == 1) // piv = 'p'
                {
                    sigma = -1;
                    for (int i = k; i < n; i++)
                    {
                        if (abs(M[i][r]) > sigma)
                        {
                            sigma = abs(M[i][r]); // σ = norm(A[k:n, r], Inf)
                        }
                    }
                    if (alpha * lambda * lambda <= abs(M[k][k]) * sigma)
                    {
                        s = 1;
                    }
                    else if ( abs(M[r][r]) >= alpha * sigma )
                    {
                        swap = true;
                        m1 = k;
                        m2 = r;
                        s = 1;
                    }
                    else
                    {
                        swap = true;
                        m1 = k + 1;
                        m2 = r;
                        s = 2;
                    }
                    if (m1 < 0 || m1 >= n || m2 < 0 || m2 >= n)
                    {
                        return false;
                    }
                    if (swap)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            // A[[m1, m2], :] = A[[m2, m1], :]
                            tmp = M[m1][i];
                            M[m1][i] = M[m2][i];
                            M[m2][i] = tmp;
                            // L[[m1, m2], :] = L[[m2, m1], :]
                            tmp = L[m1][i];
                            L[m1][i] = L[m2][i];
                            L[m2][i] = tmp;
                        }
                        for (int i = 0; i < n; i++)
                        {
                            // A[:, [m1, m2]] = A[:, [m2, m1]]
                            tmp = M[i][m1];
                            M[i][m1] = M[i][m2];
                            M[i][m2] = tmp;
                            // L[:, [m1, m2]] = L[:, [m2, m1]]
                            tmp = L[i][m1];
                            L[i][m1] = L[i][m2];
                            L[i][m2] = tmp;
                        }
                        // pp[[m1, m2]] = pp[[m2, m1]]
                        tmp = pp[m1];
                        pp[m1] = pp[m2];
                        pp[m2] = static_cast<int>(tmp);
                    }

                }
                else if (piv == 2) // piv = 'p'
                {

                }
                else
                {
                    error_msg = "Not implemented pivoting piv, only 1 and 2.";
                    return false;
                }

            }

            if (s == 1)
            {
                D[k][k] = M[k][k];
                for (int j = k + 1; j < n; j++)
                {
                    M[j][k] = M[j][k] / M[k][k];
                    L[j][k] = M[j][k];
                }

                // A[k+1:n][k+1:n] = A[k+1:n][k+1:n] - A[k+1:n, k:k] * A[k:k][k+1:n] // rank-1 update
                for (int i = 0; i < n; i++)
                {
                    Mki[i] = M[i][k]; // A[k+1:n, k:k]
                    Mik[i] = M[k][i]; // A[k:k, k+1:n]
                }
                for (int i = k + 1; i < n; i++)
                {
                    for (int j = k + 1; j < n; j++)
                    {
                        M[i][j] = M[i][j] - Mki[i] * Mik[j];
                    }
                }

                // A[k+1:n][k+1:n] = 0.5 * (A[k+1:n, k+1:n] + A[k+1:n,k+1:n]')
               for (int i = k + 1; i < n; i++)
                {
                    for (int j = k + 1; j <= i; j++)
                    {
                        M[i][j] = 0.5 * (M[i][j] + M[j][i]);
                    }
                }
                for (int i = k + 1; i < n; i++)
                {
                    for (int j = i + 1; j < n; j++)
                    {
                        M[i][j] = M[j][i];
                    }
                }

            }
            else if (s == 2)
            {
                if (k >= n)
                {
                    return false;
                }
                // E = A[k:k+1, k:k+1]
                // D[k:k+1, k:k+1] = E
                D[k][k] = M[k][k];
                D[k][k + 1] = M[k][k + 1];
                D[k + 1][k] = M[k + 1][k];
                D[k + 1][k + 1] = M[k + 1][k + 1];

                // i = k+2:n
                // C = A[i, k:k+1]
                // temp = C / E // So, temp * E = C, so temp is of size (n x 2)
                detE = D[k][k] * D[k + 1][k + 1] - D[k][k + 1] * D[k + 1][k];
                if (detE == 0)
                {
                    error_msg = "Error in pivoting: determinant of block diagonal is 0.";
                    return false;
                }

                // L[i, k:k+1] = temp // So, for all i, L[i, k:k+1] = E^{-1} * C[i]
                for (int i = k + 2; i < n; i++)
                {
                    L[i][k] = (D[k + 1][k + 1] * M[i][k] - D[k][k + 1] * M[i][k + 1]) / detE;
                    L[i][k + 1] = (- D[k + 1][k] * M[i][k] + D[k][k] * M[i][k + 1]) / detE;
                }

                // A[i, k+2:n] = A[i, k+2:n] - temp * C'
                for (int i = k + 2; i < n; i++)
                {
                    for (int j = k + 2; j < n; j++)
                    {
                        M[i][j] = M[i][j] - L[i][k] * M[j][k]- L[i][k + 1] * M[j][k + 1];
                    }
                }

                // A[i, i] = 0.5 * (A[i, i] + A[i, i]')
                for (int i = k + 2; i < n; i++)
                {
                    for (int j = k + 2; j <= i; j++)
                    {
                        M[i][j] = 0.5 * (M[i][j] + M[j][i]);
                    }
                }
                for (int i = k + 2; i < n; i++)
                {
                    for (int j = i + 1; j < n; j++)
                    {
                        M[i][j] = M[j][i];
                    }
                }
            }
            else
            {
                // Unexpected value for s
                error_msg = "Unexpected value for s = " + std::to_string(s);
                return false;
            }
            if (k + s <= n - 1)
            {
                //val,_ = findmax(abs.(A[k+s:n, k+s:n]))
               double val = -1;
               for (int i = k + s; i < n; i++)
               {
                    for (int j = k + s; j < n; j++)
                    {
                        if (abs(M[i][j]) > val)
                        {
                            val = abs(M[i][j]);
                        }
                    }
               }
               rho = std::max(rho, val);
            }
        }
        else // nothing to do
        {
            s = 1;
            D[k][k] = M[k][k];
        }

        k = k + s;

        if (k == n - 1)
        {
            D[n - 1][n - 1] = M[n - 1][n - 1];
        }
    }

    delete [] Mik;
    delete [] Mki;

    return true;
}

bool NOMAD::DiagRegularization(double ** D, int n)
{
    int s; // it is either 1 or 2
    double detE; // determinant of the 2x2 block matrix
    double traceE; // trace of the 2x2 block matrix
    double * l = new double [2]; // eigenvalues of the 2x2 block matrix

    double diag;
    double findmin = 0; // min(0, smallest_eigenvalue) <= 0
    double findmax = 1; // max(1, largest_eigenvalue) >= 0
    double findsmall = 1; // min(min(|eigenvalues|), 1) >= 0

    int j;
    int k = 0;
    while (k < n)
    {
        if (n == 1 || k == n - 1 || (D[k][k + 1] == 0 && D[k + 1][k] == 0))
        {
            s = 1;
            l[0] = D[k][k];
        }
        else
        {
            s = 2;
            detE = D[k][k] * D[k + 1][k + 1] - D[k][k + 1] * D[k + 1][k];
            traceE = D[k][k] + D[k + 1][k + 1];
            // det = l1 * l2 and traceE = l1 + l2
            l[0] = 0.5 * (traceE + sqrt(pow(traceE, 2.0) - 4 * detE));
            l[1] = 0.5 * (traceE - sqrt(pow(traceE, 2.0) - 4 * detE));
        }
        for (j = 0 ; j < s; ++j)
        {
            diag = l[j];
            if (diag <= findmin)
            {
                findmin = diag;
            }
            if (diag >= findmax)
            {
                findmax = diag;
            }
            if (fabs(diag) <= findsmall)
            {
                findsmall = fabs(diag);
            }
        }
        k += s;
    }

    double regularizer = - findmin + (findsmall + findmax) / 2;

    if (findmin < 0)
    {
        int k = 0;
        while (k < n)
        {
            if (n == 1 || k == n - 1 || (D[k][k + 1] == 0 && D[k + 1][k] == 0))
            {
                s = 1;
                diag = D[k][k];
                if (diag < 0)
                {
                    D[k][k] += regularizer;
                }
            }
            else
            {
                s = 2;
                detE = D[k][k] * D[k + 1][k + 1] - D[k][k + 1] * D[k + 1][k];
                traceE = D[k][k] + D[k + 1][k + 1];
                if ((detE <= 0) || (traceE < 0))
                {
                    D[k][k] += regularizer;
                    D[k + 1][k + 1] += regularizer;
                }
            }
            k += s;
        }
    }

    delete [] l;

    return true;
}

double NOMAD::FindSmallestEigenvalue(double ** D, int n)
{
    int s; // it is either 1 or 2
    double detE; // determinant of the 2x2 block matrix
    double traceE; // trace of the 2x2 block matrix
    double * l = new double [2]; // eigenvalues of the 2x2 block matrix

    double diag;
    double findmin = std::numeric_limits<double>::infinity();

    int j;
    int k = 0;
    while (k < n)
    {
        if (n == 1 || k == n - 1 || (D[k][k + 1] == 0 && D[k + 1][k] == 0))
        {
            s = 1;
            l[0] = D[k][k];
        }
        else
        {
            s = 2;
            detE = D[k][k] * D[k + 1][k + 1] - D[k][k + 1] * D[k + 1][k];
            traceE = D[k][k] + D[k + 1][k + 1];
            // det = l1 * l2 and traceE = l1 + l2
            l[0] = 0.5 * (traceE + sqrt(pow(traceE, 2.0) - 4 * detE));
            l[1] = 0.5 * (traceE - sqrt(pow(traceE, 2.0) - 4 * detE));
        }
        for (j = 0 ; j < s; ++j)
        {
            diag = l[j];
            if (diag <= findmin)
            {
                findmin = diag;
            }
        }
        k += s;
    }

    delete [] l;

    return findmin;
}

bool NOMAD::ldl_lsolve( double ** L, double * rhs, double * Ly, int n)
{
    for(int k = 0; k < n; k++)
    {
        Ly[k] = rhs[k];
        for(int j = 0; j < k; j++)
        {
            Ly[k] -= L[k][j] * Ly[j];
        }
    }
    return true;
}

bool NOMAD::ldl_ltsolve( double ** L, double * rhs, double * Ly, int n)
{
    for(int k = n - 1; k > -1; k--)
    {
        Ly[k] = rhs[k];
        for(int j = k + 1; j < n; j++)
        {
            Ly[k] -= L[j][k] * Ly[j];
        }
    }
    return true;
}

bool NOMAD::ldl_dsolve( double ** D,  double * rhs, double * Ly, int n)
{
    int s; // it is either 1 or 2
    double detE; // determinant of the 2x2 block matrix
    int k = 0;
    while (k < n)
    {
        if (n == 1 || k == n - 1 || (D[k][k + 1] == 0 && D[k + 1][k] == 0))
        {
            s = 1;
            if (D[k][k] == 0)
            {
                return false;
            }
            else
            {
                Ly[k] = rhs[k] / D[k][k];
            }
        }
        else
        {
            s = 2;
            detE = D[k][k] * D[k + 1][k + 1] - D[k][k + 1] * D[k + 1][k];
            if (detE == 0)
            {
                return false;
            }
            else
            {
                Ly[k] = (D[k + 1][k + 1] * rhs[k] - D[k][k + 1] * rhs[k + 1]) / detE;
                Ly[k + 1] = (- D[k + 1][k] * rhs[k] + D[k][k] * rhs[k + 1]) / detE;
            }
        }
        k += s;
    }
    return true;
}

bool NOMAD::ldl_solve( std::string & error_msg,
    double     ** D,
    double     ** L,
    double     * rhs,
    double     * sol,
    int        * pp,
    int            n
)
{
    error_msg.clear();
    bool success;

    double * prhs = new double [n];
    double * Lz = new double [n];
    for (int i = 0; i < n; i++)
    {
        prhs[i] = rhs[pp[i]];
        Lz[i] = 0.0;
    }
    success = ldl_lsolve(L, prhs, Lz, n); // Lz = rhs
    if (!success)
    {
        return false;
    }
    double * Dy = new double [n];
    for (int i = 0; i < n; i++)
    {
        Dy[i] = 0.0;
    }
    success = ldl_dsolve(D, Lz, Dy, n); // Dy = z
    if (!success)
    {
        return false;
    }

    double * psol = new double [n];
    for (int i = 0; i < n; i++)
    {
        psol[i] = 0.0;
    }
    ldl_ltsolve(L, Dy, psol, n); // L^Tx = y

    for (int i = 0; i < n; i++)
    {
        sol[i] = psol[pp[i]];
    }

    delete [] prhs;
    delete [] Lz;
    delete [] Dy;
    delete [] psol;

    return success;
}
