/*-------------------------------------------------------------------------------------*/
/*  sgtelib - A surrogate model library for derivative-free optimization               */
/*  Version 2.0.2                                                                      */
/*                                                                                     */
/*  Copyright (C) 2012-2017  Sebastien Le Digabel - Ecole Polytechnique, Montreal      */ 
/*                           Bastien Talgorn - McGill University, Montreal             */
/*                                                                                     */
/*  Author: Bastien Talgorn                                                            */
/*  email: bastientalgorn@fastmail.com                                                 */
/*                                                                                     */
/*  This program is free software: you can redistribute it and/or modify it under the  */
/*  terms of the GNU Lesser General Public License as published by the Free Software   */
/*  Foundation, either version 3 of the License, or (at your option) any later         */
/*  version.                                                                           */
/*                                                                                     */
/*  This program is distributed in the hope that it will be useful, but WITHOUT ANY    */
/*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    */
/*  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   */
/*                                                                                     */
/*  You should have received a copy of the GNU Lesser General Public License along     */
/*  with this program. If not, see <http://www.gnu.org/licenses/>.                     */
/*                                                                                     */
/*  You can find information on sgtelib at https://github.com/bastientalgorn/sgtelib   */
/*-------------------------------------------------------------------------------------*/

#ifndef __SGTELIB_MATRIX__
#define __SGTELIB_MATRIX__

#include <set>
#include <fstream>
#include <vector>
#include <list>
#include <climits>
#include <algorithm>
#include "Surrogate_Utils.hpp"
#include "Exception.hpp"
#include "Defines.hpp"

namespace SGTELIB {


  class DLL_API Matrix {

  private:

    std::string _name;

    int _nbRows; // nbRows x nbCols matrix
    int _nbCols;

    double ** _X;

  public:

    // constructor 1:
    explicit Matrix ( const std::string & name ,
             int                 nbRows    ,
             int                 nbCols      );

    // constructor 2:
    explicit Matrix ( const std::string & name ,
             int                 nbRows    ,
             int                 nbCols    ,
             double           ** A      );

    // constructor 3:
    explicit Matrix ( const std::string & file_name );

    // constructor 4:
    Matrix ( void );

    // constructor 5:
    Matrix ( double );

    // copy constructor:
    Matrix ( const Matrix & );

    // affectation operator:
    Matrix & operator = ( const Matrix & A );

    //Matrix & operator * ( const Matrix & B);



    // destructor:
    virtual ~Matrix ( void );

    // fill the matrix randomly:
    void set_random ( double l , double u , bool round = false );
    void fill (double v);

    // add rows:
    void add_rows ( const Matrix & X );
    void add_cols ( const Matrix & X );
    void add_row  ( const double * row );
    void add_rows ( const int p); // add empty rows
    void add_cols ( const int p); // add empty cols
    void remove_rows ( const int p); // remove last rows

    // GET methods:
    inline int get_nb_rows ( void ) const { return _nbRows; }
    inline int get_nb_cols ( void ) const { return _nbCols; }
    inline int get_numel   ( void ) const { return _nbRows*_nbCols; }

    // access to element (k)
    double get ( const int k ) const; 
    // access to element (i,j)
    inline double get ( const int i , const int j ) const {
      #ifdef SGTELIB_DEBUG
        if ( i < 0 || i >= _nbRows || j < 0 || j >= _nbCols ){
          display(std::cout);
          std::cout << "Error: try to access (" << i << "," << j << ") while dim is [" << _nbRows << "," << _nbCols << "]\n";
          throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::get(i,j): bad index" );
        }
      #endif
      return _X[i][j];
    }//

    const double & operator [] ( int k ) const;
    double & operator [] ( int k );



    SGTELIB::Matrix get ( const std::list<int> & list_cols , 
                          const std::list<int> & list_rows) const;

    SGTELIB::Matrix get_row (const int i) const;
    SGTELIB::Matrix get_col (const int i) const;

    SGTELIB::Matrix get_rows (const std::list<int> & list_rows) const;
    SGTELIB::Matrix get_cols (const std::list<int> & list_cols) const;

    SGTELIB::Matrix get_rows (const int i1, const int i2) const;
    SGTELIB::Matrix get_cols (const int i1, const int i2) const;
  
    void swap_rows (const int i1, const int i2);


    // count the number of different values in column j:
    int get_nb_diff_values ( int j ) const;

    // get the constant columns (constant variables):
    void get_fix_columns ( std::list<int> & fix_col ) const;

    // check symmetry
    bool is_sym ( void ) const;

    // SET methods:
    inline void set_name ( const std::string & name ) { _name = name; }
    inline std::string get_name ( void ) const { return _name; }

    void set     (const int i , const int j , const double d );
    void set_row (const SGTELIB::Matrix & T , const int i); // T is row vector
    void set_col (const SGTELIB::Matrix & T , const int j); // T is col vector
    void set_row (const double v , const int i); // T is row vector
    void set_col (const double v , const int j); // T is col vector
  
    // swap terms (i1,j1) and (i2,j2)
    void swap (const int i1 , const int j1 , const int i2 , const int j2 );

    // Multiply row
    void multiply_row (const double v , const int i); // T is row vector
    void multiply_col (const double v , const int j); // T is col vector

    // Inverse
    SGTELIB::Matrix SVD_inverse ( void ) const;

    // Inverse the diagonal terms
    SGTELIB::Matrix diag_inverse ( void ) const;

    // Build vector from a double*
    static SGTELIB::Matrix row_vector ( const double * v,
                                        const int n     );

    static SGTELIB::Matrix col_vector ( const double * v,
                                        const int n     );

    // Transpose
    SGTELIB::Matrix transpose ( void ) const;

    // Diag
    SGTELIB::Matrix diag (void ) const;

    // Trace
    double trace ( void ) const;

    // Rmse
    double rmse ( void ) const;

    // Norm
    double norm ( void ) const;
    double normsquare ( void ) const;
    void normalize_cols ( void );
    SGTELIB::Matrix col_norm ( const norm_t nt ) const;

    // Sum
    double sum ( void ) const;
    SGTELIB::Matrix sum ( const int direction ) const;

    // Mean
    double mean ( void ) const;

    // Count (number of non null values)
    int count ( void ) const;


    // Product
    static SGTELIB::Matrix product ( const SGTELIB::Matrix & A,
                                     const SGTELIB::Matrix & B);

    static SGTELIB::Matrix product ( const SGTELIB::Matrix & A,
                                     const SGTELIB::Matrix & B,
                                     const SGTELIB::Matrix & C);

    static SGTELIB::Matrix product ( const SGTELIB::Matrix & A,
                                     const SGTELIB::Matrix & B,
                                     const SGTELIB::Matrix & C,
                                     const SGTELIB::Matrix & D);

    inline void product ( const int i , const int j , const double v){ _X[i][j]*=v; };

    // Subset product, multiply
    // the p first rows and q first columns of A
    // with the q first rows and r first columns of B.
    // Result is a matrix of size p/r.
    static SGTELIB::Matrix subset_product (const SGTELIB::Matrix & A,
                                           const SGTELIB::Matrix & B,
                                           int p,
                                           int q,
                                           int r);

    static SGTELIB::Matrix diagA_product ( const SGTELIB::Matrix & A,
                                           const SGTELIB::Matrix & B);

    static SGTELIB::Matrix diagB_product ( const SGTELIB::Matrix & A,
                                           const SGTELIB::Matrix & B);

    static SGTELIB::Matrix transposeA_product ( const SGTELIB::Matrix & A,
                                                const SGTELIB::Matrix & B);


    static SGTELIB::Matrix hadamard_product ( const SGTELIB::Matrix & A,
                                              const SGTELIB::Matrix & B);

    static SGTELIB::Matrix hadamard_square  ( const SGTELIB::Matrix & A );

    static SGTELIB::Matrix hadamard_sqrt    ( const SGTELIB::Matrix & A );

    static SGTELIB::Matrix hadamard_power    ( const SGTELIB::Matrix & A ,
                                               const double e            );

    void hadamard_inverse ( void );
    void hadamard_sqrt    ( void );
    void hadamard_square  ( void );

    // Addition
    static SGTELIB::Matrix add ( const SGTELIB::Matrix & A,
                                 const SGTELIB::Matrix & B);

    // Add to the matrix itself
    void add ( const SGTELIB::Matrix & B);
    void add ( const int i , const int j , const double v){ _X[i][j]+=v; };

    // Add and fill with 0 (add two matrices of different sizes)
    static SGTELIB::Matrix add_fill ( const SGTELIB::Matrix & A,
                                      const SGTELIB::Matrix & B);

    // Substract
    static SGTELIB::Matrix sub ( const SGTELIB::Matrix & A,
                                 const SGTELIB::Matrix & B);

    void sub ( const SGTELIB::Matrix & B);
  
    // Identity matrix
    static SGTELIB::Matrix identity ( const int n );

    // ones matrix
    static SGTELIB::Matrix ones ( const int nbRows , const int nbCols );

    // Lines random permutation
    SGTELIB::Matrix random_line_permutation ( void ) const;

    // Rank of the values
    SGTELIB::Matrix rank ( void ) const;

    // Conjugate gradient
    static SGTELIB::Matrix conjugate_solve ( const SGTELIB::Matrix & A ,
                                             const SGTELIB::Matrix & b ,
                                             const SGTELIB::Matrix & x0 ,
                                             const double tol);

    // LU factorization
    SGTELIB::Matrix lu_inverse ( void ) const;
    SGTELIB::Matrix lu_inverse ( double * det ) const;

    // Cholesky
    SGTELIB::Matrix cholesky ( void ) const;
    SGTELIB::Matrix cholesky_inverse ( double * det ) const;
    SGTELIB::Matrix cholesky_inverse ( void ) const;
    static SGTELIB::Matrix cholesky_solve ( const SGTELIB::Matrix & A ,
                                            const SGTELIB::Matrix & b );



    // Triangular matrix
    static SGTELIB::Matrix tril_inverse (const SGTELIB::Matrix & L );

    static SGTELIB::Matrix triu_solve ( const SGTELIB::Matrix & U ,
                                        const SGTELIB::Matrix & b );

    static SGTELIB::Matrix tril_solve ( const SGTELIB::Matrix & L ,
                                        const SGTELIB::Matrix & b );
    // SVD decomposition:
    bool SVD_decomposition ( std::string & error_msg ,
                             SGTELIB::Matrix * &MAT_U,  // OUT, nbRows x nbCols
                             SGTELIB::Matrix * &MAT_W,  // OUT, nbCols x nbCols, diagonal
                             SGTELIB::Matrix * &MAT_V,  // OUT, nbCols x nbCols
                             int           max_mpn = 1500 ) const;

    bool SVD_decomposition ( std::string & error_msg      ,
                             double     ** U              ,  // OUT, nbRows x nbCols
                             double      * W              ,  // OUT, nbCols x nbCols, diagonal
                             double     ** V              ,  // OUT, nbCols x nbCols
                             int           max_mpn = 1500 ) const;

    // Projection matrix for linear over-determined models
    static SGTELIB::Matrix get_matrix_P     ( const SGTELIB::Matrix & Ai,
                                              const SGTELIB::Matrix & H );

    static double          get_trace_P      ( const SGTELIB::Matrix & Ai,
                                              const SGTELIB::Matrix & H );

    static SGTELIB::Matrix get_matrix_PZs   ( const SGTELIB::Matrix & Ai,
                                              const SGTELIB::Matrix & H ,
                                              const SGTELIB::Matrix & Zs);

    static SGTELIB::Matrix get_matrix_dPiPZs( const SGTELIB::Matrix & Ai,
                                              const SGTELIB::Matrix & H ,
                                              const SGTELIB::Matrix & Zs);

    static SGTELIB::Matrix get_matrix_dPiPZs( const SGTELIB::Matrix & Ai,
                                              const SGTELIB::Matrix & H ,
                                              const SGTELIB::Matrix & Zs ,
                                              const SGTELIB::Matrix & ALPHA);


    static SGTELIB::Matrix get_matrix_dPi   ( const SGTELIB::Matrix & Ai,
                                              const SGTELIB::Matrix & H );

    // Min / Max
    double max (void);
    double min (void);
    static SGTELIB::Matrix max ( const SGTELIB::Matrix & A , 
                                 const SGTELIB::Matrix & B );
    static SGTELIB::Matrix min ( const SGTELIB::Matrix & A , 
                                 const SGTELIB::Matrix & B );

    // Get min index
    int get_max_index     ( void );
    int get_min_index     ( void );
    int get_min_index_row ( const int i );
    int get_min_index_col ( const int j );

    // display:
    void display      ( std::ostream & out       ) const;
    void display_short( std::ostream & out       ) const;
    void write        ( const std::string  & file_name ) const;
    void display_size ( std::ostream & out       ) const;

    // import data in plain format
    static SGTELIB::Matrix import_data   ( const std::string & file_name );
    static SGTELIB::Matrix string_to_matrix ( std::string s );
    static SGTELIB::Matrix string_to_row ( const std::string & s , int nbCols = 0 );

    // distances
    static SGTELIB::Matrix get_distances_norm1   ( const SGTELIB::Matrix & A , 
                                                   const SGTELIB::Matrix & B );
    static SGTELIB::Matrix get_distances_norm2   ( const SGTELIB::Matrix & A , 
                                                   const SGTELIB::Matrix & B );
    static SGTELIB::Matrix get_distances_norminf ( const SGTELIB::Matrix & A , 
                                                   const SGTELIB::Matrix & B );

    int find_row (SGTELIB::Matrix & R);

    // nan 
    bool has_nan (void) const;
    bool has_inf (void) const;
    void replace_nan (double d);

    // Generate poll directions
    static SGTELIB::Matrix get_poll_directions ( const SGTELIB::Matrix& scaling,
                                                 const SGTELIB::param_domain_t * domain,
                                                 double psize );


    SGTELIB::Matrix LUPinverse (void);

  };
}

SGTELIB::Matrix operator * (const SGTELIB::Matrix & A , const double v           );
SGTELIB::Matrix operator * (const double v            , const SGTELIB::Matrix & A);
SGTELIB::Matrix operator * (const SGTELIB::Matrix & A , const SGTELIB::Matrix & B);
SGTELIB::Matrix operator + (const SGTELIB::Matrix & A , const SGTELIB::Matrix & B);
SGTELIB::Matrix operator + (const SGTELIB::Matrix & A , const double v           );
SGTELIB::Matrix operator + (const double v            , const SGTELIB::Matrix & A);
SGTELIB::Matrix operator - (const SGTELIB::Matrix & A , const SGTELIB::Matrix & B);
SGTELIB::Matrix operator - (const SGTELIB::Matrix & A , const double v);
SGTELIB::Matrix operator - (const double v            , const SGTELIB::Matrix & A);
SGTELIB::Matrix operator - (const SGTELIB::Matrix & A);
SGTELIB::Matrix operator / (const SGTELIB::Matrix & A , const double v           );

#endif
