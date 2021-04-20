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

#ifndef __SGTELIB_TRAININGSET__
#define __SGTELIB_TRAININGSET__

#include "Matrix.hpp"
#include "Defines.hpp"
namespace SGTELIB {

  /*--------------------------------------*/
  /*             TrainingSet class          */
  /*--------------------------------------*/
  class DLL_API TrainingSet {

  private:

    int _p; // number of data points in X and Z
    const int _n; // dimension -- number of variables
    const int _m; // number of outputs (includes the objective)
    bool _ready; // True if data have been processed and are ready to be used 

    // Output type
    SGTELIB::bbo_t * _bbo;
    bool _bbo_is_def;
    int _j_obj; // Index of the objective

    // f_min
    double _f_min;
    double _fs_min;
    int _i_min; // Index of the point where f_min is reached
    
    // data points:
    SGTELIB::Matrix _X; // p x n
    SGTELIB::Matrix _Z; // p x m

    // scaled matrices
    SGTELIB::Matrix _Xs; // p x n
    SGTELIB::Matrix _Zs; // p x m

    // Distance matrix
    SGTELIB::Matrix _Ds; // p x p

    // Nb of varying data
    int _nvar; // Nb of varying input
    int _mvar; // Nb of varying output
    int _pvar; // Nb of different points

    // Data
    double * _X_lb;
    double * _X_ub;
    double * _X_scaling_a;
    double * _X_scaling_b;
    double * _X_mean;
    double * _X_std;
    int    * _X_nbdiff;
    int      _X_nbdiff1;
    int      _X_nbdiff2;
    double * _Z_lb;
    double * _Z_ub;
    double * _Z_replace;
    double * _Z_scaling_a;
    double * _Z_scaling_b;
    double * _Z_mean;
    double * _Z_std;
    double * _Zs_mean;
    int    * _Z_nbdiff;

    // scaled bounds:
    //double * _X_lb_scaled;
    //double * _X_ub_scaled;

    // Mean distance between points 
    double _Ds_mean;

    // private affectation operator:
    TrainingSet & operator = ( const TrainingSet & );

    // Data preparation 

    static void compute_nbdiff   ( const SGTELIB::Matrix & MAT , 
                                   int * nbdiff ,
                                   int & njvar );

    void compute_bounds          (void);
    void compute_mean_std        (void);
    void compute_nvar_mvar       (void);
    void compute_scaling         (void);
    void compute_Ds              (void);
    void compute_scaled_matrices (void);
    void compute_f_min           (void);
    bool check_singular_data     (void);

    // FORBIDEN copy constructor:
    TrainingSet ( const TrainingSet & );


  public:
    // constructor 1:
    TrainingSet ( const SGTELIB::Matrix & X ,
                  const SGTELIB::Matrix & Z );

    // destructor:
    virtual ~TrainingSet ( void );

    // Define the bbo types
    void set_bbo_type (const std::string & s);

    // construct/process the data of the TrainingSet
    void build ( void );

    // check ready
    bool is_ready ( void ) const {return _ready;};
    void check_ready (const std::string & s) const;
    void check_ready (const std::string & file, const std::string & function, const int & i) const;
    void check_ready (void) const;

    // add points:
    bool add_points ( const SGTELIB::Matrix & Xnew ,
                      const SGTELIB::Matrix & Znew  );
      
    bool partial_reset_and_add_points( const SGTELIB::Matrix & Xnew ,
    const SGTELIB::Matrix & Znew  );
      
    bool add_point  ( const double * xnew ,
                      const double * znew  );

    // scale and unscale:
    void   X_scale      ( double * x ) const;
    void   X_unscale    ( double * y ) const;

    double X_scale      ( double x , int var_index ) const;
    double X_unscale    ( double y , int var_index ) const;

    void   Z_scale      ( double * z ) const;
    void   Z_unscale    ( double * w ) const;
    
    double Z_scale      ( double z , int output_index ) const;
    double Z_unscale    ( double w , int output_index ) const;
    double ZE_unscale   ( double w , int output_index ) const;

    void   X_scale      ( SGTELIB::Matrix &X);
    void   Z_unscale    ( SGTELIB::Matrix *Z);
    void   ZE_unscale   ( SGTELIB::Matrix *ZE);
    SGTELIB::Matrix Z_unscale  ( const SGTELIB::Matrix & Z  );
    SGTELIB::Matrix ZE_unscale ( const SGTELIB::Matrix & ZE );

    // Get data
    double get_Xs       ( const int i , const int j ) const;
    double get_Zs       ( const int i , const int j ) const;
    void   get_Xs       ( const int i , double * x  ) const;
    void   get_Zs       ( const int i , double * z  ) const;
    double get_Zs_mean  ( const int j ) const;
    int    get_X_nbdiff ( const int i ) const;
    int    get_Z_nbdiff ( const int j ) const;
    int    get_X_nbdiff1 (void) const { check_ready(); return _X_nbdiff1; };
    int    get_X_nbdiff2 (void) const { check_ready(); return _X_nbdiff2; };
    double get_Ds       ( const int i1, const int i2) const;
    const SGTELIB::Matrix get_X_nbdiff ( void ) const;
    

    SGTELIB::Matrix get_distances ( const SGTELIB::Matrix & A ,
                                    const SGTELIB::Matrix & B , 
                                    const distance_t dt = SGTELIB::DISTANCE_NORM2 ) const;

    double get_d1_over_d2 ( const SGTELIB::Matrix & XXs ) const;
    double get_d1         ( const SGTELIB::Matrix & XXs ) const;
    SGTELIB::Matrix get_exclusion_area_penalty ( const SGTELIB::Matrix & XXs , const double tc ) const;
    SGTELIB::Matrix get_distance_to_closest ( const SGTELIB::Matrix & XXs ) const;
 

    // Get basic information
    int get_nb_points      ( void ) const { return _p; };
    int get_input_dim      ( void ) const { return _n; };
    int get_output_dim     ( void ) const { return _m; };
    int get_nvar           ( void ) const { check_ready(); return _nvar; };
    int get_mvar           ( void ) const { check_ready(); return _mvar; };
    int get_pvar           ( void ) const { check_ready(); return _pvar; };

    double get_fs_min      ( void ) const { check_ready(); return _fs_min; };
    double get_f_min       ( void ) const { check_ready(); return _f_min;  };
    int    get_i_min       ( void ) const { check_ready(); return _i_min;  };
    double get_Ds_mean     ( void ) const { check_ready(); return _Ds_mean; };
    int    get_j_obj       ( void ) const { check_ready(); return _j_obj; };
    //double get_Ds_min      ( void ) const ;

    SGTELIB::bbo_t get_bbo ( int j) const { check_ready(); return _bbo[j]; };

    double get_X_scaling_a ( int j) const { check_ready(); return _X_scaling_a[j]; };

    // Return the design matrix
    const SGTELIB::Matrix & get_matrix_Xs ( void ) const { check_ready(); return _Xs; };
    const SGTELIB::Matrix & get_matrix_Zs ( void ) const { check_ready(); return _Zs; };
    const SGTELIB::Matrix & get_matrix_Ds ( void ) const { check_ready(); return _Ds; };

    const SGTELIB::Matrix & get_matrix_X ( void ) const { return _X; };
      
    // display:
    void display ( std::ostream & out ) const;
    void info (void) const;

    // Selection of points
    std::list<int> select_greedy ( const SGTELIB::Matrix & X,
                                   const int imin,
                                   const int pS,
                                   const double lambda0, 
                                   const distance_t dt = SGTELIB::DISTANCE_NORM2);

  };
}

#endif
