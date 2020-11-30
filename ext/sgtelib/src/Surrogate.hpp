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

#ifndef __SGTELIB_SURROGATE__
#define __SGTELIB_SURROGATE__

#include "Matrix.hpp"
#include "TrainingSet.hpp"
#include "Kernel.hpp"
#include "Surrogate_Parameters.hpp"
#include <map>
namespace SGTELIB {

  /*--------------------------------------*/
  /*             Surrogate class          */
  /*--------------------------------------*/
  class DLL_API Surrogate {

  // Surrogate Ensemble is a friend, so that it can access to private and protected 
  // prediction methods of other derived classed of Surrogate_Ensemble
  friend class Surrogate_Ensemble;

  protected:

    // TrainingSet containing the data 
    // (may be shared between several surrogates)
    SGTELIB::TrainingSet & _trainingset;
    // Parameters
    SGTELIB::Surrogate_Parameters _param;

    // Input dim
    const int _n;
    // Ouptut dim
    const int _m;

    // Number of data points of the training set
    int _p_ts;    
    int _p_ts_old;    
    // Number of data points used in the model
    int _p; 
    int _p_old;

    // Is the surrogate ready to perform predictions ?
    bool _ready; 

    // Predictions
    // _Zhs: Prediction of the model in the training points
    // (used to compute emax, rmse, eotp, linv)
    SGTELIB::Matrix * _Zhs; 
    SGTELIB::Matrix * _Shs; // Predictive std on the training points

    // _Zvs: Cross-Validation prediction of the model in the training points
    // (used to compute rmsecv, oecv)
    SGTELIB::Matrix * _Zvs; 
    SGTELIB::Matrix * _Svs; // Cross-validation std on the training points

    // List of points used to build the model
    std::list<int> _selected_points;

    // Map to store all the metrics
    std::map< metric_t , SGTELIB::Matrix > _metrics;

    // psize_max : Larger value of psize that led to a success
    // in the previous parameter optimization.
    double _psize_max; 

    // Output stream
    std::ofstream _out;
    bool _display;

    // private affectation operator:
    // Surrogate & operator = ( const Surrogate & );

    // build model (private):
    virtual bool build_private (void) = 0;
    virtual bool init_private  (void);

    // Compute metrics
    bool compute_metric ( const metric_t mt );
    virtual void compute_metric_linv (void);

    // Function used to compute "_metric_oe" and "_metric_oecv"
    SGTELIB::Matrix compute_order_error  ( const SGTELIB::Matrix & Zpred );
    double compute_aggregate_order_error ( const SGTELIB::Matrix & Zpred );

    SGTELIB::Matrix compute_efi( const SGTELIB::Matrix & Zs,
                                 const SGTELIB::Matrix & Ss );

    SGTELIB::Matrix compute_fh ( const SGTELIB::Matrix & Zs );

    // predict model (private):
    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZs,
                                         SGTELIB::Matrix * std, 
                                         SGTELIB::Matrix * ei ,
                                         SGTELIB::Matrix * cdf ); 
 
    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZs) = 0; 

    // Display private 
    virtual void display_private ( std::ostream & out ) const = 0;

    // get matrices (these matrices are unscaled before being returned)
    // (That's why these functions cant be public)
    const SGTELIB::Matrix get_matrix_Xs (void);
    const SGTELIB::Matrix get_matrix_Zs (void);
    const SGTELIB::Matrix get_matrix_Ds (void);

    // Compute scaled data
    // Compute the cross-validation matrix
    virtual const SGTELIB::Matrix * get_matrix_Zvs (void) = 0;
    virtual const SGTELIB::Matrix * get_matrix_Zhs (void);
    virtual const SGTELIB::Matrix * get_matrix_Shs (void);
    virtual const SGTELIB::Matrix * get_matrix_Svs (void);

  public:
    // constructor:

    Surrogate ( SGTELIB::TrainingSet & trainingset,
                const SGTELIB::Surrogate_Parameters& param);

    Surrogate ( SGTELIB::TrainingSet & trainingset,
                const SGTELIB::model_t& mt );

    Surrogate ( SGTELIB::TrainingSet & trainingset,
                const std::string & s );

    // destructor:
    virtual ~Surrogate ( void );

    void reset_metrics ( void );
    void info ( void ) const ;

    // check ready
    void check_ready (const std::string & s) const;
    void check_ready (const std::string & file, const std::string & function, const int & i) const;
    void check_ready (void) const;

    // Get metrics
    SGTELIB::Matrix get_metric (SGTELIB::metric_t mt);
    double          get_metric (SGTELIB::metric_t mt , int j);
    bool is_defined(const SGTELIB::metric_t mt);
    bool is_defined(const SGTELIB::metric_t mt, const int j);

    // construct:
    bool build (void);

    // predict:
    void predict ( const SGTELIB::Matrix & XX ,
                         SGTELIB::Matrix * ZZ , // nb : ZZ is a ptr
                         SGTELIB::Matrix * std, 
                         SGTELIB::Matrix * ei , 
                         SGTELIB::Matrix * cdf); 

    void predict ( const SGTELIB::Matrix & XX ,
                         SGTELIB::Matrix * ZZ ); 

    // Compute unscaled data
    const SGTELIB::Matrix get_matrix_Zh (void);
    const SGTELIB::Matrix get_matrix_Sh (void);
    const SGTELIB::Matrix get_matrix_Zv (void);
    const SGTELIB::Matrix get_matrix_Sv (void);

    // add points:
    bool add_points ( const SGTELIB::Matrix & Xnew ,
                      const SGTELIB::Matrix & Znew  );
    bool add_point  ( const double * xnew ,
                      const double * znew  );

    // exclusion_area
    SGTELIB::Matrix get_exclusion_area_penalty ( const SGTELIB::Matrix & XX , const double tc ) const;
    SGTELIB::Matrix get_distance_to_closest ( const SGTELIB::Matrix & XX ) const;

    bool is_ready                (void) const {return _ready;};
    void display_trainingset     (void) const {_trainingset.build();_trainingset.display(std::cout);};
    SGTELIB::model_t get_type    (void) const {return _param.get_type();};
    std::string get_string       (void) const {return _param.get_string();};
    std::string get_short_string (void) const {return _param.get_short_string();};
    void display                 ( std::ostream & out ) const;
    const SGTELIB::Surrogate_Parameters get_param (void) const {return _param;};

    // Set:
    void set_kernel_coef (double v)    { _param.set_kernel_coef(v); };

    // Parameter optimization
    bool optimize_parameters ( void );
    double eval_objective ( void );

  };
}

#endif

