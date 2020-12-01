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

#ifndef __SGTELIB_SURROGATE_PARAMETERS__
#define __SGTELIB_SURROGATE_PARAMETERS__

#include "Defines.hpp"
#include "Kernel.hpp"
#include "Metrics.hpp"
#include <string>

namespace SGTELIB {

  /*--------------------------------------*/
  /*         Surrogate_Parameters class        */
  /*--------------------------------------*/
  class DLL_API Surrogate_Parameters  {

  protected:

    // Type of the model
    const model_t _type;

    // Degree
    int _degree; 
    SGTELIB::param_status_t _degree_status;
    // Kernel type 
    kernel_t _kernel_type;
    SGTELIB::param_status_t _kernel_type_status;
    // Kernel coef
    double _kernel_coef;
    SGTELIB::param_status_t _kernel_coef_status;
    // Use ridge
    double _ridge;
    SGTELIB::param_status_t _ridge_status;
    // Distance 
    distance_t _distance_type;
    SGTELIB::param_status_t _distance_type_status;
    // Ensemble weights
    SGTELIB::Matrix _weight; // Stores the weight values
    weight_t _weight_type; // Indicate the method (WTA1, WTA2, WTA3, or OPTIM)
    SGTELIB::param_status_t _weight_status;
    // Kriging covariance coefficients
    SGTELIB::Matrix _covariance_coef;
    SGTELIB::param_status_t _covariance_coef_status;
    // Metric used for weight calculation
    metric_t _metric_type;
    // Preset
    std::string _preset;
    // Output file
    std::string _output;
    // Optimization budget
    int _budget;

    // Nb of parameters that are optimized
    int _nb_parameter_optimization;

  public:

    // Constructors
    explicit Surrogate_Parameters ( const model_t mt);
    explicit Surrogate_Parameters ( const std::string & s);
    //Surrogate_Parameters ( const SGTELIB::Surrogate_Parameters & p );

    // Defaults
    void set_defaults (void);

    // Read strings
    bool authorized_field  (const std::string & field) const;
    static bool authorized_optim  (const std::string & field);
    void read_string (const std::string & model_description);
    static SGTELIB::model_t read_model_type ( const std::string & model_description);
    static std::string to_standard_field_name (const std::string& field);

    // Check    
    void check (void);

    // Get
    model_t         get_type            (void) const {return _type;};
    int             get_degree          (void) const {return _degree;};
    kernel_t        get_kernel_type     (void) const {return _kernel_type;};
    double          get_kernel_coef     (void) const {return _kernel_coef;};
    double          get_ridge           (void) const {return _ridge;};
    SGTELIB::Matrix get_weight          (void) const {return _weight;};
    weight_t        get_weight_type     (void) const {return _weight_type;};
    metric_t        get_metric_type     (void) const {return _metric_type;};
    std::string     get_metric_type_str (void) const {return SGTELIB::metric_type_to_str(_metric_type);};
    distance_t      get_distance_type   (void) const {return _distance_type;};
    std::string     get_preset          (void) const {return _preset;};
    std::string     get_output          (void) const {return _output;};
    SGTELIB::Matrix get_covariance_coef (void) const {return _covariance_coef;};
    int             get_budget          (void) const {return _budget;};
    // Get the distance type (return OPTIM if the distance type has to be optimized).
    std::string get_distance_type_str (void) const {
      if (_distance_type_status==SGTELIB::STATUS_OPTIM) return "OPTIM";
      else return SGTELIB::distance_type_to_str(_distance_type);
    };


    // Set
    void set_kernel_coef     ( const double v          ) { _kernel_coef = v;      };
    void set_weight_type     ( const weight_t wt       ) { _weight_type = wt;     };
    void set_weight          ( const SGTELIB::Matrix& W ) { _weight = W;           };
    void update_covariance_coef ( const int nvar );
    //void set_distance_type   ( distance_t dt ) { _distance_type = dt; };
    //void set_degree          ( int d         ) { _degree = d;         };

    // destructor:
    virtual ~Surrogate_Parameters ( void );
    std::string get_string(void) const;
    std::string get_short_string(void) const;
    void display ( std::ostream & out ) const;

    // Information and optimization
    SGTELIB::Matrix get_x ( void );
    double get_x_penalty ( void );
    bool check_x ( void );
    void display_x ( std::ostream & out );
    void set_x ( const SGTELIB::Matrix& X );
    void get_x_bounds ( SGTELIB::Matrix * LB ,
                        SGTELIB::Matrix * UB ,
                        SGTELIB::param_domain_t * domain,
                        bool * logscale );
    int get_nb_parameter_optimization (void) { return _nb_parameter_optimization; };

  };
}

#endif
