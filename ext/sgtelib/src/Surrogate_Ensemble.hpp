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

#ifndef __SGTELIB_SURROGATE_ENSEMBLE__
#define __SGTELIB_SURROGATE_ENSEMBLE__

#include "Surrogate.hpp"
#include "Surrogate_Factory.hpp"


//#include <time.h>


namespace SGTELIB {

  const double wta3_alpha = 0.05;
  const double wta3_beta  = -1;

  /*--------------------------------------*/
  /*         Surrogate_Ensemble class        */
  /*--------------------------------------*/
  class DLL_API Surrogate_Ensemble : public SGTELIB::Surrogate {

    /*--------------------------------------------------------*/
    /*  these members are defined in the Surrogate superclass */
    /*--------------------------------------------------------*/
    // int _p; // number of data points in X and Z
    // int _n; // dimension -- number of variables
    // int _m; // number of outputs (includes the objective)

  protected:

    int _kmax; // Nb of surrogates in the ensemble
    int _kready; // Nb of surrogates READY in the ensemble
    //SGTELIB::Matrix _W; // Weight vector
    std::vector<SGTELIB::Surrogate *>  _surrogates; // List des surrogates
    bool * _active; // Array of boolean. Is _active[k] is true if surrogate k is ready
                    // AND if there is a j such that W(k,j)!=0 
                    // ie: the weight in k is non null for at least one output
    double * _metric; // Value of the metric for the Ensemble

    // build model (private):
    virtual bool build_private (void) override;
    virtual bool init_private  (void) override;


    // Compute metrics
    virtual const SGTELIB::Matrix * get_matrix_Zhs (void) override;
    virtual const SGTELIB::Matrix * get_matrix_Shs (void) override;
    virtual const SGTELIB::Matrix * get_matrix_Zvs (void) override;

    void compute_W_by_select(void);
    void compute_W_by_wta1  (void);
    void compute_W_by_wta3  (void);

    // predict model (private):
    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZ ,
                                         SGTELIB::Matrix * std, 
                                         SGTELIB::Matrix * ei ,
                                         SGTELIB::Matrix * cdf ) override;
 
    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZ ) override;


  public:

    // Constructor
    Surrogate_Ensemble ( SGTELIB::TrainingSet & trainingset ,   
                         SGTELIB::Surrogate_Parameters param) ;

    /*
    Surrogate_Ensemble ( SGTELIB::TrainingSet & trainingset ,   
                         const std::string & s) ;
    */

    // destructor:
    virtual ~Surrogate_Ensemble ( void );

    virtual void display_private ( std::ostream & out ) const override;
    void display ( std::ostream & out , const int k ) const {_surrogates.at(k)->display(out);};

    // ==============================================//
    // Method for inspection of the basic surrogates //
    // ==============================================//

    // Test if basic model k is ready.
    bool is_ready (const int k) const;

    // Compute the boolean array _active
    void compute_active_models ( void ) ;
    // Check the weight vector
    bool check_weight_vector ( void ) const;



    // ==============================================//
    // Method to define the model_list //
    // ==============================================//
    void model_list_display        ( std::ostream & out );
    void model_list_preset         ( const std::string & preset );
    void model_list_remove_all ( void );
    void model_list_add ( const std::string & definition );

  };
}

#endif
