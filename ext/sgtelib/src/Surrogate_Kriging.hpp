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

#ifndef __SGTELIB_SURROGATE_KRIGING__
#define __SGTELIB_SURROGATE_KRIGING__

#include "Surrogate.hpp"

namespace SGTELIB {

  /*--------------------------------------*/
  /*         Surrogate_Kriging class        */
  /*--------------------------------------*/
  class Surrogate_Kriging : public SGTELIB::Surrogate {

    /*--------------------------------------------------------*/
    /*  these members are defined in the Surrogate superclass */
    /*--------------------------------------------------------*/
    // int _p; // number of data points in X and Z
    // int _n; // dimension -- number of variables
    // int _m; // number of outputs (includes the objective)

  private:

    /*--------------------------------------*/
    /*          Attributes                  */
    /*--------------------------------------*/
    SGTELIB::Matrix _R; // Covariance Matrix
    SGTELIB::Matrix _Ri; // Inverte of _R
    SGTELIB::Matrix _H; // Polynomial terms
    SGTELIB::Matrix _alpha;
    SGTELIB::Matrix _beta;
    SGTELIB::Matrix _var;
    double _detR;

    /*--------------------------------------*/
    /*          Building methods            */
    /*--------------------------------------*/
    const SGTELIB::Matrix compute_covariance_matrix ( const SGTELIB::Matrix & XXs ); 

    /*--------------------------------------*/
    /*          Build model                 */
    /*--------------------------------------*/
    virtual bool build_private (void) override;
    virtual bool init_private  (void) override;
    //void init_covariance_coef (void);
    virtual void compute_metric_linv (void) override;
    bool compute_cv_values (void);

    /*--------------------------------------*/
    /*          predict                     */
    /*--------------------------------------*/ 
    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZs) override;

    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZs,
                                         SGTELIB::Matrix * std, 
                                         SGTELIB::Matrix * ei ,
                                         SGTELIB::Matrix * cdf ) override;

    /*--------------------------------------*/
    /*          Compute matrices            */
    /*--------------------------------------*/
    virtual const SGTELIB::Matrix * get_matrix_Zvs (void) override;
    virtual const SGTELIB::Matrix * get_matrix_Svs (void) override;

  public:

    /*--------------------------------------*/
    /*          Constructor                 */
    /*--------------------------------------*/
    Surrogate_Kriging ( SGTELIB::TrainingSet & trainingset ,   
                        SGTELIB::Surrogate_Parameters param) ;

    /*--------------------------------------*/
    /*          Destructor                  */
    /*--------------------------------------*/
    virtual ~Surrogate_Kriging ( void );

    /*--------------------------------------*/
    /*          Misc                        */
    /*--------------------------------------*/
    virtual void display_private ( std::ostream & out ) const override;

  };
}

#endif
