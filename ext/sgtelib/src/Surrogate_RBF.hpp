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

#ifndef __SGTELIB_SURROGATE_RBF__
#define __SGTELIB_SURROGATE_RBF__

#include "Surrogate.hpp"

namespace SGTELIB {

  /*--------------------------------------*/
  /*         Surrogate_RBF class        */
  /*--------------------------------------*/
  class Surrogate_RBF : public SGTELIB::Surrogate {

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
    int _q; // Nb of basis function
    int _qrbf; // Nb of RBF basis function
    int _qprs; // Nb of PRS basis function
    SGTELIB::Matrix _H; // Design matrix
    SGTELIB::Matrix _HtH; // H'*H
    SGTELIB::Matrix _HtZ; // H'*Z
    SGTELIB::Matrix _Ai; // inverse of H or Ht*H+r*J
    SGTELIB::Matrix _Alpha; // Coefficients

    std::list<int> _selected_kernel;


    /*--------------------------------------*/
    /*          Building methods            */
    /*--------------------------------------*/
    const SGTELIB::Matrix compute_design_matrix ( const SGTELIB::Matrix & XXs , const bool constraints );

    /*--------------------------------------*/
    /*          Build model                 */
    /*--------------------------------------*/
    bool select_kernels ( void ); 
    virtual bool build_private (void) override;
    virtual bool init_private  (void) override;
    //SGTELIB::Matrix get_bumpiness (void);

    /*--------------------------------------*/
    /*          predict                     */
    /*--------------------------------------*/ 
    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZs) override;

    /*--------------------------------------*/
    /*          Compute matrices            */
    /*--------------------------------------*/
    virtual const SGTELIB::Matrix * get_matrix_Zvs (void) override;

  public:

    /*--------------------------------------*/
    /*          Constructor                 */
    /*--------------------------------------*/
    Surrogate_RBF ( SGTELIB::TrainingSet & trainingset ,   
                     SGTELIB::Surrogate_Parameters param) ;

    /*--------------------------------------*/
    /*          Destructor                  */
    /*--------------------------------------*/
    virtual ~Surrogate_RBF ( void );

    /*--------------------------------------*/
    /*          Misc                        */
    /*--------------------------------------*/
    virtual void display_private ( std::ostream & out ) const override;

  };
}

#endif
