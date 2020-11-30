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

#ifndef __SGTELIB_SURROGATE_PRS__
#define __SGTELIB_SURROGATE_PRS__

#include "Surrogate.hpp"

namespace SGTELIB {

  /*--------------------------------------*/
  /*         Surrogate_PRS class        */
  /*--------------------------------------*/
  class Surrogate_PRS : public SGTELIB::Surrogate {

    /*--------------------------------------------------------*/
    /*  these members are defined in the Surrogate superclass */
    /*--------------------------------------------------------*/
    // int _p; // number of data points in X and Z
    // int _n; // dimension -- number of variables
    // int _m; // number of outputs (includes the objective)

  protected:

    int _q; // Nb of basis function
    SGTELIB::Matrix _M; // Monomes
    SGTELIB::Matrix _H; // Design matrix
    SGTELIB::Matrix _Ai; // Inverse of Ht*H
    SGTELIB::Matrix _alpha; // Coefficients

    virtual const SGTELIB::Matrix compute_design_matrix ( const SGTELIB::Matrix& Monomes, 
                                                          const SGTELIB::Matrix & Xs );

    // build model (private):
    virtual bool build_private (void) override;

    void predict_private ( const SGTELIB::Matrix & XXs,
                                 SGTELIB::Matrix * ZZs) override;


    // Compute metrics
    const SGTELIB::Matrix * get_matrix_Zvs (void) override;

    bool compute_alpha ( void );

  public:

    // Constructor
    Surrogate_PRS ( SGTELIB::TrainingSet & trainingset ,   
                    SGTELIB::Surrogate_Parameters param) ;

    // destructor:
    virtual ~Surrogate_PRS ( void );

    // Build the monome exponents
    static int get_nb_PRS_monomes(const int nvar, const int degree);
    static SGTELIB::Matrix get_PRS_monomes(const int nvar, const int degree);
    virtual void display_private ( std::ostream & out ) const override;

  };
}

#endif
