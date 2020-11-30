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

#ifndef __SGTELIB_SURROGATE_PRS_CAT__
#define __SGTELIB_SURROGATE_PRS_CAT__

#include "Surrogate_PRS.hpp"

namespace SGTELIB {

  /*--------------------------------------*/
  /*         Surrogate_PRS_CAT class      */
  /*--------------------------------------*/
  class Surrogate_PRS_CAT : public SGTELIB::Surrogate_PRS {

  protected:
    std::set<double> _cat; // Categories
    int _nb_cat; // Number of categories

    virtual const SGTELIB::Matrix compute_design_matrix ( const SGTELIB::Matrix& Monomes, 
                                                          const SGTELIB::Matrix & Xs ) override;
    // build model (private):
    virtual bool build_private (void) override;
    virtual bool init_private  (void) override;
  public:

    // Constructor
    Surrogate_PRS_CAT ( SGTELIB::TrainingSet & trainingset ,   
                        SGTELIB::Surrogate_Parameters param) ;

    // destructor:
    virtual ~Surrogate_PRS_CAT ( void );

    virtual void display_private ( std::ostream & out ) const override;

  };
}

#endif
