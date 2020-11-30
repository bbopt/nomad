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

#ifndef __SGTELIB_SURROGATE_CN__
#define __SGTELIB_SURROGATE_CN__

#include "Surrogate.hpp"

namespace SGTELIB {

  /*--------------------------------------*/
  /*         Surrogate_CN class        */
  /*--------------------------------------*/
  class Surrogate_CN : public SGTELIB::Surrogate {

  private:
  
    virtual bool build_private (void) override;
 
    virtual void predict_private ( const SGTELIB::Matrix & XXs,
                                         SGTELIB::Matrix * ZZs) override;

    // Compute metrics
    virtual const SGTELIB::Matrix * get_matrix_Zvs (void) override;
    virtual const SGTELIB::Matrix * get_matrix_Zhs (void) override;
    virtual const SGTELIB::Matrix * get_matrix_Svs (void) override;
    virtual const SGTELIB::Matrix * get_matrix_Shs (void) override;

    bool compute_cv_values (void);

  public:

    // Constructor
    Surrogate_CN ( SGTELIB::TrainingSet & trainingset ,   
                   SGTELIB::Surrogate_Parameters param) ;

    // destructor:
    virtual ~Surrogate_CN ( void );

    virtual void display_private ( std::ostream & out ) const override;

  };
}

#endif

