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

#ifndef __SURROGATE_FACTORY__
#define __SURROGATE_FACTORY__

#include "Defines.hpp"
#include "Exception.hpp"
#include "Surrogate.hpp"
#include "Surrogate_KS.hpp"
#include "Surrogate_CN.hpp"
#include "Surrogate_RBF.hpp"
#include "Surrogate_PRS.hpp"
#include "Surrogate_PRS_EDGE.hpp"
#include "Surrogate_PRS_CAT.hpp"
#include "Surrogate_Ensemble.hpp"
#include "Surrogate_LOWESS.hpp"
#include "Surrogate_Kriging.hpp"

namespace SGTELIB {

DLL_API SGTELIB::Surrogate * Surrogate_Factory ( SGTELIB::TrainingSet    & C,
                                         const std::string & s );

DLL_API SGTELIB::Surrogate * Surrogate_Factory ( SGTELIB::Matrix & X0,
                                         SGTELIB::Matrix & Z0,
                                         const std::string & s );

DLL_API void surrogate_delete ( SGTELIB::Surrogate * S );

}

#endif
