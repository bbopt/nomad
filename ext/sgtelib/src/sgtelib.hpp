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

// SGTELIB 2014-02-07

#ifndef __SGTELIB__
#define __SGTELIB__

#include "Surrogate_Factory.hpp"
#include "Tests.hpp"
#include "Matrix.hpp"
#include "Defines.hpp"
#include "Surrogate_Utils.hpp"
#include "sgtelib_help.hpp"

namespace SGTELIB {
  void sgtelib_server ( const std::string & model , const bool verbose );
  void sgtelib_predict ( const std::string & file_list , const std::string & model );
  void sgtelib_best    ( const std::string & file_list , const bool verbose );
  void sgtelib_help ( std::string word="GENERAL" );
  void sgtelib_test ( void );
}

#endif
