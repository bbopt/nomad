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

#ifndef __KERNEL__
#define __KERNEL__

#include "Defines.hpp"
#include "Exception.hpp"
#include "Matrix.hpp"
namespace SGTELIB {

  // kernel type:
  enum kernel_t {
    KERNEL_D1 ,
    KERNEL_D2 ,
    KERNEL_D3 ,
    KERNEL_D4 ,
    KERNEL_D5 ,
    KERNEL_D6 ,
    KERNEL_D7 ,
    KERNEL_I0 ,
    KERNEL_I1 ,
    KERNEL_I2 ,
    KERNEL_I3 ,
    KERNEL_I4
  };
  const int NB_KERNEL_TYPES = 11;
  const int NB_DECREASING_KERNEL_TYPES = 6;

  // kernel
  double kernel ( const SGTELIB::kernel_t kt , const double ks ,const double r );
  SGTELIB::Matrix kernel ( const SGTELIB::kernel_t kt , const double ks , SGTELIB::Matrix R );
  // kernel is decreasing ?
  bool kernel_is_decreasing ( const SGTELIB::kernel_t kt );
  // kernel has a shape parameter ?
  bool kernel_has_parameter ( const SGTELIB::kernel_t kt );
  // kernel has a shape parameter ?
  int kernel_dmin ( const SGTELIB::kernel_t kt );
  // Kernel to str
  std::string kernel_type_to_str  ( SGTELIB::kernel_t       );
  // string to kernel type
  SGTELIB::kernel_t str_to_kernel_type ( const std::string & s );
  SGTELIB::kernel_t int_to_kernel_type ( const int i );
}

#endif

