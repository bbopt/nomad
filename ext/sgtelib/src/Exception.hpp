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

#ifndef __SGTELIB_EXCEPTION__
#define __SGTELIB_EXCEPTION__

#include <iostream>
#include <sstream>
#include <cstdlib>

namespace SGTELIB {

  class Exception : public std::exception {

  private:

    std::string _file;
    int         _line;
    std::string _err_msg;

    mutable std::string _tmp;

  public:

    virtual const char * what ( void ) const throw() {
      std::ostringstream oss;
      oss << _file << ":" << _line << " (" << _err_msg << ")";
      _tmp = oss.str();
      return _tmp.c_str();
    }

    Exception ( const std::string & file    ,
		int                 line    ,
		const std::string & err_msg   )
      : _file(file) , _line(line) , _err_msg(err_msg) {}
    
    virtual ~Exception ( void ) throw() {}
  };
}

// usage: throw SGTELIB::Exception ( "file.cpp" , __LINE__ , "error message" );

#endif
