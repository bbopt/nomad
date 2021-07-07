/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/
/**
 \file   Exception.hpp
 \brief  Custom class for exceptions (headers)
 \author Sebastien Le Digabel
 \date   2010-03-29
 \see    Exception.cpp
 */
#ifndef __NOMAD_4_0_EXCEPTION__
#define __NOMAD_4_0_EXCEPTION__

#include <sstream>

#include "../nomad_nsbegin.hpp"

/// Custom class for exceptions.
/**
 * NOMAD uses this type of exceptions.
 * It indicates the file and line number at which a throw is made.
 *
 * \b Example
 *
 * \code
 * throw Exception(__FILE__, __LINE__, "an error message");
 * \endcode
 */
class Exception : public std::exception
{

private:

    mutable std::string _what;  ///< Error message.
    std::string         _file;  ///< File where the exception is thrown.
    size_t              _line;  ///< Line number at which the exception is thrown.

protected:
    std::string         _typeMsg;   ///< Basic exception message indicating the type of exception

public:

    /// Constructor.
    /**
     \param file A string corresponding to the file where the
     exception is thrown -- \b IN
     \param line An integer corresponding to the line number
     at which the exception is thrown -- \b IN.
     \param msg  A string corresponding to the error message -- \b IN.
     */
    Exception(const std::string& file, const size_t line, const std::string & msg)
      : _what(msg),
        _file(file),
        _line(line),
        _typeMsg("")
    {}

    /// Destructor.
    virtual ~Exception() throw() {}

    /// Access to the error message.
    /**
     \return A string with the error message.
     */
    const char * what() const throw();
};

class UserTerminateException : public Exception
{
public:
    // Constructor
    UserTerminateException(const std::string &file, const int line, const std::string &msg)
      : Exception(file, line, msg)
    {}
};
#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_EXCEPTION__
