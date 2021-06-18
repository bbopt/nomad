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
 \file   DirectionType.hpp
 \brief  Types for Poll direction : Ortho Mads (2n, n+1 Uni, n+1 neg), LT_MADS
 \author Christophe Tribes
 \date   May 2019
 \see    DirectionType.cpp
 */

#ifndef __NOMAD_4_0_DIRECTION_TYPE__
#define __NOMAD_4_0_DIRECTION_TYPE__

#include <list>
#include <sstream>
#include <vector>

#include "../nomad_nsbegin.hpp"

// Direction type
enum class DirectionType
{
    ORTHO_2N,
    ORTHO_NP1_NEG,
    ORTHO_NP1_QUAD,
    NP1_UNI,
    SINGLE,
    DOUBLE,
    LT_2N,
    LT_1,
    LT_2,
    LT_NP1,
    GPS_2N_STATIC,
    GPS_2N_RAND,
    GPS_BINARY,
    GPS_NP1_STATIC,
    GPS_NP1_STATIC_UNIFORM,
    GPS_NP1_RAND,
    GPS_NP1_RAND_UNIFORM,
    GPS_1_STATIC,
    UNDEFINED_DIRECTION
    ///< DirectionType is mandatory
};

typedef std::vector<DirectionType> DirectionTypeList;

/// Convert a list of strings (ex "ORTHO 2N", "ORTHO NP1") to a DirectionType.
DirectionType stringToDirectionType(const std::list<std::string> & ls);

/// Convert a string (ex "ORTHO 2N", "ORTHO NP1") to a DirectionType.
DirectionType stringToDirectionType(const std::string & s);

/// Convert an DirectionType to a string
std::string directionTypeToString(const DirectionType& dT);
/// Convert a DirectionTypeList to a string
std::string directionTypeListToString(const DirectionTypeList& dirTypeList);

inline std::ostream& operator<<(std::ostream& out, const DirectionType &directionType)
{
    out << directionTypeToString(directionType);
    return out;
}


inline std::ostream& operator<<(std::ostream& out, const DirectionTypeList &dirTypeList)
{
    out << directionTypeListToString(dirTypeList);
    return out;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD_4_0_DIRECTION_TYPE__
