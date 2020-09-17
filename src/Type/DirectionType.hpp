/**
 \file   DirectionType.hpp
 \brief  Types for Poll direction : Ortho Mads (2n, n+1 Uni, n+1 neg), LT_MADS
 \author Christophe Tribes
 \date   May 2019
 \see    DirectionType.cpp
 */

#ifndef __NOMAD400_DIRECTION_TYPE__
#define __NOMAD400_DIRECTION_TYPE__

#include <list>
#include <sstream>

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


/// Convert a list of strings (ex "ORTHO 2N", "ORTHO NP1")  to a DirectionType.
DirectionType stringToDirectionType(const std::list<std::string> & ls);

/// Convert a string (ex "ORTHO 2N", "ORTHO NP1")  to a DirectionType.
DirectionType stringToDirectionType(const std::string & s);

/// Convert an EvalType to a string
std::string directionTypeToString (const DirectionType& dT);

inline std::ostream& operator<<(std::ostream& out, const DirectionType &directionType)
{
    out << directionTypeToString(directionType);
    return out;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD400_DIRECTION_TYPE__
