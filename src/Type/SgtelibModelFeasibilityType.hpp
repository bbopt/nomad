/**
 \file   SgtelibModelFeasibilityType.hpp
 \brief  types for parameter SGTELIB_MODEL_FEASIBILITY
 \author Viviane Rochon Montplaisir
 \date   July 2019
 \see    SgtelibModel.hpp
 */
#ifndef __NOMAD400_SGTELIB_MODEL_FEASIBILITY_TYPE__
#define __NOMAD400_SGTELIB_MODEL_FEASIBILITY_TYPE__

#include <string>
#include <sstream>

#include "../nomad_nsbegin.hpp"

// Feasibility types for sgtelib model Search
enum class SgtelibModelFeasibilityType
{
    C         , /// one model for each constraint
    H         , /// one model for H
    B         , /// one binary model
    M         , /// one model of the max of (c_j)
    UNDEFINED   /// Undefined
};


// Convert a string (ex "C", "H", "B"...)
// to a SgtelibModelFeasibilityType.
SgtelibModelFeasibilityType stringToSgtelibModelFeasibilityType(const std::string &s);

std::string SgtelibModelFeasibilityTypeToString(const SgtelibModelFeasibilityType &smft);

inline std::ostream& operator<<(std::ostream& os, const SgtelibModelFeasibilityType &smft)
{
    switch (smft)
    {
        case SgtelibModelFeasibilityType::C:
            os << "C";
            break;
        case SgtelibModelFeasibilityType::H:
            os << "H";
            break;
        case SgtelibModelFeasibilityType::B:
            os << "B";
            break;
        case SgtelibModelFeasibilityType::M:
            os << "M";
            break;
        default:
            return os << "UNDEFINED";
            break;
    }

    return os;
}



#include "../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIB_MODEL_FEASIBILITY_TYPE__
