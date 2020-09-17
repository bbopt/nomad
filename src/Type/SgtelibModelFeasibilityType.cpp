/**
 \file   SgtelibModelFeasibilityType.cpp
 \brief  types for parameter SGTELIB_MODEL_FEASIBILITY (implementation)
 \author Viviane Rochon Montplaisir
 \date   December 2018
 \see    SgtelibModelFeasibilityType.hpp
 */

#include "../Type/SgtelibModelFeasibilityType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string (ex "C", "H", "B"...)
// to a NOMAD::SgtelibModelFeasibilityType.
NOMAD::SgtelibModelFeasibilityType NOMAD::stringToSgtelibModelFeasibilityType(const std::string &sConst)
{
    auto ret = NOMAD::SgtelibModelFeasibilityType::UNDEFINED;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "C")
    {
        ret = NOMAD::SgtelibModelFeasibilityType::C;
    }
    else if (s == "H")
    {
        ret = NOMAD::SgtelibModelFeasibilityType::H;
    }
    else if (s == "B")
    {
        ret = NOMAD::SgtelibModelFeasibilityType::B;
    }
    else if (s == "M")
    {
        ret = NOMAD::SgtelibModelFeasibilityType::M;
    }
    else if (s == "UNDEFINED")
    {
        ret = NOMAD::SgtelibModelFeasibilityType::UNDEFINED;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::SgtelibModelFeasibilityType: " + s);
    }

    return ret;
}


std::string NOMAD::SgtelibModelFeasibilityTypeToString(const NOMAD::SgtelibModelFeasibilityType &smft)
{
    std::ostringstream oss;
    oss << smft;

    return oss.str();
}




