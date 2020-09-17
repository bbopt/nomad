/**
 \file   SgtelibModelFormulationType.cpp
 \brief  types for parameter SGTELIB_MODEL_FORMULATION (implementation)
 \author Viviane Rochon Montplaisir
 \date   December 2018
 \see    SgtelibModelFormulationType.hpp
 */

#include "../Type/SgtelibModelFormulationType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string (ex "FS", "FSP", "EIS"...)
// to a NOMAD::SgtelibModelFormulationType.
NOMAD::SgtelibModelFormulationType NOMAD::stringToSgtelibModelFormulationType(const std::string &sConst)
{
    auto ret = NOMAD::SgtelibModelFormulationType::UNDEFINED;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "FS")
    {
        ret = NOMAD::SgtelibModelFormulationType::FS;
    }
    else if (s == "FSP")
    {
        ret = NOMAD::SgtelibModelFormulationType::FSP;
    }
    else if (s == "EIS")
    {
        ret = NOMAD::SgtelibModelFormulationType::EIS;
    }
    else if (s == "EFI")
    {
        ret = NOMAD::SgtelibModelFormulationType::EFI;
    }
    else if (s == "EFIS")
    {
        ret = NOMAD::SgtelibModelFormulationType::EFIS;
    }
    else if (s == "EFIM")
    {
        ret = NOMAD::SgtelibModelFormulationType::EFIM;
    }
    else if (s == "EFIC")
    {
        ret = NOMAD::SgtelibModelFormulationType::EFIC;
    }
    else if (s == "PFI")
    {
        ret = NOMAD::SgtelibModelFormulationType::PFI;
    }
    else if (s == "D")
    {
        ret = NOMAD::SgtelibModelFormulationType::D;
    }
    else if (s == "EXTERN")
    {
        ret = NOMAD::SgtelibModelFormulationType::EXTERN;
    }
    else if (s == "UNDEFINED")
    {
        ret = NOMAD::SgtelibModelFormulationType::UNDEFINED;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::SgtelibModelFormulationType: " + s);
    }

    return ret;
}


std::string NOMAD::SgtelibModelFormulationTypeToString(const NOMAD::SgtelibModelFormulationType &smft)
{
    std::ostringstream oss;
    oss << smft;

    return oss.str();
}

