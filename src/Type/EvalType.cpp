/**
 \file   EvalType.cpp
 \brief  types for Eval (implementation)
 \author Viviane Rochon Montplaisir
 \date   November 2019
 \see    EvalType.hpp
 */

#include "../Type/EvalType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string ("BB", "SGTE")
// to a NOMAD::EvalType.
// "UNDEFINED" throws an exception, as well as any value other than "BB" or "SGTE".
NOMAD::EvalType NOMAD::stringToEvalType(const std::string &sConst)
{
    NOMAD::EvalType ret;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "BB")
    {
        ret = NOMAD::EvalType::BB;
    }
    else if (s == "SGTE")
    {
        ret = NOMAD::EvalType::SGTE;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::EvalType: " + s);
    }

    return ret;
}


// Convert a NOMAD::EvalType to a string.
// NOMAD::EvalType::UNDEFINED returns "UNDEFINED".
// An unrecognized eval type returns an exception.
std::string NOMAD::evalTypeToString(const NOMAD::EvalType& evalType)
{
    std::string s;

    switch(evalType)
    {
        case NOMAD::EvalType::BB:
            s = "BB";
            break;
        case NOMAD::EvalType::SGTE:
            s = "SGTE";
            break;
        case NOMAD::EvalType::UNDEFINED:
            s = "UNDEFINED";
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized NOMAD::EvalType " + std::to_string((int)evalType));
            break;
    }

    return s;
}

