/**
 \file   BBInputType.cpp
 \brief  types for BBInput (implementation)
 \author Viviane Rochon Montplaisir
 \date   December 2018
 \see    BBInputType.hpp
 */

#include "../Type/BBInputType.hpp"
#include "../Util/ArrayOfString.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"

// Convert a string ("R", "I", "B") to a NOMAD::BBInputType.
NOMAD::BBInputType NOMAD::stringToBBInputType(const std::string &sConst)
{
    NOMAD::BBInputType ret = NOMAD::BBInputType::CONTINUOUS;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "R")
    {
        ret = NOMAD::BBInputType::CONTINUOUS;
    }
    else if (s == "*R")
    {
        ret = NOMAD::BBInputType::ALL_CONTINUOUS;
    }
    else if (s == "I")
    {
        ret = NOMAD::BBInputType::INTEGER;
    }
    else if (s == "*I")
    {
        ret = NOMAD::BBInputType::ALL_INTEGER;
    }
    else if (s == "B")
    {
        ret = NOMAD::BBInputType::BINARY;
    }
    else if (s == "*B")
    {
        ret = NOMAD::BBInputType::ALL_BINARY;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::BBInputType: " + s);
    }

    return ret;
}


// Convert a string containing multiple BBInputTypes (ex "( R I B R ) or *I") to a NOMAD::BBInputTypeList.
// Supporting both classic version with parenthesis and modern version without parenthesis.
NOMAD::BBInputTypeList NOMAD::stringToBBInputTypeList(const std::string &s)
{
    NOMAD::BBInputTypeList bbInputType;
    NOMAD::ArrayOfString aos(s);
    std::size_t arraysize = aos.size();
    if (arraysize >= 2 && aos[0] == "(" && aos[arraysize-1] == ")")
    {
        // * is not supported inside the vector notation
        if (s.find("*") < std::string::npos)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::BBInputType: " + s);
        }

        aos.erase(arraysize-1);
        aos.erase(0);
        arraysize -= 2;

        for (size_t i = 0; i < arraysize; i++)
        {
            bbInputType.push_back(NOMAD::stringToBBInputType(aos[i]));
        }
    }

    // Manage the 'all of the same type (*)' situation
    if (s.find("*") < std::string::npos)
    {
        // Concatenate all strings together before interpretation to BBInputType
        std::string ss;
        for (size_t i = 0; i < arraysize; i++)
        {
            ss+=aos[i];
        }

        bbInputType.push_back(NOMAD::stringToBBInputType(ss));
    }

    if (arraysize > 0 && bbInputType.size() == 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::BBInputType: " + s);
    }

    return bbInputType;
}

