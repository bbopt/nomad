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

