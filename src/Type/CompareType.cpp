/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
 \file   CompareType.cpp
 \brief  Types for Comparison of objective vector f and h (implementation)
 \author Ludovic Salomon
 \date   February 2022
 \see    CompareType.hpp
 */

#include "../Type/CompareType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string (Ex. "EQUAL", "INDIFFERENT", "DOMINATED", "DOMINATING")
// to a NOMAD::CompareType.
// "UNDEFINED" throws an exception, as well as any value other than "EQUAL", "INDIFFERENT",
// "DOMINATED", "DOMINATING".
NOMAD::CompareType NOMAD::stringToCompareType(const std::string &sConst)
{
    NOMAD::CompareType ret;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "EQUAL")
    {
        ret = NOMAD::CompareType::EQUAL;
    }
    else if (s == "INDIFFERENT")
    {
        ret = NOMAD::CompareType::INDIFFERENT;
    }
    else if (s == "DOMINATED")
    {
        ret = NOMAD::CompareType::DOMINATED;
    }
    else if (s == "DOMINATING")
    {
        ret = NOMAD::CompareType::DOMINATING;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::CompareType: " + s);
    }

    return ret;
}


// Convert a NOMAD::CompareType to a string.
// NOMAD::CompareType::UNDEFINED returns "UNDEFINED".
// An unrecognized compute type returns an exception.
std::string NOMAD::compareTypeToString(NOMAD::CompareType compareType)
{
    std::string s;

    switch(compareType)
    {
        case NOMAD::CompareType::EQUAL:
            s = "EQUAL";
            break;
        case NOMAD::CompareType::INDIFFERENT:
            s = "INDIFFERENT";
            break;
        case NOMAD::CompareType::DOMINATED:
            s = "DOMINATED";
            break;
        case NOMAD::CompareType::DOMINATING:
            s = "DOMINATING";
            break;
        case NOMAD::CompareType::UNDEFINED:
            s = "UNDEFINED";
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized NOMAD::CompareType " + std::to_string((int)compareType));
            break;
    }

    return s;
}

