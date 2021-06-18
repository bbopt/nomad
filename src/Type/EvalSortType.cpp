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
 \file   EvalSortType.cpp
 \brief  types for sorting EvalPoints (implementation)
 \author Viviane Rochon Montplaisir
 \date   April 2021
 \see    EvalSortType.hpp
 */

#include "../Type/EvalSortType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string ("LEXICOGRAPHICAL", "DIR_LAST_SUCCESS", "RANDOM", "SURROGATE")
// to a NOMAD::EvalSortType.
NOMAD::EvalSortType NOMAD::stringToEvalSortType(const std::string &sConst)
{
    NOMAD::EvalSortType ret;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "DIR_LAST_SUCCESS")
    {
        ret = NOMAD::EvalSortType::DIR_LAST_SUCCESS;
    }
    else if (s == "LEXICOGRAPHICAL")
    {
        ret = NOMAD::EvalSortType::LEXICOGRAPHICAL;
    }
    else if (s == "RANDOM")
    {
        ret = NOMAD::EvalSortType::RANDOM;
    }
    else if (s == "SURROGATE")
    {
        ret = NOMAD::EvalSortType::SURROGATE;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::EvalSortType: " + s);
    }

    return ret;
}


// Convert a NOMAD::EvalSortType to a string.
// An unrecognized EvalSortType returns an exception.
std::string NOMAD::evalSortTypeToString(const NOMAD::EvalSortType& evalSortType)
{
    std::string s;

    switch(evalSortType)
    {
        case NOMAD::EvalSortType::DIR_LAST_SUCCESS:
            s = "DIR_LAST_SUCCESS";
            break;
        case NOMAD::EvalSortType::LEXICOGRAPHICAL:
            s = "LEXICOGRAPHICAL";
            break;
        case NOMAD::EvalSortType::RANDOM:
            s = "RANDOM";
            break;
        case NOMAD::EvalSortType::SURROGATE:
            s = "SURROGATE";
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized NOMAD::EvalSortType " + std::to_string((int)evalSortType));
            break;
    }

    return s;
}

