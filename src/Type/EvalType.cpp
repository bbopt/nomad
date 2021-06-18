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
 \file   EvalType.cpp
 \brief  types for Eval (implementation)
 \author Viviane Rochon Montplaisir
 \date   November 2019
 \see    EvalType.hpp
 */

#include "../Type/EvalType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string ("BB", "MODEL", "SURROGATE")
// to a NOMAD::EvalType.
// "UNDEFINED" or "LAST" throws an exception, as well as any value other than "BB", "MODEL".
NOMAD::EvalType NOMAD::stringToEvalType(const std::string &sConst)
{
    NOMAD::EvalType ret;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "BB")
    {
        ret = NOMAD::EvalType::BB;
    }
    else if (s == "MODEL")
    {
        ret = NOMAD::EvalType::MODEL;
    }
    else if (s == "SURROGATE")
    {
        ret = NOMAD::EvalType::SURROGATE;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::EvalType: " + s);
    }

    return ret;
}


// Convert a NOMAD::EvalType to a string.
// NOMAD::EvalType::UNDEFINED returns "UNDEFINED".
// NOMAD::EvalType::LAST throws an exception.
// An unrecognized eval type throws an exception.
std::string NOMAD::evalTypeToString(const NOMAD::EvalType& evalType)
{
    std::string s;

    switch(evalType)
    {
        case NOMAD::EvalType::BB:
            s = "BB";
            break;
        case NOMAD::EvalType::MODEL:
            s = "MODEL";
            break;
        case NOMAD::EvalType::SURROGATE:
            s = "SURROGATE";
            break;
        case NOMAD::EvalType::UNDEFINED:
            s = "UNDEFINED";
            break;
        case NOMAD::EvalType::LAST:
        default:
            throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized NOMAD::EvalType " + std::to_string((int)evalType));
            break;
    }

    return s;
}

