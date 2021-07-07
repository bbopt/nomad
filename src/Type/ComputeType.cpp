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
 \file   ComputeType.cpp
 \brief  Types for Computation of f and h (implementation)
 \author Viviane Rochon Montplaisir
 \date   February 2021
 \see    ComputeType.hpp
 */

#include "../Type/ComputeType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string (Ex. "STANDARD", "PHASE_ONE", "USER")
// to a NOMAD::ComputeType.
// "UNDEFINED" throws an exception, as well as any value other than "STANDARD", "PHASE_ONE", "USER".
NOMAD::ComputeType NOMAD::stringToComputeType(const std::string &sConst)
{
    NOMAD::ComputeType ret;
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "STANDARD")
    {
        ret = NOMAD::ComputeType::STANDARD;
    }
    else if (s == "PHASE_ONE")
    {
        ret = NOMAD::ComputeType::PHASE_ONE;
    }
    else if (s == "USER")
    {
        ret = NOMAD::ComputeType::USER;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::ComputeType: " + s);
    }

    return ret;
}


// Convert a NOMAD::ComputeType to a string.
// NOMAD::ComputeType::UNDEFINED returns "UNDEFINED".
// An unrecognized compute type returns an exception.
std::string NOMAD::computeTypeToString(const NOMAD::ComputeType& computeType)
{
    std::string s;

    switch(computeType)
    {
        case NOMAD::ComputeType::STANDARD:
            s = "STANDARD";
            break;
        case NOMAD::ComputeType::PHASE_ONE:
            s = "PHASE_ONE";
            break;
        case NOMAD::ComputeType::USER:
            s = "USER";
            break;
        case NOMAD::ComputeType::UNDEFINED:
            s = "UNDEFINED";
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized NOMAD::ComputeType " + std::to_string((int)computeType));
            break;
    }

    return s;
}

