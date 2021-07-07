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
    NOMAD::trim(s);

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




