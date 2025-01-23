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
#include "DMultiMadsSearchStrategyType.hpp"

#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"

// Convert a string ("DOM", "MULTI") to a NOMAD::DMultiMadsNMSearchType.
NOMAD::DMultiMadsNMSearchType NOMAD::stringToDMultiMadsNMSearchType(const std::string& sConst)
{
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "DOM")
    {
        return NOMAD::DMultiMadsNMSearchType::DOM;
    }
    else if (s == "MULTI")
    {
        return NOMAD::DMultiMadsNMSearchType::MULTI;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::DMultiMadsNMSearchType: " + s);
    }

    return NOMAD::DMultiMadsNMSearchType::DOM;
}

// Convert a NOMAD::DMultiMadsNMSearchType to a string
std::string NOMAD::DMultiMadsNMSearchTypeToString(NOMAD::DMultiMadsNMSearchType NMSearchStrategy)
{
    if (NMSearchStrategy == DMultiMadsNMSearchType::DOM)
    {
        return "DOM";
    }
    else if (NMSearchStrategy == DMultiMadsNMSearchType::MULTI)
    {
        return "MULTI";
    }
    else
    {
        return "UNDEFINED";
    }
}

// Convert a string ("DOM", "MULTI") to a NOMAD::DMultiMadsQuadSearchType.
NOMAD::DMultiMadsQuadSearchType NOMAD::stringToDMultiMadsQuadSearchType(const std::string& sConst)
{
    std::string s = sConst;
    NOMAD::toupper(s);

    if (s == "DMS")
    {
        return NOMAD::DMultiMadsQuadSearchType::DMS;
    }
    else if (s == "DOM")
    {
        return NOMAD::DMultiMadsQuadSearchType::DOM;
    }
    else if (s == "MULTI")
    {
        return NOMAD::DMultiMadsQuadSearchType::MULTI;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::DMultiMadsQuadSearchType: " + s);
    }

    return NOMAD::DMultiMadsQuadSearchType::DMS;
}

// Convert a NOMAD::DMultiMadsNMSearchType to a string
std::string NOMAD::DMultiMadsQuadSearchTypeToString(NOMAD::DMultiMadsQuadSearchType quadSearchStrategy)
{
    if (quadSearchStrategy == DMultiMadsQuadSearchType::DMS)
    {
        return "DMS";
    }
    else if (quadSearchStrategy == DMultiMadsQuadSearchType::DOM)
    {
        return "DOM";
    }
    else if (quadSearchStrategy == DMultiMadsQuadSearchType::MULTI)
    {
        return "MULTI";
    }
    else
    {
        return "UNDEFINED";
    }
}
