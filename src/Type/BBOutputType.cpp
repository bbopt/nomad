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
 \file   BBOutputType.cpp
 \brief  types for BBOutput (implementation)
 \author Viviane Rochon Montplaisir
 \date   December 2018
 \see    BBOutputType.hpp
 */

#include "../Type/BBOutputType.hpp"
#include "../Util/ArrayOfString.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


// Convert a string (ex "OBJ", "EB", "PB"...)
// to a NOMAD::BBOutputType.
NOMAD::BBOutputType::BBOutputType(const std::string &sConst)
{
    _type = NOMAD::BBOutputType::Type::BBO_UNDEFINED;
    std::string s = sConst;
    NOMAD::toupper(s);

    // Detect if output is revealing (DiscoMads Algorithm)
    std::size_t pos = s.find("-R");      // position of R in str
    if (pos!=std::string::npos)
    {   
        // keep only substring defining output type    
        s = s.substr(0,pos);
        _isRevealing = true;
    }
    
    if (s == "OBJ")
    {
        _type = NOMAD::BBOutputType::Type::OBJ;
    }
    else if (s == "EB")
    {
        _type = NOMAD::BBOutputType::Type::EB;
    }
    else if (s == "PB" || s == "CSTR")
    {
        _type = NOMAD::BBOutputType::Type::PB;
    }
    else if (s == "RPB" )
    {
        _type = NOMAD::BBOutputType::Type::RPB;
    }
    else if (s == "CNT_EVAL")
    {
        _type = NOMAD::BBOutputType::Type::CNT_EVAL;
    }
    else if (s == "EXTRA_O" || s == "NOTHING" || s == "-" || s == "BBO_UNDEFINED")
    {
        _type = NOMAD::BBOutputType::Type::BBO_UNDEFINED;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for NOMAD::BBOutputType: " + s);
    }

}

// Convert a string containing multiple BBOutputTypes (ex "OBJ EB PB PB")
// to a NOMAD::BBOutputTypeList.
NOMAD::BBOutputTypeList NOMAD::stringToBBOutputTypeList(const std::string &s)
{
    NOMAD::BBOutputTypeList list;
    NOMAD::ArrayOfString aos(s);
    for (size_t i = 0; i < aos.size(); i++)
    {
        list.emplace_back(aos[i]);
    }
    return list;
}


// Convert a NOMAD::BBOutputTypeList into a string.
std::string NOMAD::BBOutputTypeListToString( const BBOutputTypeList & bbotList )
{
    std::ostringstream oss;
    for (const auto& bbot : bbotList )
    {
        oss << bbot;
        if (bbot._isRevealing)
        {
            oss << "-R";
        }
        oss << " ";
    }
    return oss.str();
}


// Count the number of constraints
size_t NOMAD::getNbConstraints(const BBOutputTypeList& bbotList)
{
    size_t nbConstraints = 0;
    for (const auto& bbot: bbotList)
    {
        if (bbot.isConstraint())
        {
            nbConstraints++;
        }
    }

    return nbConstraints;
}




// Count the number of objectives
size_t NOMAD::getNbObj(const BBOutputTypeList& bbotList)
{
    size_t nbObj = 0;
    for (const auto& bbot: bbotList)
    {
        if (bbot.NOMAD::BBOutputType::isObjective())
        {
            nbObj++;
        }
    }

    return nbObj;
}

// Count the number of revealing output
size_t NOMAD::getNbRevealing(const BBOutputTypeList& bbotList)
{
    size_t nbReveal = 0;
    for (const auto& bbot: bbotList)
    {
        if (bbot.NOMAD::BBOutputType::isRevealing())
        {
            nbReveal++;
        }
    }

    return nbReveal;
}

std::istream& NOMAD::operator>>(std::istream& is, NOMAD::BBOutputTypeList &bbOutputTypeList)
{
    std::string s;

    while (is >> s)
    {
        bbOutputTypeList.emplace_back(s);
    }

    return is;
}


std::string NOMAD::BBOutputType::display() const
{
    std::string s = "BBO_UNDEFINED";
    switch (_type)
    {
        case BBOutputType::Type::OBJ:
            s = "OBJ";
            break;
        case BBOutputType::Type::PB:
            s = "PB";
            break;
        case BBOutputType::Type::RPB:
            s = "RPB";
            break;
        case BBOutputType::Type::EB:
            s = "EB";
            break;
        case BBOutputType::Type::CNT_EVAL:
            s = "CNT_EVAL";
            break;
        case BBOutputType::Type::BBO_UNDEFINED:
        default:
            break;
    }
    
    if (_isRevealing)
    {
        s = s+"-R";
    }
    
    return s;
}
