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
 \file   BBInputType.hpp
 \brief  Types for BBInput
 \author Viviane Rochon Montplaisir
 \date   December 2018
 \see    BBInputType.cpp
 */


#ifndef __NOMAD_4_0_BB_INPUT_TYPE__
#define __NOMAD_4_0_BB_INPUT_TYPE__

#include <sstream>
#include <vector>

#include "../nomad_nsbegin.hpp"

/// Enum for blackbox input type
/**
 Related to problem formulation
 \note Categorical variables not supported yet.
*/
enum class BBInputType
{
    CONTINUOUS  ,     ///< Continuous variable (default) (R)
    ALL_CONTINUOUS  , ///< All variables are continuous variable (default, *R). Need a checkAndComply to set the BBInputTypeList.
    INTEGER     ,     ///< Integer variable (I)
    ALL_INTEGER ,     ///< All variables are integer (*I). Need a checkAndComply to set the BBInputTypeList.
    //CATEGORICAL ,   ///< Categorical variable          (C)
    BINARY         ,  ///< Binary variable               (B)
    ALL_BINARY       ///< All variables are binary (*B). Need a checkAndComply to set the BBInputTypeList.
};


typedef std::vector<BBInputType> BBInputTypeList;
typedef std::vector<BBInputType>::const_iterator BBInputTypeListIt;


/// Utility for BBInputTypes.
/**
 Convert a string ("R", "*R", "I", "*I", "B", "*B") to a blackbox input type.
 */
BBInputType stringToBBInputType(const std::string &s);

/// Utility for BBInputTypes
/**
 * Convert a string containing multiple BBInputTypes (ex "( R I B R )") to a BBInputTypeList or a single 'all of the same type (*)' BBInputType to BBInputTypeList of a single element. When have a single 'all of the same type (*)', the list of full dimension is created when calling PbParameters::checkAndComply. \n
 *
 *\todo Support the syntax 1-4 I
 *
 */
BBInputTypeList stringToBBInputTypeList(const std::string &s);


inline std::ostream& operator<<(std::ostream& out, const BBInputType &bbinputtype)
{
    switch(bbinputtype)
    {
        case BBInputType::CONTINUOUS:
            out << "R";
            break;
        case BBInputType::INTEGER:
            out << "I";
            break;
        case BBInputType::BINARY:
            out << "B";
            break;
        default:
            out << "R"; // Unrecognized, output default: CONTINUOUS
            break;
    }
    return out;
}

/// Display a list of blackbox input type.
inline std::ostream& operator<<(std::ostream& out, const BBInputTypeList &bbinputtypelist)
{
    BBInputTypeListIt it;
    bool first = true;
    for (it = bbinputtypelist.begin(); it != bbinputtypelist.end(); ++it)
    {
        if (!first)
        {
            out << " ";
        }
        out << *it;
        first = false;
    }
    return out;
}

#include "../nomad_nsend.hpp"
#endif  // __NOMAD_4_0_BB_INPUT_TYPE__
