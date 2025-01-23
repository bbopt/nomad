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
 \file   BBOutputType.hpp
 \brief  types for BBOutput
 \author Viviane Rochon Montplaisir
 \date   September 2018
 \see    BBOutput.hpp
 */
#ifndef __NOMAD_4_5_BB_OUTPUT_TYPE__
#define __NOMAD_4_5_BB_OUTPUT_TYPE__

#include <string>
#include <sstream>
#include <vector>

#include "../nomad_platform.hpp"
#include "../nomad_nsbegin.hpp"

struct DLL_UTIL_API BBOutputType
{
    
    /// Blackbox output types
    enum Type
    {
        OBJ,        ///< Objective value
        EB,         ///< Extreme barrier constraint
        PB,         ///< Progressive barrier constraint
        RPB,         ///<  Revealed constraint (obtained from DiscoMads revealing mechanism, not from blackbox evaluation). Similar to BBO_Undefined but it is considered in H calculation.
        CNT_EVAL,   ///< Output set to 0 or 1 to count the blackbox evaluation or not
        BBO_UNDEFINED ///< Output ignored
    };

    Type _type;
    
    // flag for a Revealing bbo.
    bool _isRevealing = false;
    
    // Implicit constructor
    BBOutputType()
    {
        _type = NOMAD::BBOutputType::Type::BBO_UNDEFINED;
    }
    
    /**
     Convert a string (ex "OBJ", "EB", "PB"...)
     to a BBOutputType.
     */
    BBOutputType(const std::string & s) ;
    
    /**
     Make a BBOutputType from its type.
     */
    BBOutputType(const BBOutputType::Type & t, bool isRevealing = false)
    {
        _type = t;
        _isRevealing = isRevealing;
    }
    
    // Assignment operator
    BBOutputType(const BBOutputType & b)
    {
        _type = b._type;
        _isRevealing = b._isRevealing;
    }
    
    /// Verify if the BBOutputType defines a constraint
    bool isConstraint() const
    {
        return (_type == NOMAD::BBOutputType::Type::EB || _type == NOMAD::BBOutputType::Type::PB || _type == NOMAD::BBOutputType::Type::RPB);
    }
    
    
    bool isObjective() const
    {
        return (_type == NOMAD::BBOutputType::Type::OBJ);
    }

    bool isExtraOutput() const
    {
        return (_type == NOMAD::BBOutputType::Type::BBO_UNDEFINED);
    }
    
    bool isRevealing() const
    {
        return(_isRevealing);
    }
    
    void setRevealingStatus(bool isRevealing)
    {
        _isRevealing = isRevealing;
    }
    
    std::string display() const;
    
    
    bool operator==(const BBOutputType& bbo) const
    {
        return (this->_type == bbo._type);

    };
    bool operator!=(const BBOutputType& bbo) const
    {
        return this->_type != bbo._type;

    };
    

    bool operator==(const BBOutputType::Type& type) const
    {
        return (this->_type == type);

    };
    bool operator!=(const BBOutputType::Type& type) const
    {
        return this->_type != type;

    };
    
};



/// Definition for the list of blackbox output types
typedef std::vector<BBOutputType> BBOutputTypeList;

typedef BBOutputTypeList::const_iterator BBOutputTypeListIt;


/// Utility for BBOutputType
/**
 Convert a string containing multiple BBOutputTypes (ex "OBJ EB PB PB")
 to a BBOutputTypeList.
 */
DLL_UTIL_API BBOutputTypeList stringToBBOutputTypeList(const std::string &s);

/// Utility for BBOutputType
/**
 Convert a BBOutputTypeList into a string
 */
DLL_UTIL_API std::string BBOutputTypeListToString ( const BBOutputTypeList & bbotList );

///// Helper to test if a BBOutputType is a constraint (PB, EB, ....)
DLL_UTIL_API bool BBOutputTypeIsConstraint(const BBOutputType & bbotType);

/// Count the number of constraints
DLL_UTIL_API size_t getNbConstraints(const BBOutputTypeList& bbotList);

/// Count the number of objectives
DLL_UTIL_API size_t getNbObj(const BBOutputTypeList& bbotList);

/// Count the number of revealing output (for discomads)
DLL_UTIL_API size_t getNbRevealing(const BBOutputTypeList& bbotList);

/// Read and interpret BBOutputType
inline std::ostream& operator<<(std::ostream& os, const BBOutputType &bbot)
{
    os << bbot.display();

    return os;
}

/// Display BBOutputType
inline std::ostream& operator<<(std::ostream& out, const BBOutputTypeList &bboutputtypelist)
{
    BBOutputTypeListIt it;
    bool first = true;
    for (it = bboutputtypelist.begin(); it != bboutputtypelist.end(); ++it)
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


DLL_UTIL_API std::istream& operator>>(std::istream& is, BBOutputTypeList& bbOutputTypeList);


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_5_BB_OUTPUT_TYPE__
