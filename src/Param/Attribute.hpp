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
//
//  Attribute.hpp
//  nomad
//
//  Created by Christophe Tribes on 2017-12-07.
//  Copyright (c) 2017 GERAD. All rights reserved.
//

#ifndef __NOMAD_4_0_ABSTRACTATTRIBUTE__
#define __NOMAD_4_0_ABSTRACTATTRIBUTE__

#include "../Util/defines.hpp"

#include "../nomad_nsbegin.hpp"

/// An attribute contains all meta data about a Nomad parameter except its type,  default and current value.
/**
 TypeAttribute is the templated class that derives from an Attribute and completes the values (current and initial) and type for a Nomad parameter. A specific type of attribute is obtained by calling the templated AttributeFactory::Create function.\n

 A Nomad parameter has a name, some help and info, some  keywords and some additional flags :
 - Attribute::_uniqueEntry, a flag for a parameter that cannot appear more than once in a file.
 - Attribute::_algoCompatibilityCheck, a flag for a parameter that can be used to check if two sets of parameters can produce the same run.
 - Attribute::_restartAttribute, a flag for a parameter that can be changed when doing a restart.
 */
class Attribute {
public:

    virtual const std::string & getName(){ return _name; }
    virtual const std::string & getShortInfo(){ return _shortInfo; }
    virtual const std::string & getHelpInfo(){ return _helpInfo; }
    virtual const std::string & getKeywords(){ return _keywords; }
    virtual bool isForAlgoCompatibilityCheck() { return _algoCompatibilityCheck; }
    virtual bool isRestartAttribute() { return _restartAttribute; }
    virtual bool getParamFromUniqueEntry() { return _uniqueEntry; }
    virtual bool isInternal() { return _internal; }

    void setShortInfo(const std::string& s) { _shortInfo = s; }
    void setHelpInfo(const std::string& s) { _helpInfo = s; }
    void setKeywords(const std::string& s) { _keywords = s; }

    bool hasEmptyInfo () const
    {
        if ( _helpInfo.empty() && _shortInfo.empty() )
            return true;
        else
            return false;
    }

    virtual ~Attribute() {}

    virtual void resetToDefaultValue() = 0;

    Attribute (const std::string& Name, bool algoCompatibilityCheck,
               bool restartAttribute, bool uniqueEntry,
               const std::string& ShortInfo,
               const std::string& HelpInfo, const std::string& Keywords)
        : _name(Name),
          _shortInfo(ShortInfo),
          _helpInfo(HelpInfo),
          _keywords(Keywords),
          _algoCompatibilityCheck(algoCompatibilityCheck),
          _restartAttribute(restartAttribute),
          _uniqueEntry(uniqueEntry),
          _internal(false)
    { if (Keywords.find("internal") != std::string::npos)
        _internal=true ;
    }

    Attribute (const std::string& Name, bool algoCompatibilityCheck,
               bool restartAttribute, bool uniqueEntry,
               const std::string& ShortInfo, const std::string& HelpInfo)
        : _name(Name),
          _shortInfo(ShortInfo),
          _helpInfo(HelpInfo),
          _algoCompatibilityCheck(algoCompatibilityCheck),
          _restartAttribute(restartAttribute),
          _uniqueEntry(uniqueEntry),
          _internal(false)
    {}

    Attribute (const std::string& Name, bool algoCompatibilityCheck,
               bool restartAttribute, bool uniqueEntry, const std::string& ShortInfo)
        : _name(Name),
          _shortInfo(ShortInfo),
          _algoCompatibilityCheck(algoCompatibilityCheck),
          _restartAttribute(restartAttribute),
          _uniqueEntry(uniqueEntry),
          _internal(false)
    {}

    Attribute(const std::string& Name)
        : _name(Name),
        _algoCompatibilityCheck(false),
        _restartAttribute(false),
        _uniqueEntry(true),
        _internal(false)
    {}

    virtual void display( std::ostream& os , bool flagShortInfo = true ) const
    {
        os << _name << " " ;
        if ( flagShortInfo && _shortInfo.size() > 0 )
        {
            os << " (" << _shortInfo << ")";
        }
    }

protected:

    std::string _name; ///< The name of a parameter is used to get access to the value.
    std::string _shortInfo;
    std::string _helpInfo;
    std::string _keywords;  ///< Registered keywords
    bool        _algoCompatibilityCheck; ///< Flag for parameter that can be used to check if two sets of parameters can produce the same run.
    bool        _restartAttribute; ///< Flag for parameter that can be changed when doing a restart
    bool        _uniqueEntry; ///< Flag for a parameter that cannot appear more than once in a file
    bool        _internal; ///< Flag for interanl parameter
};


inline std::ostream & operator << ( std::ostream & os, const Attribute & att)
{
    att.display(os);
    return os;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD_4_0_ABSTRACTATTRIBUTE__
