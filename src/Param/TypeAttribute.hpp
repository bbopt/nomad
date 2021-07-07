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
//  TypeAttribute.hpp
//  nomad
//
//  Created by Christophe Tribes on 2017-12-08.
//  Copyright (c) 2017 GERAD. All rights reserved.
//

#ifndef __NOMAD_4_0_TYPEATTRIBUTE__
#define __NOMAD_4_0_TYPEATTRIBUTE__

#include "../Param/Attribute.hpp"
#include "../Util/defines.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for an Attribute of a templated type.
/**
 TypeAttribute is the templated class that derives from an Attribute and completes the values (current and initial) and type for a Nomad parameter. A specific type of attribute is obtained by calling the templated AttributeFactory::Create function.
 */
template <typename T>
class TypeAttribute : public Attribute {
private:
    T _value;
    const T _initValue;

public:

    TypeAttribute ( ) = delete;

    /// Constructor call by the AttributeFactory::Create
    template<typename ... ARGS>TypeAttribute(std::string Name,
                                             T initValue,
                                             bool algoCompatibilityCheck,
                                             bool restartAttribute,
                                             bool uniqueEntry,
                                             ARGS&& ... infoArgs):
    Attribute(Name,
              algoCompatibilityCheck,
              restartAttribute,
              uniqueEntry,
              std::forward<ARGS>(infoArgs)...),
    _value(initValue),
    _initValue(initValue)
    {}


    virtual ~TypeAttribute() {}

    const T& getValue() const { return _value; }
    const T& getInitValue() const { return _initValue; }
    void setValue(T v) { _value=v; }
    bool uniqueEntry() const { return _uniqueEntry; }

    /**
     \return \c true if current value equals initial value, \c false otherwise.
     */
    bool isDefaultValue() const { return (_value==_initValue) ; }

    void display( std::ostream & out , bool flagShortInfo = true ) const
    {

        out << _name << " " << _value;
        if ( flagShortInfo && _shortInfo.size() > 0 )
            out << " (" << _shortInfo << ")";
    }

    void resetToDefaultValue() noexcept { _value = _initValue;}

};


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_TYPEATTRIBUTE__
