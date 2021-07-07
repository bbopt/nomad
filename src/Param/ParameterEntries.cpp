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
 \file   ParameterEntries.cpp
 \brief  Parameter entries (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-05
 \see    ParameterEntries.hpp
 */
#include "../Param/ParameterEntries.hpp"

/*--------------------------------------------*/
/*                 destructor                 */
/*--------------------------------------------*/
NOMAD::ParameterEntries::~ParameterEntries ( void )
{
    // No need to delete anything since smart pointers are used.
}

/*--------------------------------------------*/
/*      finds a specific entry in the set     */
/*--------------------------------------------*/
std::shared_ptr<NOMAD::ParameterEntry> NOMAD::ParameterEntries::find ( const std::string & name ) const
{
    auto it = _entries.find ( std::make_shared<NOMAD::ParameterEntry>(name) );
    if ( it != _entries.end() )
        return (*it);
    return nullptr;
}

/*----------------------------------------*/
/*      inserts an entry into the set     */
/*----------------------------------------*/
void NOMAD::ParameterEntries::insert(std::shared_ptr<NOMAD::ParameterEntry> entry )
{
    std::shared_ptr<NOMAD::ParameterEntry> cur = find ( entry->getName() );
    if ( cur )
    {
        entry->setUnique ( false );
        cur->setUnique   ( false );
        while ( cur->getNext() )
            cur = cur->getNext();
        cur->setNext ( entry );
    }
    _entries.insert ( entry );
}


/*----------------------------------------*/
/*       erase an entry from the set      */
/*----------------------------------------*/
void NOMAD::ParameterEntries::erase(std::shared_ptr<NOMAD::ParameterEntry> entry )
{
    _entries.erase ( entry );
}

/*----------------------------------------*/
/*       erase all entries                */
/*----------------------------------------*/
void NOMAD::ParameterEntries::eraseAll ( )

{
    _entries.clear();
}

/*----------------------------------------*/
/*       find a non-interpreted entry     */
/*----------------------------------------*/
std::shared_ptr<NOMAD::ParameterEntry> NOMAD::ParameterEntries::findNonInterpreted ( void ) const
{
    for (auto it : _entries)
    {
        if ( !it->hasBeenInterpreted() )
        {
            return it;
        }
    }
    return nullptr;
}


std::vector<std::shared_ptr<NOMAD::ParameterEntry>> NOMAD::ParameterEntries::findAllNonInterpreted() const
{
    std::vector<std::shared_ptr<NOMAD::ParameterEntry>> allNonInterp;
    for (auto it : _entries)
    {
        if ( !it->hasBeenInterpreted() )
        {
            allNonInterp.push_back(it);
        }
    }
    return allNonInterp;
}


/*--------------------------------------------*/
/*                   display                  */
/*--------------------------------------------*/
void NOMAD::ParameterEntries::display(std::ostream &out) const
{
    auto end = _entries.end();
    for (auto it = _entries.begin() ; it != end ; ++it )
    {
        out << **it << std::endl;
    }
}
