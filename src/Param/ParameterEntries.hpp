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
 \file   ParameterEntries.hpp
 \brief  Parameter entries (headers)
 \author Sebastien Le Digabel
 \date   2010-04-05
 \see    ParameterEntries.cpp
 */
#ifndef __NOMAD_4_0_PARAMETERENTRIES__
#define __NOMAD_4_0_PARAMETERENTRIES__

#include <set>
#include <vector>
#include "../Param/ParameterEntry.hpp"
#include "../Util/Uncopyable.hpp"

#include "../nomad_nsbegin.hpp"

/// Store and manage ParameterEntry.
/**
 - All the ParameterEntry objects are stored as a multiset with a comparison operator based on name.
 - ParameterEntries stores all the parameters provided in a file.
 */
class ParameterEntries : private Uncopyable {

private:

    /// List of ParameterEntry objects (the entries).
    std::multiset<std::shared_ptr<ParameterEntry>, ParameterEntryComp> _entries;

public:

    /// Constructor.
    explicit ParameterEntries ( void ) {}

    /// Destructor.
    virtual ~ParameterEntries ( void );

    /// Find a specific entry in a set.
    /**
     \param  name The name of the wanted ParameterEntry object -- \b IN.
     \return      A pointer to the ParameterEntry object if it
     has been found in the list of entries,
     or \c nullptr otherwise.
     */
    std::shared_ptr<ParameterEntry> find ( const std::string & name ) const;

    /// Insert a new entry in the list of entries.
    /**
     \param entry A pointer to the new ParameterEntry object -- \b IN.
     */
    void insert(std::shared_ptr<ParameterEntry> entry);

    /// Erase an entry from the list of entries.
    /**
     \param entry A pointer to the ParameterEntry object to erase -- \b IN.
     */
    void erase(std::shared_ptr<ParameterEntry> entry);

    /// Erase all entries.
    void eraseAll();

    /// Find a non-interpreted entry.
    /**
     \return A pointer to the first ParameterEntry that has not been
     interpreted so far,
     or \c nullptr if all entries have already been interpreted.
     */
    std::shared_ptr<ParameterEntry> findNonInterpreted ( void ) const;

    /// Find all non-interpreted entries
    std::vector<std::shared_ptr<ParameterEntry>> findAllNonInterpreted() const;

    /// Display.
    /**
     \param out The std::ostream object -- \b IN.
     */
    void display(std::ostream &out) const;
};

/// Display a ParameterEntries object.
/**
 \param out The std::ostream object -- \b IN.
 \param e   The ParameterEntries object to be displayed -- \b IN.
 \return    The std::ostream object.
 */
inline std::ostream& operator<< (std::ostream &out,
                                 const ParameterEntries &e)
{
    e.display(out);
    return out;
}

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_PARAMETERENTRIES__
