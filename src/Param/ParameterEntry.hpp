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
  \file   ParameterEntry.hpp
  \brief  Parameter entry (headers)
  \author Sebastien Le Digabel
  \date   2010-04-05
  \see    ParameterEntry.cpp
*/
#ifndef __NOMAD_4_0_ParameterEntry__
#define __NOMAD_4_0_ParameterEntry__

#include "../Util/utils.hpp"

#include "../nomad_nsbegin.hpp"


/// Parameter entry.
/**
    - Describes the data relative to a parameter in a parameters file.
    - Objects of this class are stored in a ParameterEntries object.
*/
class ParameterEntry {

private:

    std::string            _name;   ///< Name of the parameter.
    std::list<std::string> _values; ///< List of values for the parameter.
    bool                   _ok;     ///< If the parameter is valid.
    bool                   _unique; ///< If the parameter is unique.
    std::shared_ptr<ParameterEntry> _next;   ///< Acces to the next parameter.

    std::string            _paramFile; ///< File from which this parameter was read
    int                    _line;       ///< Line for this parameter in _paramFile

    /// If the parameter has been interpreted.
    bool _hasBeenInterpreted;

public:

    /// Constructor.
    /**
       Ignores all entries after \c '#'.
       \param entry           A string describing the parameter entry -- \b IN.
       \param removeComments  A boolean equal to \c true if entries after
                              \c '#' are ignored -- \b IN
                              (Opt) (default = \c true).
    */
    explicit ParameterEntry(const std::string & entry, bool removeComments = true);

    /// Destructor.
    virtual ~ParameterEntry ( void ) {}

    /*---------------*/
    /*  GET methods  */
    /*---------------*/

    /// Access to the name of the parameter.
    /**
       \return The name.
    */
    const std::string & getName ( void ) const { return _name; }

    /// Access to the parameter values.
    /**
       \return The parameter values as a list of strings.
    */
    const std::list<std::string> & getValues ( void ) const { return _values; }

    std::string getAllValues ( void ) const;

    /// Access to the number of values of the parameter.
    /**
       \return The number of values.
    */
    size_t getNbValues ( void ) const { return _values.size(); }

    /// Access to the \c _ok flag.
    /**
       This flag is equal to \c true if the parameter entry is well defined.
       \return A boolean equal to \c true if the parameter is valid.
    */
    bool isOk ( void ) const { return _ok; }

    /// Access to the \c _unique flag.
    /**
       This flag is decided when a parameters file is read.
       \return A boolean equal to \c true if the parameter is unique
       in a parameters file.
    */
    bool isUnique ( void ) const { return _unique; }

    /// Access to another ParameterEntry.
    /**
       ParameterEntry objects are stored in a ParameterEntries
       object. The link between elements is assumed by the \c _next member
       returned by this function.
       \return A pointer to the next entry.
    */
    std::shared_ptr<ParameterEntry> getNext ( void ) const { return _next; }

    /// Access to the \c _hasBeenInterpreted flag.
    /**
       \return A boolean equal to \c true if the parameter has already
               been interpreted.
    */
    bool hasBeenInterpreted ( void ) const { return _hasBeenInterpreted; }

    /// Access to the parameter file of the parameter.
    /**
       \return The parameter file where this parameter was read.
    */
    const std::string & getParamFile ( void ) const { return _paramFile; }

    /// Access to the line number for this parameter in the parameter file.
    /**
       \return The line number at which this parameter can be found in the parameter file.
    */
    const int & getLine ( void ) const { return _line; }

    /*---------------*/
    /*  SET methods  */
    /*---------------*/

    /// Set the \c _next pointer.
    /**
       \param p A pointer to the next ParameterEntry to be inserted -- \b IN.
    */
    void setNext   ( std::shared_ptr<ParameterEntry> p ) { _next = p; }

    /// Set the \c _unique flag.
    /**
       \param u Value of the flag -- \b IN.
    */
    void setUnique ( bool u ) { _unique = u; }

    /// Set the \c _hasBeenInterpreted flag. to \c true.
    void setHasBeenInterpreted  ( void ) { _hasBeenInterpreted = true; }

    /// Set the name of the parameter file \c _paramFile
    void setParamFile(const std::string& paramFile ) { _paramFile = paramFile; }

    /// Set the line \c _line for this parameter in the parameter file
    void setLine( int line ) { _line = line; }

    /// Comparison with another entry.
    /**
       The comparison is based on the parameter name.

       \param p The right-hand side object -- \b IN.
       \return A boolean equal to \c true if \c this->_name \c < \c p._name.
    */
    bool operator < ( const ParameterEntry & p ) const { return _name < p._name; }

    /// Display.
    /**
       \param out The std::ostream object -- \b IN.
    */
    void display(std::ostream &out) const;
};

 /// Allows the comparison of two ParameterEntry objects.
 struct ParameterEntryComp {
    /// Comparison of two ParameterEntry objects.
    /**
       \param  p1 Pointer to the first ParameterEntry  -- \b IN.
       \param  p2 Pointer to the second ParameterEntry -- \b IN.
       \return A boolean equal to \c true if \c *p1 \c < \c *p2.
    */
    bool operator() ( const std::shared_ptr<ParameterEntry> p1 , const std::shared_ptr<ParameterEntry> p2 ) const
    {
        return (*p1 < *p2);
    }
};

/// Display a ParameterEntry object.
/**
     \param out The std::ostream object -- \b IN.
     \param e   The ParameterEntry object to be displayed -- \b IN.
     \return    The std::ostream object.
*/
inline std::ostream& operator<< (std::ostream &out,
                                   const ParameterEntry &e)
{
    e.display(out);
    return out;
}

#include "../nomad_nsend.hpp"


#endif  // __NOMAD_4_0_ParameterEntry__
