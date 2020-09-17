/**
 \file   LHSearchType.hpp
 \brief  Type for parameter LH_SEARCH
 \author Viviane Rochon Montplaisir
 \date   January 2019
 \see    LHSearchType.cpp
 */


#ifndef __NOMAD400_LH_SEARCH_TYPE__
#define __NOMAD400_LH_SEARCH_TYPE__

#include <string>
#include <sstream>

#include "../nomad_nsbegin.hpp"

/// Class for handling the types of LH search.
/**
  *  An initial LH search can be used instead of providing X0. \n
  *  The secondary LH searches can be done at each iteration if no
    other searches has obtained a success.
 */
class LHSearchType
{
private:
    bool    _enable;    ///< Flag for calling LH search.
    size_t  _lhsearch0; ///< Number of points for initial LH search.
    size_t  _lhsearch1; ///< Number of points for secondary LH searches.

public:
    /// Constructor
    /**
     \param entries The parameters for initial search and secondary searches given as a \c string -- \b IN.
     */
    LHSearchType(const std::string& entries = "") ;


    /// Test if LH search is enabled
    /**
     \return \c true if LH search is enabled \c false if not.
     */
    bool isEnabled() const { return _enable; }

    /// Number of initial LH search points
    /**
     \return The number of points for initial LH search.
     */
    size_t getNbInitial() const { return _lhsearch0; }

    /// Number of secondary LH search points at each iteration
    /**
     \return The number of points for secondary LH searches.
     */
    size_t getNbIteration() const { return _lhsearch1; }

    /// Comparison operator \c ==.
    /**
     \param lhst    The right-hand side object -- \b IN.
     \return        A boolean equal to \c true if  \c *this \c == \c array.
     */
    bool operator== (const LHSearchType & lhst) const
    {
        return (( _enable == lhst._enable) && ( _lhsearch0 == lhst._lhsearch0 ) && ( _lhsearch1 == lhst._lhsearch1 )) ;
    }

    /// Comparison operator \c !=.
    /**
     \param lhst    The right-hand side object -- \b IN.
     \return        A boolean equal to \c true if  \c *this \c != \c array.
     */
    bool operator!= (const LHSearchType &lhst ) const { return !(*this == lhst); }


};


inline std::ostream& operator<<(std::ostream& os, const LHSearchType &lhsearch)
{
    os << lhsearch.getNbInitial() << " " << lhsearch.getNbIteration();

    return os;
}


#include "../nomad_nsend.hpp"

#endif  // __NOMAD400_LH_SEARCH_TYPE__
