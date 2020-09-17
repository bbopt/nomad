/**
 \file   ListOfVariableGroup.hpp
 \brief  List of VariableGroup
 \author Christophe Tribes
 \date   May 2020
 \see    ListOfVariableGroup.cpp
 */

#ifndef __NOMAD400_LIST_OF_VARIABLE_GROUP__
#define __NOMAD400_LIST_OF_VARIABLE_GROUP__

#include <set>
#include <sstream>
#include <list>

#include "../nomad_nsbegin.hpp"

/// Type definition for the representation of a variable group (set of indices).
typedef std::set<size_t> VariableGroup;

std::ostream& operator<<(std::ostream& out, const VariableGroup& vg);



/// Type definition for the representation of a vector of points.
typedef std::list<VariableGroup> ListOfVariableGroup;

std::ostream& operator<<(std::ostream& out, const ListOfVariableGroup& lvg);

// Not implemented:
//std::istream& operator>>(std::istream& in, ListOfVariableGroup& aop);



#include "../nomad_nsend.hpp"
#endif // __NOMAD400_LIST_OF_VARIABLE_GROUP__
