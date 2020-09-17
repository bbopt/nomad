/**
 \file   ListOfVariableGroup.cpp
 \brief  List of NOMAD::VariableGroup (implementation)
 \author Christophe Tribes
 \date   May 2020
 \see    ListOfVariableGroup.hpp
 */
#include "../Type/ListOfVariableGroup.hpp"

std::ostream& NOMAD::operator<< (std::ostream& out, const NOMAD::VariableGroup& vg)
{
    out << " ( " ;
    for (auto index : vg)
    {
        out << index;
    }
    out << " ) ";
    return out;
}

std::ostream& NOMAD::operator<< (std::ostream& out, const NOMAD::ListOfVariableGroup& lvg)
{
    size_t i=0;
    for (auto vg : lvg)
    {
        if (i > 0)
        {
            out << " ";
        }
        out << vg;
        i++;
    }
    return out;
}









