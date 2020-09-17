/**
 \file   ArrayOfPoint.cpp
 \brief  Vector of NOMAD::Points (implementation)
 \author Viviane Rochon Montplaisir
 \date   August 2019
 \see    ArrayOfPoint.hpp
 */
#include "../Math/ArrayOfPoint.hpp"



std::ostream& NOMAD::operator<< (std::ostream& out, const NOMAD::ArrayOfPoint& aop)
{
    for (size_t i = 0; i < aop.size(); i++)
    {
        if (i > 0)
        {
            out << " ";
        }
        out << aop[i].display();
    }
    return out;
}









