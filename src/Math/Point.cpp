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
 \file   Point.cpp
 \brief  Custom class for points (implementation)
 \author Sebastien Le Digabel and Viviane Rochon Montplaisir
 \date   March 2017
 \see    Point.hpp
 */
#include "../Math/Point.hpp"

NOMAD::Point& NOMAD::Point::operator=(const NOMAD::Point &point)
{
    NOMAD::ArrayOfDouble::operator=(point);
    return *this;
}


NOMAD::Point& NOMAD::Point::operator=(const NOMAD::ArrayOfDouble &aod)
{
    NOMAD::ArrayOfDouble::operator=(aod);
    return *this;
}


/*---------*/
/* Display */
/*---------*/
// Regular display() has parenthesis around the coordinates. Ex:
// ( 3.46 6.85 5.72 5.85 )
// displayNoPar() ditches the parenthesis. Ex:
// 3.46 6.85 5.72 5.85
std::string NOMAD::Point::display(const NOMAD::ArrayOfDouble &prec) const
{
    return NOMAD::ArrayOfDouble::pStart + " " + NOMAD::ArrayOfDouble::display(prec) + " " + NOMAD::ArrayOfDouble::pEnd;
}


std::string NOMAD::Point::displayNoPar(const NOMAD::ArrayOfDouble &format) const
{
    return NOMAD::ArrayOfDouble::display(format);
}


/*------------------------------------------------------------------*/
/* Comparison function for use in the cache. Defines a weak order.  */
/* The additional condition for weak ordering is:                   */
/* For any X, Y, Z, if X < Y, then either X < Z or Y < Z.           */
/* See also: Comparison operator '<'.                               */
/* Note: if weakLess(lhs, rhs), then lhs < rhs.                    */
/* We don't use weakLess in operator() because it implies more     */
/* operations (truncations).                                        */
/*------------------------------------------------------------------*/
bool NOMAD::Point::weakLess(const NOMAD::Point &lhs, const NOMAD::Point &rhs)
{
    if (&lhs == &rhs)
    {
        return false;
    }

    if (lhs._n < rhs._n)
    {
        return true;
    }
    if (lhs._n > rhs._n)
    {
        return false;
    }

    for (size_t i = 0 ; i < lhs._n ; i++)
    {
        if (NOMAD::Double::weakLess(lhs[i], rhs[i]))
        {
            return true;
        }

        if (NOMAD::Double::weakLess(rhs[i], lhs[i]))
        {
            return false;
        }
    }

    return false;
}


/*----------------------------------------*/
/* Vector going from Point X to Point Y.  */
/*----------------------------------------*/
NOMAD::Direction NOMAD::Point::vectorize(const NOMAD::Point& X, const NOMAD::Point& Y)
{
    size_t n = X.size();
    if (n != Y.size())
    {
        throw NOMAD::Exception (__FILE__, __LINE__, "Cannot vectorize 2 points of different dimensions");
    }
    NOMAD::Direction Z(n);

    for (size_t i = 0; i < n; i++)
    {
        Z[i] = Y[i]-X[i];
    }
    return Z;
}


/* Point P resulting from adding Direction to this Point */
NOMAD::Point NOMAD::Point::operator+(const NOMAD::Direction& dir) const
{
    size_t n = size();
    if (n != dir.size())
    {
        throw NOMAD::Exception (__FILE__, __LINE__, "Cannot add a dimension to a point of different dimension");
    }

    NOMAD::Point P(n);
    for (size_t i = 0; i < n; i++)
    {
        P[i] = _array[i] + dir[i];
    }

    return P;
}


/*---------------------------------------*/
/* Distance between point X and point Y. */
/* Using norm L2.                        */
/*---------------------------------------*/
NOMAD::Double NOMAD::Point::dist(const NOMAD::Point& X, const NOMAD::Point& Y)
{
    NOMAD::Direction V = vectorize(X,Y);
    return V.norm();
}


NOMAD::Point NOMAD::Point::makeFullSpacePointFromFixed(const NOMAD::Point &fixedVariable) const
{
    size_t nbFixed = fixedVariable.nbDefined();
    size_t fullSpaceDim = fixedVariable.size();
    size_t subSpaceDim = fullSpaceDim - nbFixed;

    if (size() != subSpaceDim)
    {
        std::string s = "Error converting point " + this->NOMAD::Point::display();
        s += " (size " + std::to_string(this->size()) + ")";
        s += " to full space defined by fixed variable " + fixedVariable.NOMAD::Point::display();
        s += " (size " + std::to_string(fixedVariable.size()) + ")";
        s += ": point should be of size " + std::to_string(fixedVariable.size());
        s += " - " + std::to_string(nbFixed) + " = " + std::to_string(fixedVariable.size()-nbFixed);
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    NOMAD::Point fullSpacePoint = fixedVariable;

    if (0 == nbFixed)
    {
        // Fallback case.
        // no fixed variable defined - full space point is equal to this.
        fullSpacePoint = *this;
    }
    else
    {
        size_t iSub = 0;
        for (size_t i = 0; i < fullSpacePoint.size() && iSub < _n; i++)
        {
            if (!fullSpacePoint[i].isDefined())
            {
                fullSpacePoint[i] = _array[iSub];
                iSub++;
            }
        }
    }

    return fullSpacePoint;
}


NOMAD::Point NOMAD::Point::makeSubSpacePointFromFixed(const NOMAD::Point &fixedVariable, const bool verifyValues) const
{
    size_t nbFixed = fixedVariable.nbDefined();
    size_t fullSpaceDim = fixedVariable.size();
    size_t subSpaceDim = fullSpaceDim - nbFixed;

    if (size() != fullSpaceDim)
    {
        std::string s = "Error converting point " + this->NOMAD::Point::display();
        s += " (size " + std::to_string(this->size()) + ")";
        s += " to subspace defined by fixed variable " + fixedVariable.NOMAD::Point::display();
        s += " (size " + std::to_string(fixedVariable.size()) + ")";
        s += ": they should have the same size.";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    NOMAD::Point subSpacePoint(subSpaceDim);

    if (0 == nbFixed)
    {
        // Fallback case.
        // no fixed variable defined - sub space point is equal to this.
        subSpacePoint = *this;
    }
    else
    {
        size_t iSub = 0;
        for (size_t i = 0; i < fullSpaceDim && i < _n; i++)
        {
            if (i >= fixedVariable.size() || !fixedVariable[i].isDefined())
            {
                subSpacePoint[iSub] = _array[i];
                iSub++;
            }
            else if (verifyValues && fixedVariable[i].isDefined() && _array[i] != fixedVariable[i])
            {
                std::string s = "Error converting point " + this->display();
                s += " to subspace defined by fixed variable " + fixedVariable.display();
                throw NOMAD::Exception(__FILE__,__LINE__,s);
            }
        }
    }

    return subSpacePoint;
}


NOMAD::Point NOMAD::Point::projectPointToSubspace(const NOMAD::Point &fixedVariable) const
{
    return makeSubSpacePointFromFixed(fixedVariable, false);
}


bool NOMAD::Point::hasFixed(const NOMAD::Point &fixedVariable) const
{
    bool compatible = true;

    for (size_t i = 0; i < fixedVariable.size() && i < _n; i++)
    {
        if (fixedVariable[i].isDefined() && fixedVariable[i] != _array[i])
        {
            compatible = false;
            break;
        }
    }

    return compatible;
}


std::ostream& NOMAD::operator<< (std::ostream& out, const NOMAD::Point& point)
{
    out << point.display();
    return out;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::Point& point)
{
    int pointSize = 0;
    point.resize(pointSize);

    std::string s;
    is >> s;
    if (NOMAD::ArrayOfDouble::pStart != s)
    {
        is.setstate(std::ios::failbit);
        std::string err = "Expecting \"" + NOMAD::ArrayOfDouble::pStart + "\", got \"" + s + "\"";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    while (is >> s && NOMAD::ArrayOfDouble::pEnd != s)
    {
        pointSize++;
        point.resize(pointSize);
        point[pointSize-1].atof(s);
    }
    if (NOMAD::ArrayOfDouble::pEnd != s)
    {
        is.setstate(std::ios::failbit);
        std::string err = "Expecting \"" + NOMAD::ArrayOfDouble::pEnd + "\", got \"" + s + "\"";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }


    return is;

}








