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
 \file   Point.hpp
 \brief  Point
 \author Viviane Rochon Montplaisir
 \date   March 2017
 \see    Point.cpp
 */

#ifndef __NOMAD_4_0_POINT__
#define __NOMAD_4_0_POINT__

#include <numeric>
#include "../Math/ArrayOfDouble.hpp"
#include "../Math/Direction.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for the representation of a point.
/**
 A point is defined by its size and its coordinates from ArrayOfDouble. Addition and comparison (< and weak less) operators are provided.
 A point can be converted from a sub space to a full space and the reverse (Point::makeFullSpacePointFromFixed and Point::makeSubSpacePointFromFixed).
*/
class Point : public ArrayOfDouble
{
public:
    /*--------------*/
    /* Constructors */
    /*--------------*/
    /**
     \param n   Dimension of the point -- \b IN.
     \param val Initial value for all coordinates -- \b IN.
     */
    explicit Point(const size_t n = 0, const Double& val = Double())
      : ArrayOfDouble(n, val)
    {}

    /**
     \param v Initial value for all coordinates given as a vector of double - \b IN.
     */
    explicit Point(const std::vector<double> & v)
      : ArrayOfDouble(v)
    {}
    
    /// Copy constructors.
    /**
     \param pt The point to copy -- \b IN.
     */
    Point(const Point &pt)
      : ArrayOfDouble(pt)
    {}

    /// Assignment operator
    /**
     \param pt The point to assign -- \b IN.
     */
    Point& operator=(const Point& pt);

    /// Assignment operator
    /**
     \param aod The array of double to assign -- \b IN.
     */
    Point& operator=(const ArrayOfDouble& aod);

    /// Destructor.
    virtual ~Point() {}

    /*---------------*/
    /* Class methods */
    /*---------------*/

    /// Formated point display.
    /**
     * Put parenthesis around the coordinates. Ex: ( 3.46 6.85 5.72 5.85 )
     *
     \param prec Precision for display -- \b IN
     \return     Formated string.
     */
    std::string display(const ArrayOfDouble &prec = ArrayOfDouble()) const override;
    //void display(std::ostream& out) const override;

    /// Formated point display.
    /**
     * No parenthesis around the coordinates. Ex: 3.46 6.85 5.72 5.85
     *
     \param prec Precision for display -- \b IN.
     \return     Formated string.
     */
    std::string displayNoPar(const ArrayOfDouble &prec = ArrayOfDouble()) const;

    //void displayNoPar(std::ostream& out) const;

    /*------------*/
    /* Comparison */
    /*------------*/
    /// Comparison operator \c <.
    /**
     * Ends up to lexicographic order. Not useful. Throw an exception.
     * Ignore warning about unused variable point
     */
    bool operator<(const Point &NOMAD_UNUSED(point)) const
    {
        throw Exception(__FILE__,__LINE__,"Error: Attempting to use lexicographical order on a Point.");
    }

    /// Weak comparison operator.
    /**
     * This propriety must be met:
     * For a given d3, if weakLess(d1, d2), then either weakLess(d1, d3), or weakLess(d3,d1).
     *
     \param lhs  Left element of comparison
     \param rhs  Right element of comparison
     \return     \c bool for comparison.
     */
    static bool weakLess(const Point &lhs, const Point &rhs);

    /// Addition point = point + direction
    /**
     The current object \c *this is not modified.

     \param dir The direction to add -- \b IN.
     \return    A \c point equal to \c *this \c + \c dir.
     */
    Point operator+(const Direction& dir) const;

    /*----------*/
    /* Distance */
    /*----------*/
    /// Euclidean distance between two points.
    /**
     \param X   First point -- \b IN.
     \param Y   Second point -- \b IN.
     \return    The distance.
     */
    static Double dist(const Point& X, const Point& Y);

    /*--------*/
    /* Vector */
    /*--------*/

    /// Create a direction by substracting 2 points.
    /**
     \param X   Right point -- \b IN.
     \param Y   Left point -- \b IN.
     \return    A direction = Y-X.
     */
    static Direction vectorize(const Point& X, const Point& Y);


     /// Convert a point from sub space to full space using fixed variables.
    /**
     \param fixedVariable   Fixed values given as a point,
     \return                Full space \c Point.
     */
    Point makeFullSpacePointFromFixed(const Point &fixedVariable) const;


    /// Convert a point from full space to sub space using fixed variables.
    /**
     \param fixedVariable   Fixed values given as a point,
     \param verifyValues    If true, the Point's values must be the same as the ones defined by fixedVariable.
     \return                Sub space point.
     \see projectPointToSubspace
     */
    Point makeSubSpacePointFromFixed(const Point &fixedVariable, const bool verifyValues = true) const;


    /// Project a point from full space to sub space using fixed variables.
    /**
     \param fixedVariable   Fixed values given as a point,
     \return                Sub space point.
     \note the Point's values are not verified.
     \see makeSubSpacePointFromFixed
     */
    Point projectPointToSubspace(const Point &fixedVariable) const;


    /// Verify if a point is part of the sub-space defined by fixed variable
    /**
     \param fixedVariable   Fixed values given as a point,
     \return                \c true if \c *this is part of sub-space
     */
    bool hasFixed(const Point &fixedVariable) const;

};

/// Display.
/**
 \param pt   Object to be displayed -- \b IN.
 \param out  An output stream -- \b IN.
 \return     Reference to the modified output stream
 */
std::ostream& operator<<(std::ostream& out, const Point& pt);


/// Input.
/**
 * - Allows the input of point objects with operator \c >>.
 \param in      A \c std::istream object (can be a file) -- \b IN/OUT.
 \param pt      A point object to be read -- \b OUT.
 \return        The modified \c std::istream object.
 */
std::istream& operator>>(std::istream& in, Point& pt);



#include "../nomad_nsend.hpp"
#endif // __NOMAD_4_0_POINT__
