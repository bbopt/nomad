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
 \file   ArrayOfDouble.hpp
 \brief  Base class to be used for an array of n \c Double values, as base class for \c Point and \c Direction.
 \author Viviane Rochon Montplaisir
 \date   October 2017
 \see    ArrayOfDouble.cpp
 */

#ifndef __NOMAD_4_0_ARRAYOFDOUBLE__
#define __NOMAD_4_0_ARRAYOFDOUBLE__

#include <numeric>
#include "../Math/Double.hpp"
#include "../Util/ArrayOfString.hpp"

#include "../nomad_nsbegin.hpp"

/// \brief Class for the representation of an array of n values.
/**
 An array of n values is defined by its size and its coordinates.
*/
class ArrayOfDouble {

public:

    static const std::string pStart; ///< Static variable used for array delimitation.
    static const std::string pEnd; ///< Static variable used for array delimitation.

protected:
    /*---------*/
    /* Members */
    /*---------*/

    size_t _n;          ///< Dimension of the array.
    Double* _array;     ///< Values of the array.

public:
    /*-------------*/
    /* Constructor */
    /*-------------*/

    /// Constructor
    /**
     \param n   Dimension of the array -- \b IN (Opt) (default = 0).
     \param val Initial value for all elements of the array -- \b IN (Opt) (default = undefined real).
     */
    explicit ArrayOfDouble(const size_t n = 0, const Double &val = Double());
    
    /// Constructor #2 -- used by PyNomad interface
    /**
     \param v vector of double -- \b IN .
     */
    explicit ArrayOfDouble(const std::vector<double> & v);
    

    /// Copy constructor.
    /**
     \param coords Array object to be copied -- \b IN.
     */
    ArrayOfDouble(const ArrayOfDouble &coords);

    /// Affectation operator.
    /**
     \param coords  Right-hand side object -- \b IN.
     \return        Reference to \c *this as the result of the affectation.
     */
    ArrayOfDouble& operator= (const ArrayOfDouble &coords);

    /// Destructor.
    virtual ~ArrayOfDouble();

    /*---------*/
    /* Get/Set */
    /*---------*/
    /// Operator \c [].
    /**
     \param i   Index (0 for the first element) -- \b IN.
     \return    \c (i+1)th coordinate of the array.
     */
    Double& operator[](size_t i) const;

    /// Access to the dimension of the array.
    /**
     \return Dimension of the array.
     */
    size_t size() const { return _n; }

    /***********************/
    /* Other class methods */
    /***********************/

    /// Change the \c ArrayOfDouble dimension, and set all coordinates to d.
    /**
     \param n   New dimension -- \b IN (Opt) (default = 0).
     \param d   Initial value for all coordinates -- \b IN (Opt) (default = undefined Double).
     */
    void reset(size_t n = 0, const Double &d = Double());

    /// Change the dimensionof the array. The values are kept.
    /**
     \param n New dimension of the array -- \b IN.
     \param d Initial value for new coordinates, if enlarging the array.
     */
    void resize(size_t n, const Double &d = Double());

    /// Test if empty
    /**
       \return A \c bool equal to \c true if the array is empty, /c false if not.
    */
    bool isEmpty() const { return _n == 0; }

    /// Set all the coordinates to a specific value.
    /**
       \param d The value for all coordinates -- \b IN.
    */
    void set(const Double &d) const
    { std::fill(_array, _array + _n , d); }

    /// Set a single coordinate to a specific value. If the flag relative is
    /// true, use the bounds to set the value (0 <= d <= 1).
    /**
     \param d           The value for all coordinates -- \b IN.
     \param index       The index of the coordinate -- \b IN.
     \param relative    The flag to set the coordinate relative to the bounds -- \b IN (Opt).
     \param lb          The lower bound for relative value -- \b IN (Opt) (only if relative=true)
     \param ub          The upper bound for relative value -- \b IN (Opt) (only if relative=true)
     */
    void set(size_t index,
             const Double &d,
             bool relative = false,
             const Double &lb = Double(),
             const Double &ub = Double());


    /// Set the coordinates with an array of reals.
    /**
       \param n Dimension of the array -- \b IN.
       \param a Array of size \c n of reals -- \b IN.
    */
    void set(size_t n, const Double* a);

    /// Check if all the coordinates are defined.
    /**
       \return A \c bool equal to \c true if all the coordinates are defined, \c false if not.
    */
    bool isComplete() const;

    /// Check if at least one coordinate is defined.
    /**
       \return A \c bool equal to \c true if at least one coordinate is defined, \c false if not.
    */
    virtual bool isDefined() const;

    /// Count the number of defined values.
    /**
       \return The number of values that are defined.
    */
    size_t nbDefined() const;

    // Return A \c bool equal to \c true if at least one value to be defined, \c false if not.
    bool toBeDefined() const;

    /// Snap an array to the bounds. Ignore mesh.
    void snapToBounds(const ArrayOfDouble &lb,
                      const ArrayOfDouble &ub);

    /// Verify if the array is inside the bounds. Ignores undefined bounds.
    bool inBounds(const ArrayOfDouble &lowerBound,
                  const ArrayOfDouble &upperBound) const;

    /// Read values and fill the array with corresponding double values.
    void readValuesAsArray(const ArrayOfString& valueString);

    /// Return an ArrayOfDouble which holds the absolute value of
    /// all defined values
    ArrayOfDouble abs() const;

    /// Return max of all defined values.
    Double max() const;

    /// Mutiplication with a scalar.
    /**
     - This implements \c *this \c = \c d \c * \c *this.
     - The current object \c *this is modified.
     \param d   The scalar -- \b IN.
     \return    The AOD times \c d.
     */
    const ArrayOfDouble & operator *= ( const Double & d );

    /// Division by a scalar.
    /**
     - This implements \c *this \c = \c d \c / \c *this.
     - The current object \c *this is modified.
     \param d   The scalar -- \b IN.
     \return    The AOD times \c d.
     */
    const ArrayOfDouble & operator /= ( const Double & d );

    /// Addition with another array.
    /**
     The current object \c *this is not modified.
     \param p   The other AOD -- \b IN.
     \return    A third AOD equal to \c *this \c + \c p.
     */
    const ArrayOfDouble operator+(const ArrayOfDouble& p) const;

    /// Substraction with an other array.
    /**
     The current object \c *this is not modified.
     \param p   The other AOD -- \b IN.
     \return    A third AOD equal to \c *this \c - \c p.
     */
    const ArrayOfDouble operator-(const ArrayOfDouble& p) const;

    /**
     * Verify that all elements of this array are multiples of granularity,
     * for non-zero values.
     * If that is the case, index is set to -1.
     * If that is not the case, index is set to the first index violating the granularity.
     \param gran    The granularity to test the array -- \b IN.
     \param index   Index of array violating condition -- \b OUT.
     \return        A \c bool equal to \c true if array is multiple of granularity, \c false if not.
     */
    bool isMultipleOf(const ArrayOfDouble &gran, int &index) const;

    /**
    * Round all elements to their precision given as number of decimals.
    \param precision    The number of decimals for each element -- \b IN.
    \return        A \c bool equal to \c true if a single element has been rounded.
    */
    bool roundToPrecision(const NOMAD::ArrayOfDouble & precision);
    
    
    /*------------*/
    /* Comparison */
    /*------------*/
    /// Comparison operator \c ==.
    /**
     \param coords  The right-hand side object -- \b IN.
     \return        A \c bool equal to \c true if  \c *this \c == \c coords. If any value is not defined, return \c false.
     */
    bool operator== (const ArrayOfDouble &coords) const;

    /// Comparison operator \c !=.
    /**
     \param coords The right-hand side object -- \b IN.
     \return       A \c bool equal to \c true if  \c *this \c != \c coords.
     */
    bool operator!= (const ArrayOfDouble &coords) const { return !(*this == coords); }

    /// Comparison operator \c <=.
    /**
     * Throws an exception if the sizes do not match.
     \param coords  The right-hand side object -- \b IN.
     \return        A \c bool equal to \c true if each coordinate of \c *this is inferior or equal to the value of \c coords at the same index.
     */
    virtual bool operator<= (const ArrayOfDouble &coords) const;

    /// Comparison operator \c <.
    /**
     * Throws an exception if the sizes do not match.
     \param coords  The right-hand side object -- \b IN.
     \return        A \c bool equal to \c true if each coordinate of \c *this is inferior or equal to the value of \c coords at the same index, and at least one coordinate is strictly inferior.
     */
    virtual bool operator< (const ArrayOfDouble &coords) const;

    /// Display with a given precision
    virtual std::string display(const ArrayOfDouble &prec = ArrayOfDouble()) const;

protected:
    //

    /// Helper function to verify that n1 == n2
    void verifySizesMatch(size_t n1, size_t n2, std::string filename, size_t linenum) const;

    /// Helper function to compare arrays
    /**
     \param aod             Array for comparison -- \b IN.
     \param isInf           Result of comparison -- \b OUT.
     \param isStrictInf     Result of comparison -- \b OUT.
     */
    void compare(const ArrayOfDouble& aod,
                 bool &isInf,
                 bool &isStrictInf) const;

};


/// Display.
/**
 \param aod   Object to be displayed -- \b IN.
 \param out   An output stream -- \b IN.
 \return      Reference to the modified output stream
 */
std::ostream& operator<<(std::ostream& out, const ArrayOfDouble& aod);

/// Input.
/**
 * - Allows the input of \c ArrayOfDouble objects with operator \c >>.
 \param in      A \c std::istream object (can be a file) -- \b IN/OUT.
 \param aod     The \c ArrayOfDouble object to be read -- \b OUT.
 \return        The modified \c std::istream object.
 */
std::istream& operator>>(std::istream& in, ArrayOfDouble& aod);

#include "../nomad_nsend.hpp"
#endif // __NOMAD_4_0_ARRAYOFDOUBLE__
