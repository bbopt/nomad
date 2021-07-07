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
/*-------------------------------------------------------------------------------------*/
/**
 \file   Double.hpp
 \brief  Custom class for double-precision reals (headers)
 \author Sebastien Le Digabel
 \date   2010-04-02
 \see    Double.cpp
 */
#ifndef __NOMAD_4_0_DOUBLE__
#define __NOMAD_4_0_DOUBLE__

#include <cmath>

#include "../nomad_platform.hpp"
#include "../Util/defines.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"

#include "../nomad_nsbegin.hpp"

    /// Custom class for double-precision reals.
    /**
     - Allows comparisons on reals with custom precision.
     - Deals with undefined values.
     - Use \c todouble() to access the true double value.
     */
    class Double {

    private:
        double        _value;   ///< Value of the number.
        bool          _defined; ///< \c true if the number has a defined value.

        // \todo Make these local static objects
        DLL_UTIL_API static double      _epsilon;    ///< Desired precision on comparisons.
        DLL_UTIL_API static std::string _infStr;     ///< Infinity string.
        DLL_UTIL_API static std::string _undefStr;   ///< Undefined value string.

    public:

        /*-------------------------------------------------------------------*/

        /// Exception for undefined Double objects.
        class NotDefined : public Exception {
        public:
            /// Constructor.
            NotDefined ( const std::string & file ,
                         int                 line ,
                         const std::string & msg    )
            : Exception ( file , line , msg ) {}
        };

        /// Exception for divisions by zero with double objects.
        class InvalidValue : public Exception {
        public:
            /// Constructor.
            InvalidValue ( const std::string & file ,
                           int                 line ,
                           const std::string & msg    )
            : Exception ( file , line , msg ) {}
        };

        /*-------------------------------------------------------------------*/
#ifdef MEMORY_DEBUG
        /// Access to the number of double objects in memory.
        /**
         \return Number of double objects in memory.
         */
        static int getCardinality() { return Double::_cardinality; }

        /// Access to the max number of double objects in memory.
        /**
         \return Max number of double objects in memory.
         */
        static int getMaxCardinality() { return Double::_maxCardinality; }
#endif

        /// Constructor #1.
        explicit Double();

        /// Constructor #2.
        /**
         From a double.

         \param v The double to be copied into \c *this -- \b IN.
         */
        Double(const double & v);

        /// Copy constructor.
        /**
         \param d The double to be copied -- \b IN.
         */
        Double(const Double& d);

        /// Destructor.
        virtual ~Double();

        /// Conversion from a string to a double.
        /**
         \param s   The string to be converted -- \b IN.
         \return    \c true if the string is valid, \c false if not.
         */
        bool atof ( const std::string & s );

        /// Conversion from a string.
        /**
         * Value is determined by a string that may begin with \c 'r' to
         * indicate a proportion (relative value).
         *
         \param s      The string to be converted -- \b IN.
         \param rel    A flag indicating if the conversion is relative -- \b OUT.
         \return       \c true if the string is valid, \c false if not.
         */
        bool relativeAtof ( const std::string & s , bool & rel );

        /// Reset the double.
        void clear() { _value = 0.0; _defined = false; }

        /// Reset the double.
        void reset() { clear(); }

        /// Affectation operator #1.
        /**
         \param d   Right-hand side object -- \b IN.
         \return    Reference to \c *this as the result of the affectation.
         */
        Double & operator = ( const Double & d );

        /// Affectation operator #2.
        /**
         \param r   Right-hand side object -- \b IN.
         \return    Reference to \c *this as the result of the affectation.
         */
        Double & operator = ( double r );

        /// Access to the double value.
        const double & todouble() const;

        /// Get the double value, truncated with respect to epsilon.
        double trunk() const;
        
        /// Return the number of decimals of a double.
        std::size_t nbDecimals() const;

        /// Access to the double value.
        const std::string tostring() const;

        /// Is the value defined ?
        bool isDefined() const { return _defined; }

        /// Special way to set a double
        /**
         * Special trick: _defined is \c false, but _value is 1.0 (instead of default 0.0).
         * This means the value is to be set to some other value that we do not have access to immediately.
         * Normally, we should never access _value if _defined is \c false.
        */
         void setToBeDefined() { _defined = false; _value = 1.0; }


        /// Special way to assess if double is defined
        /**
         * Special trick: _defined is \c false, but _value is 1.0 (instead of default 0.0).
         * This means the value is to be set to some other value that we do not have
         * access to immediately.
         * Normally, we should never access _value if _defined is \c false.

         \return c true if \c *this is defined, \c false if not.
         */
        bool toBeDefined() const { return (!_defined && 1.0 == _value); }

        /// Is the value an integer ?
        bool isInteger() const;

        /// Is the value binary ?
        bool isBinary() const;

        /// Access to the precision.
        static double getEpsilon()  { return Double::_epsilon; }

        /// Set the precision.
        static void setEpsilon ( double eps );

        /// Access to the undefined value string.
        static std::string getUndefStr() { return Double::_undefStr; }

        /// Set the undefined value string.
        static void setUndefStr ( const std::string & undefStr )
        {
            Double::_undefStr = undefStr;
        }

        /// Access to the infinity string.
        static std::string getInfStr() { return Double::_infStr; }

        /// Set the infinity string.
        static void setInfStr ( const std::string & infStr )
        {
            Double::_infStr = infStr;
        }

        /// Rounding to the nearest integer.
        int round() const;


        /// Rounding to the nearest integer.
        const Double roundd() const;

        /**
        * Round the current value to given precision (number of decimals).
        \return \c true if rounding is done defined, \c false if not.
        */
        bool roundToPrecision(const NOMAD::Double & precision) ;
        
        /// Rounding upward to an integer.
        const Double ceil() const;

        /// Rounding downward to an integer.
        const Double floor() const;


        /// Absolute value.
        /**
         \return Max{\c -*this,\c *this}.
         */
        const Double abs() const;

        /// Square.
        /**
         \return \c *this \c * \c *this.
         */
        const Double pow2() const;

        /// Square root.
        /**
         \return \c (*this)^0.5.
         */
        const Double sqrt() const;

        /// Relative error with another double.
        /**
         * This error computation is based on:
         *
         * Ziv. Relative distance–an error measure in round-off error analysis.
         * Mathematics of Computation, 39(160):563–569, 1982. doi:10.1090/S0025-5718-1982-0669649-2.
         * The error will be in [0;2]

         \param  x  Value to compare -- b IN.
         \return    Relative error value in \c [0;1].
         */
        const Double relErr ( const Double & x ) const;

        /// Is the double a multiple of the granularity?
        /**
         * Both this value and granularity must be defined.
         * Granularity must be positive.
         * Otherwise, return \c false.

         \param  granularity    Granularity to compare
         \return                \c true if the value is a multiple of granularity.
         */
        bool isMultipleOf(const Double &granularity) const;

        /// Next multiple
        /**
         * Return the first double larger or equal to value
         * that is a multiple of the granularity.
         * Examples:
         * If granularity = 2.2 and value = 1.5, return 2.2.
         * If granularity = 2.2 and value = -1.5, return 0.0.
         * If granularity = 2.2 and value = 40, return 41.8.
         * If granularity = 2.2 and value = 39.6, return 39.6.

         \param granularity Granularity to compare
         \return            Next multiple Double
         */
        const Double nextMult(const Double &granularity) const;

        /// Previous multiple
        /**
          * Returns the first double lesser or equal to value
          * that is a multiple of the granularity.
          \param granularity Granularity to compare
          \return            Previous multiple Double
         */
        const Double previousMult(const Double &granularity) const;

        /// Operator \c ++ (prefix position).
        /**
         Allows \c ++d;

         \return    Reference to \c *this after incrementation by one.
         */
        Double & operator++ ();

        /// Operator \c ++ (suffix position).
        /**
         Allows \c d++;

         \param n   Increment.
         \return    Copy with value incremented by \c n.
         */
        Double operator++ ( int n );

        /// Operator \c -- (prefix position).
        /**
         Allows \c --d;

         \return    Reference to \c *this after decrementing by one.
         */    Double & operator-- ();

        /// Operator \c -- (suffix position).
        /**
         Allows \c d--;

         \return    Copy with value decremented by \c n.
         */
        Double operator-- ( int n );

        /// Operator \c +=.
        /**
         Allows \c d \c += \c d1.

         \param d1  Increment -- \b IN.
         \return    Reference to \c *this after decrement by \c d1.
         */
        const Double & operator += ( const Double & d1 );

        /// Operator \c -=.
        /**
         Allows \c d \c -= \c d1.

         \param d1  Decrement -- \b IN.
         \return    Reference to \c *this after decrement by \c d1.
         */
        const Double & operator -= ( const Double & d1 );

        /// Operator \c *=.
        /**
         Allows \c d \c *= \c d1.

         \param d1  Multiplicative factor -- \b IN.
         \return    Reference to \c *this after multiplication by \c d1.
         */
        const Double & operator *= ( const Double & d1 );

        /// Operator \c /=.
        /**
         Allows \c d \c /= \c d1. Throws a Exception::InvalidValue if \c d1==0.

         \param d1  Divisor -- \b IN.
         \return    Reference to \c *this after dividing by \c d1.
         */
        const Double & operator /= ( const Double & d1 );

        /// Weak comparison operator.
        /**
         * This propriety must be met:
         * If weakLess(d1, d2), then either weakLess(d1, d3), or weakLess(d3,d1).

         \param d1  First element of comparison
         \param d2  Second element of comparison
         \return    \c true if \c d1.trunk() < \c d2.trunk(), \c false if not.
         */
        static bool weakLess(const Double &d1, const Double &d2);

        /// Display.
        /**
        \param d   Object to be displayed -- \b IN.
        \param os  Reference to stream -- \b IN.
        \return    Reference to stream.
        */
        friend std::ostream& operator<< ( std::ostream& os, const Double& d );

        /// Display with a precision.
        /**
         \param prec      Display precision -- \b IN.
         \param refWidth  Width -- \b IN.
         \return          String.
         */
        std::string display(const int prec = DISPLAY_PRECISION_STD,
                            const size_t refWidth = 0) const;

        /// Display with a C style format
        /**
          * %f      w=-1 prec=-1 c='f'
          * %4.5f   w= 4 prec= 5 c='f'
          * %4f     w= 4 prec= 1 c='f'
          * %.5f    w=-1 prec= 5 c='f'
          * %.f     w=-1 prec= 0 c='f'
          *
          * c may be in 'e', 'E', 'f', 'g', 'G', 'd', or 'i'
          *
          * e Scientific notation (mantise/exponent) using e character 3.9265e+2
          * E Scientific notation (mantise/exponent) using E character 3.9265E+2
          * f Decimal floating point                                   392.65
          * g Use the shorter of %e or %f                              392.65
          * G Use the shorter of %E or %f                              392.65
          * d or i Integer rounded value                               393
          */
        std::string display(const std::string& format) const;

    };


    /*---------------------------------------------------------------------------*/


    /// Input.
    /**
     * - Allows the input of double objects with operator \c >>.
     * - Can read undefined values (parameter \c UNDEF_STR with default \c "-".)
     * - Example:
     * \code
     * Double d1 , d2;
     * std::cout << "Enter d1 and d2: ";
     * std::cin  >> d1 >> d2;
     * std::cout << "d1 and d2 are equal to " << d1 << " and " << d2 << std::endl;
     * \endcode

     \param in      A \c std::istream object (can be a file) -- \b IN/OUT.
     \param d       The double object to be read -- \b OUT.
     \return        The modified \c std::istream object.
     */
    std::istream & operator >> ( std::istream & in , Double & d );

    /// Inverse operator.
    /**
     Allows operations such as \c d \c = \c -d.

     \param d   Value to invert-- \b IN.
     \return    Inverted value.
     */
    inline const Double operator - ( const Double & d )
    {
        return Double (-d.todouble());
    }

    /// Addition operator \c +.
    /**
     Allows operations such as \c d \c = \c d1 \c + \c d2.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    Sum of the two elements.
     */
    inline const Double operator + ( const Double & d1 , const Double & d2 )
    {
        return Double ( d1.todouble() + d2.todouble() );
    }

    /// Substraction operator \c -.
    /**
     Allows operations such as \c d \c = \c d1 \c - \c d2.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    Difference between first and secont element.
     */
    inline const Double operator - ( const Double & d1 , const Double & d2 )
    {
        return Double (d1.todouble() - d2.todouble());
    }

    /// Multiplication operator \c *.
    /**
     Allows operations such as \c d \c = \c d1 \c * \c d2.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    Product of the two elements.
     */
    inline const Double operator * ( const Double & d1 , const Double & d2 )
    {
        return Double ( d1.todouble() * d2.todouble() );
    }

    /// Division operator \c /.
    /**
     Allows operations such as \c d \c = \c d1 \c / \c d2.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    Ratio of the first element by second element.
     */
    const Double operator / ( const Double & d1 , const Double & d2 );

    /// Comparison operator \c ==.
    /**
     Allows the comparison \c d1 \c == \c d2.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c true if \c d1 \c == \c d2, \c false otherwise.
     */
    inline bool operator == ( const Double & d1 , const Double & d2 )
    {
        return fabs ( d1.todouble() - d2.todouble() ) < Double::getEpsilon();
    }

    /// Comparison operator \c !=.
    /**
     Allows the comparison \c d1 \c != \c d2.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c true if \c d1 \c != \c d2, \c false otherwise.
     */
    inline bool operator != ( const Double & d1 , const Double & d2 )
    {
        return !(d1==d2);
    }

    /// Comparison operator \c <.
    /**
     Allows the comparison \c d1 \c < \c d2 accounting for \c Double precision.
     \note Defines a partial order, not a weak order.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c true if \c d1 \c < \c d2, \c false otherwise.
     */
    inline bool operator < ( const Double & d1 , const Double & d2 )
    {
        return d1.todouble() < d2.todouble() - Double::getEpsilon();
    }

    /// Comparison operator \c >.
    /**
     Allows the comparison \c d1 \c > \c d2 accounting for \c Double precision.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c true if \c d1 \c > \c d2, \c false otherwise.
     */
    inline bool operator > ( const Double & d1 , const Double & d2 )
    {
        return d1.todouble() > d2.todouble() + Double::getEpsilon();
    }

    /// Comparison operator \c <=.
    /**
     Allows the comparison \c d1 \c <= \c d2 accounting for \c Double precision.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c true if \c d1 \c <= \c d2, \c false otherwise.
     */
    inline bool operator <= ( const Double & d1 , const Double & d2 ) { return !(d1>d2); }

    /// Comparison operator \c >=.
    /**
     Allows the comparison \c d1 \c >= \c d2 accounting for \c Double precision.

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c true if \c d1 \c >= \c d2, \c false otherwise.
     */
    inline bool operator >= ( const Double & d1 , const Double & d2 )
    {
        return !(d1<d2);
    }

    /// Largest of two values \c >=.
    /**
     Return the largest of two \c Double, accounting for \c Double precision

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c max(d1,d2)
     */
    inline Double max (const Double d1 , const Double d2 ) { return (d1>d2)?d1:d2; }

    /// Smallest of two values \c >=.
    /**
     Return the smallest of two Double, accounting for \c Double precision

     \param d1  First element -- \b IN.
     \param d2  Second element -- \b IN.
     \return    \c min(d1,d2)
     */
    inline Double min ( const Double d1 , const Double d2 ) { return (d1<d2)?d1:d2; }




#include "../nomad_nsend.hpp"
#endif // __NOMAD_4_0_DOUBLE__
