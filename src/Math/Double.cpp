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
 \file   Double.cpp
 \brief  Custom class for double-precision reals (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-02
 \see    Double.hpp
 */
#include <cctype>   // for toupper
#include <iomanip>  // For std::setprecision
#include "../Math/Double.hpp"
#include "../Util/defines.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
double      NOMAD::Double::_epsilon         = NOMAD::DEFAULT_EPSILON;
std::string NOMAD::Double::_infStr          = NOMAD::DEFAULT_INF_STR;
std::string NOMAD::Double::_undefStr        = NOMAD::DEFAULT_UNDEF_STR;
#ifdef MEMORY_DEBUG
int         NOMAD::Double::_cardinality     = 0;
int         NOMAD::Double::_maxCardinality  = 0;
#endif

/*-----------------------------------------------*/
/*                  Constructor 1                */
/*-----------------------------------------------*/
NOMAD::Double::Double()
  : _value(0.0),
    _defined(false)
{
#ifdef MEMORY_DEBUG
    ++NOMAD::Double::_cardinality;
    if ( NOMAD::Double::_cardinality > NOMAD::Double::_maxCardinality )
        ++NOMAD::Double::_maxCardinality;
#endif
}

/*-----------------------------------------------*/
/*                  Constructor 2                */
/*-----------------------------------------------*/
NOMAD::Double::Double(const double & v)
  : _value(v),
    _defined(true)
{
#ifdef MEMORY_DEBUG
    ++NOMAD::Double::_cardinality;
    if (NOMAD::Double::_cardinality > NOMAD::Double::_maxCardinality)
        ++NOMAD::Double::_maxCardinality;
#endif
}

/*-----------------------------------------------*/
/*                  Copy Constructor             */
/*-----------------------------------------------*/
NOMAD::Double::Double(const NOMAD::Double &d)
  : _value(d._value),
    _defined(d._defined)
{
#ifdef MEMORY_DEBUG
    ++NOMAD::Double::_cardinality;
    if (NOMAD::Double::_cardinality > NOMAD::Double::_maxCardinality)
        ++NOMAD::Double::_maxCardinality;
#endif
}

/*-----------------------------------------------*/
/*                    destructor                 */
/*-----------------------------------------------*/
NOMAD::Double::~Double()
{
#ifdef MEMORY_DEBUG
    --NOMAD::Double::_cardinality;
#endif
}

/*-----------------------------------------------*/
/*               set epsilon (static)            */
/*-----------------------------------------------*/
void NOMAD::Double::setEpsilon (double eps)
{
    if (eps <= 0.0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "NOMAD::Double::setEpsilon(): invalid epsilon");
    }
    else if (eps < std::numeric_limits<double>::epsilon())
    {
        std::ostringstream oss;
        oss << "NOMAD::Double::setEpsilon(): minimum value for epsilon is std epsilon = ";
        oss << std::numeric_limits<double>::epsilon();
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }
    NOMAD::Double::_epsilon = eps;
}


/*-----------------------------------------------*/
/*            get the value as double            */
/*-----------------------------------------------*/
const double & NOMAD::Double::todouble() const
{
    if (! _defined)
    {
        throw NotDefined(__FILE__, __LINE__,
                          "NOMAD::Double::todouble(): value not defined");
    }
    return _value;
}


/*-----------------------------------------------------*/
/*            get the truncated double value           */
/*             relative to current epsilon             */
/*-----------------------------------------------------*/
double NOMAD::Double::trunk() const
{
    if (! _defined)
    {
        throw NotDefined(__FILE__, __LINE__,
                          "NOMAD::Double::trunk(): value not defined");
    }

    double trunk =  _epsilon * std::floor(_value / _epsilon);
    return trunk;

}

bool NOMAD::Double::roundToPrecision(const NOMAD::Double & precision )
{
    if (! _defined)
    {
        throw NotDefined(__FILE__, __LINE__,
                          "NOMAD::Double::roundToPrecision(): value not defined");
    }
    
    bool modif = false;
    if (precision.isDefined())
    {
        if (precision > 0)
        {
            double powprec = std::pow(10,precision.round());
            _value = std::round(_value * powprec) / powprec;
        }
        else
        {
            // Integer, binary (and categorical) have no decimal
            _value = std::round(_value);
        }
        modif = true;
    }
    
    return modif;

}

bool NOMAD::Double::weakLess(const NOMAD::Double &d1, const NOMAD::Double &d2)
{
    return (d1.trunk() < d2.trunk());
}



/*-----------------------------------------------*/
/*            get the value as string            */
/*-----------------------------------------------*/
const std::string NOMAD::Double::tostring() const
{
    return display(NOMAD::DISPLAY_PRECISION_STD);
}


/*------------------------------------------*/
/*                    input                 */
/*------------------------------------------*/
std::istream & NOMAD::operator>> ( std::istream & in , NOMAD::Double & d )
{
    std::string s;
    in >> s;

    if ( !in.fail() && !d.atof (s) )
    {
        in.setstate ( std::ios::failbit );
    }

    return in;
}


/*-----------------------------------------------*/
/*      atof: value determined by a string       */
/*-----------------------------------------------*/
bool NOMAD::Double::atof(const std::string &ss)
{
    std::string s = ss;
    NOMAD::toupper(s);

    if ( ss == NOMAD::Double::_undefStr
        || ss == NOMAD::DEFAULT_UNDEF_STR_1
        || ss == NOMAD::DEFAULT_UNDEF_STR_HYPHEN
        || ss == "-" + NOMAD::Double::_undefStr
        || ss == "-" + NOMAD::DEFAULT_UNDEF_STR_1 )
    {
        _value   = 0.0;
        _defined = false;
        return true;
    }

    if ( s == "INF" ||  s == "+INF" ||
        s == "NOMAD::INF" || s == "+NOMAD::INF" ||
        ss == NOMAD::Double::_infStr ||
        ss == ("+" + NOMAD::Double::_infStr) )
    {
        _value   = NOMAD::INF;
        _defined = true;
        return true;
    }

    if ( s == "-INF" || s == "-NOMAD::INF" || ss == ("-" + NOMAD::Double::_infStr) )
    {
        _value   = -NOMAD::INF;
        _defined = true;
        return true;
    }

    if ( s.empty() || (s.size() == 1 && !isdigit(s[0])) )
    {
        return false;
    }

    if ( !isdigit(s[0]) && s[0] != '+' && s[0] != '-' && s[0] != '.' )
    {
        return false;
    }

    size_t n = s.size();
    for ( size_t k = 1 ; k < n ; ++k )
    {
        if ( !isdigit(s[k]) && s[k] != '.' )
        {
            if ( s[k] == 'E' )
            {
                if ( s.size() == k+1 )
                {
                    return false;
                }
                ++k;
                if ( !isdigit(s[k]) && s[k] != '+' && s[k] != '-' )
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }

    *this = std::atof ( s.c_str() );
    return true;
}

/*-------------------------------------------------*/
/*  atof from a string that can begin with 'r' to  */
/*  indicate a proportion (relative value)         */
/*-------------------------------------------------*/
bool NOMAD::Double::relativeAtof ( const std::string & s , bool & relative )
{
    if ( std::toupper(s[0]) == 'R' )
    {
        relative  = true;
        std::string ss = s;
        ss.erase(ss.begin());
        if ( !atof(ss) )
            return false;
        return ( *this >= 0.0 );
    }
    relative = false;
    return atof(s);
}

/*-----------------------------------------------*/
/*            is the value an integer?           */
/*-----------------------------------------------*/
bool NOMAD::Double::isInteger ( void ) const
{
    if ( !_defined )
        return false;
    return ( NOMAD::Double(std::floor(_value))) == ( NOMAD::Double(std::ceil(_value)) );
}

/*-----------------------------------------------*/
/*             is the value binary ?             */
/*-----------------------------------------------*/
bool NOMAD::Double::isBinary ( void ) const
{
    if ( !_defined )
        return false;
    return ( NOMAD::Double(_value) == 0.0 || NOMAD::Double(_value) == 1.0 );
}

/*-------------------------------------*/
/*               d = d1/d2             */
/*-------------------------------------*/
const NOMAD::Double NOMAD::operator / ( const NOMAD::Double & d1 ,
                                       const NOMAD::Double & d2   )
{
    if ( !d1.isDefined() || !d2.isDefined() )
        throw NOMAD::Double::NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double: d1 / d2: d1 or d2 not defined" );
    if ( d2.todouble() == 0.0 )
        throw NOMAD::Double::InvalidValue ( "Double.cpp" , __LINE__ ,
                                            "NOMAD::Double: d1 / d2: division by zero" );
    return NOMAD::Double ( d1.todouble() / d2.todouble() );
}

/*-------------------------------------*/
/*                d1 += d2             */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator += ( const NOMAD::Double & d2 )
{
    if ( !_defined || !d2._defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double: d1 += d2: d1 or d2 not defined" );
    _value += d2._value;
    return *this;
}

/*-------------------------------------*/
/*               d1 -= d2              */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator -= ( const NOMAD::Double & d2 )
{
    if ( !_defined || !d2._defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double: d1 -= d2: d1 or d2 not defined" );
    _value -= d2._value;
    return *this;
}

/*-------------------------------------*/
/*                d1 *= d2             */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator *= ( const NOMAD::Double & d2 )
{
    if ( !_defined || !d2._defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double: d1 *= d2: d1 or d2 not defined" );
    _value *= d2._value;
    return *this;
}

/*-------------------------------------*/
/*               d1 /= d2              */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator /= ( const NOMAD::Double & d2 )
{
    if ( !_defined || !d2._defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double: d1 /= d2: d1 or d2 not defined" );
    if ( d2._value == 0.0 )
        throw InvalidValue ( "Double.cpp" , __LINE__ ,
                             "NOMAD::Double: d1 /= d2: division by zero" );
    _value /= d2._value;
    return *this;
}

/*-------------------------------------*/
/*                  ++d                */
/*-------------------------------------*/
NOMAD::Double & NOMAD::Double::operator++ ( void )
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ , "NOMAD::Double: ++d: d not defined" );
    _value += 1;
    return *this;
}

/*-------------------------------------*/
/*                  d++                */
/*-------------------------------------*/
NOMAD::Double NOMAD::Double::operator++ ( int n )
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ , "NOMAD::Double: d++: d not defined" );
    NOMAD::Double tmp = *this;
    if( n <= 0 )
        n = 1;
    _value += n;
    return tmp;
}

/*-------------------------------------*/
/*                --d                  */
/*-------------------------------------*/
NOMAD::Double & NOMAD::Double::operator-- ( void )
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ , "NOMAD::Double: --d: d not defined" );
    _value -= 1;
    return *this;
}

/*-------------------------------------*/
/*                  d--                */
/*-------------------------------------*/
NOMAD::Double NOMAD::Double::operator-- ( int n )
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double: d--: d not defined" );
    NOMAD::Double tmp = *this;
    if ( n <= 0 )
        n = 1;
    _value -= n;
    return tmp;
}

/*-------------------------------------*/
/*              operators =            */
/*-------------------------------------*/
NOMAD::Double & NOMAD::Double::operator= ( const NOMAD::Double & d )
{
    _value   = d._value;
    _defined = d._defined;
    return *this;
}

NOMAD::Double & NOMAD::Double::operator= ( double r )
{
    _value   = r;
    _defined = true;
    return *this;
}


#include "../nomad_nsbegin.hpp"
std::ostream& operator<<(std::ostream& os, const NOMAD::Double& d)
{
    if ( d.isDefined() )
    {
        double value = d.todouble();
        if ( value == NOMAD::INF )
        {
            os << NOMAD::Double::getInfStr();
        }
        else if ( value == -NOMAD::INF )
        {
            os << "-" << NOMAD::Double::getInfStr();
        }
        else if ( std::floor(value) == std::ceil(value) && fabs(value) < INT_MAX-1 )
        {
            os << static_cast<int>(value);
        }
        else
        {
            os << d.display(-1);
        }
    }
    else
    {
        os << NOMAD::Double::getUndefStr();
    }

    return os;
}
#include "../nomad_nsend.hpp"


std::string NOMAD::Double::display(const int prec, const size_t refWidth) const
{
    std::ostringstream oss;

    if(NOMAD::INF == _value)
    {
        return NOMAD::DEFAULT_INF_STR;
    }
    else if(NOMAD::INF == -_value)
    {
        std::string str = "-" + NOMAD::DEFAULT_INF_STR;
        return str;
    }

    // Set the number of digits after the point (ignore if prec < 0).
    if (prec >= 0)
    {
        // Fixed
        oss.setf(std::ios::fixed, std::ios::floatfield);
        size_t width = 0;
        std::ostringstream osstemp;
        if (_defined)
        {
            osstemp.precision(NOMAD::DISPLAY_PRECISION_FULL);
            osstemp << _value;
        }
        else
        {
            osstemp << NOMAD::DEFAULT_UNDEF_STR_HYPHEN;
        }
        std::string s = osstemp.str();

        size_t spacePadding = 0;
        NOMAD::getFormat(s, prec, width, spacePadding);

        if (refWidth != 0)
        {
            // Override computed width
            width = refWidth;
        }

        // If the number of decimals in _value is greater then prec, then
        // output it as is so it gets truncated.
        // Ex: 447.000774493 -> 447.000774
        // If it is smaller, use the string and complete with space padding.
        // Ex. -1878.99 -> "-1878.99    "
        size_t nbDec = NOMAD::nbDecimals(s);
        if (nbDec >= (size_t)prec)
        {
            oss << std::setprecision(prec) << std::setw(static_cast<int>(width)) << _value;
        }
        else
        {
            // Add extra spaces to s before adding s to oss.
            for (size_t i = 0; i < spacePadding && i < width; i++)
            {
                s += " ";
            }
            oss << std::setw(static_cast<int>(width)) << s;
        }

        // Replace superfluous 0's with spaces
        size_t pos0 = oss.str().find_last_not_of("0");
        if (std::string::npos != pos0 && nbDec > 0)
        {
            s = oss.str();
            if ('.' == s[pos0]) { pos0++; } // Leave an extra '0' after the decimal point
            pos0++; // Start replacing from first non-0 char
            size_t nbRep = s.size() - pos0;
            std::string sSpaces(nbRep, ' ');
            s.replace(pos0, nbRep, sSpaces);
            oss.str(s);
        }
    }
    else if (_defined)
    {
        // Just output value.
        oss << _value;
    }
    else
    {
        oss << NOMAD::DEFAULT_UNDEF_STR_HYPHEN;
    }

    return oss.str();
}


/*------------------------------------------*/
/*              display with format         */
/*------------------------------------------*/
// This method is taken from NOMAD 3.
std::string NOMAD::Double::display(const std::string& format) const
{
    std::ostringstream oss;

    // interpret the format:
    // ---------------------

    // %f      w=-1 prec=-1 c='f'
    // %4.5f   w= 4 prec= 5 c='f'
    // %4f     w= 4 prec= 1 c='f'
    // %.5f    w=-1 prec= 5 c='f'
    // %.f     w=-1 prec= 0 c='f'

    // c may be in 'e', 'E', 'f', 'g', 'G', 'd', or 'i'

    // e Scientific notation (mantise/exponent) using e character 3.9265e+2
    // E Scientific notation (mantise/exponent) using E character 3.9265E+2
    // f Decimal floating point                                   392.65
    // g Use the shorter of %e or %f                              392.65
    // G Use the shorter of %E or %f                              392.65
    // d or i Integer rounded value                               393

    std::string format2 = format;

    int  w    = -1;
    int  prec = -1;
    char c    =  0;

    if (!format2.empty() && format2[0] == '%')
    {
        size_t n = format2.size();

        c = format2[n-1];

        if ( c != 'e' && c != 'E' && c != 'f' && c != 'g' && c != 'G' && c != 'd' && c != 'i' )
        {
            c = ( std::floor(_value) == std::ceil(_value) && fabs(_value) < INT_MAX-1 )
                    ? 'd' : 'f';
            format2.push_back(c);
            ++n;
        }

        if (n > 2)
        {
            std::string sw, sprec;

            size_t k = format2.find(".");

            if ( k > 0 && k < n-1 )
            {
                if ( n==3 )
                {
                    sprec = "0";
                }
                else
                {
                    if ( k > 1 )
                    {
                        sw = format2.substr ( 1 , k-1 );
                    }
                    sprec = format2.substr ( k+1 , n-k-2 );
                }
            }
            else
            {
                sw = format2.substr ( 1 , n-2 );
            }
            if ( !NOMAD::atoi ( sw , w ) )
            {
                w = -1;
            }

            if ( !NOMAD::atoi ( sprec , prec ) )
            {
                prec = -1;
            }
        }

        if ( c=='d' || c=='i' )
        {
            prec = 0;
        }
    }

    // display the value:
    oss << std::setw(w);
    if (_defined)
    {
        if ( _value == NOMAD::INF )
        {
            oss << NOMAD::DEFAULT_INF_STR;
        }
        else if ( c=='d' || c=='i' ||
                 ( format2.empty() &&
                  std::floor(_value) == std::ceil(_value) && fabs(_value) < INT_MAX-1 ) )
        {
            oss << roundd();
        }
        else
        {
            int                     old_prec  = static_cast<int>(oss.precision());
            std::ios_base::fmtflags old_flags = oss.flags();

            if ( prec >= 0 )
            {
                oss.precision(prec);
            }

            if ( c == 'f' )
            {
                oss.setf(std::ios::fixed);
                oss << _value;
            }
            else if ( c == 'e' )
            {
                oss.unsetf ( std::ios::fixed );
                oss.setf   ( std::ios::scientific );
                oss << _value;
            }
            else if ( c == 'E' )
            {
                oss.unsetf ( std::ios::fixed );
                oss.setf   ( std::ios::scientific | std::ios::uppercase );
                oss << _value;
            }

            else if ( c == 'g' )
            {
                // use ostringstream for some format because the oss << _value is not working in some situations.
                std::ostringstream streamS,streamF;
                streamS.precision ( prec );
                streamF.precision ( prec );
                streamF.unsetf(std::ios::scientific);
                streamF.setf( std::ios::fixed );
                streamS.unsetf(std::ios::fixed);
                streamS.setf( std::ios::scientific);
                streamS << _value;
                streamF << _value;

                if ( streamS.str().length() < streamF.str().length() )
                {
                    oss << streamS.str();
                }
                else
                {
                    oss << streamF.str();
                }
            }

            else if ( c == 'G' )
            {
                std::ostringstream streamS,streamF;
                streamS.precision ( prec );
                streamF.precision ( prec );
                streamF.unsetf(std::ios::scientific);
                streamF.setf( std::ios::fixed | std::ios::uppercase );
                streamS.unsetf(std::ios::fixed);
                streamS.setf( std::ios::scientific | std::ios::uppercase );
                streamS << _value ;
                streamF << _value ;
                if (streamS.str().length() < streamF.str().length())
                {
                    oss << streamS.str();
                }
                else
                {
                    oss << streamF.str();
                }
            }

            oss.precision ( old_prec  );
            oss.flags     ( old_flags );
        }
    }
    else
    {
        oss << NOMAD::DEFAULT_UNDEF_STR;
    }

    return oss.str();
}


/*------------------------------------------*/
/*                round to int              */
/*------------------------------------------*/
int NOMAD::Double::round ( void ) const
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::round(): value not defined" );

    double d = (_value < 0.0 ? -std::floor(.5-_value) : std::floor(.5+_value));

    if ( d > NOMAD::P_INF_INT || d < NOMAD::M_INF_INT )
        throw InvalidValue ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::round(): value cannot be rounded to integer because it is outside of range" );

    return static_cast<int> (d);
}

// Return the number of decimals of a positive Double.
// Ex. 123.4567 -> 4
//     123 -> 0
//     0.000 -> 3
std::size_t NOMAD::Double::nbDecimals() const
{
    std::size_t nbDec;

    if (_value < _epsilon)
    {
        std::string str = "Error: nbDecimals of number smaller than EPSILON is not supported";
        throw NOMAD::Exception(__FILE__, __LINE__, str);
    }

    NOMAD::Double rem( _value );
    int dec = (int)std::floor(log10(rem.todouble()));
    rem -= pow(10, dec);
    while (rem._value >= _epsilon)
    {
        dec = (int)std::floor(log10(rem.todouble()));
        rem -= pow(10, dec);
    }
    if (dec > 0)
    {
        nbDec = 0;
    }
    else
    {
        nbDec = -dec;
    }

    return nbDec;
}


/*------------------------------------------*/
/*              round to double             */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::roundd () const
{
    if ( !_defined )
    {
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::round(): value not defined" );
    }

    return (_value < 0.0 ? -std::floor(.5-_value) : std::floor(.5+_value));

}


/*------------------------------------------*/
/*                  Ceil                  */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::ceil ( void ) const
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::ceil(): value not defined" );
    return NOMAD::Double( std::ceil(_value) );
}

/*------------------------------------------*/
/*                  Floor                  */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::floor ( void ) const
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::floor(): value not defined" );
    return NOMAD::Double( std::floor(_value) );
}

/*------------------------------------------*/
/*                    abs                   */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::abs ( void ) const
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::abs(): value not defined" );
    return fabs ( _value );
}

/*------------------------------------------*/
/*                  square                  */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::pow2 ( void ) const
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::pow2(): value not defined" );
    return pow ( _value , 2 );
}

/*------------------------------------------*/
/*                square root               */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::sqrt ( void ) const
{
    if ( !_defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::sqrt(): value not defined" );
    if ( *this < 0.0 )
        throw NOMAD::Double::InvalidValue ( "Double.cpp" , __LINE__ ,
                                            "NOMAD::Double::sqrt(x): x < 0" );

    return std::sqrt ( _value );
}

/*---------------------------------------------*/
/*  relative error with another NOMAD::Double  */
/*---------------------------------------------*/
//
// This error computation is based on:
//
// A. Ziv. Relative distance–an error measure in round-off error analysis. Mathematics
// of Computation, 39(160):563–569, 1982. doi:10.1090/S0025-5718-1982-0669649-2.
//
// The error will be in [0;2]
const NOMAD::Double NOMAD::Double::relErr ( const NOMAD::Double & x ) const
{
    if ( !_defined || !x._defined )
        throw NotDefined ( "Double.cpp" , __LINE__ ,
                           "NOMAD::Double::rel_err(): one of the values is not defined" );

    // 1. test if x==y:
    if ( this == &x || _value == x._value )
        return 0.0;

    double diff = fabs ( _value - x._value );

    // 2. test if one of the values is zero:
    if ( _value == 0.0 || x._value == 0.0 )
    {

        // we return min{2,|x-y|} (instead of 1):
        if ( diff > 2.0 )
            return 2.0;
        return diff;
    }

    // 3. compute the original error:
    double a   = fabs ( _value   );
    double b   = fabs ( x._value );
    double err = diff / ( (a<b) ? b : a );

    // 4. test if we have opposite signs:
    if ( _value * x._value < 0.0 )
    {

        // the original error gives err in ]1;2] : we check if |x-y| < 1
        // and if so we return |x-y| :
        if ( diff < 1.0 )
            return diff;
    }

    // we return the original error:
    return err;
}


bool NOMAD::Double::isMultipleOf(const NOMAD::Double &granularity) const
{
    bool isMult = false;

    if (!isDefined())
    {
        // Double is not defined - Return false.
        isMult = false;
    }
    else if (isDefined() && abs().todouble() <= _epsilon)
    {
        // Value 0 is always a multiple of granularity.
        isMult = true;
    }

    else if (granularity.isDefined() && (0.0 < granularity.todouble()))
    {
        if (isDefined())
        {
            // Verify that:
            //     round (value / granularity) * granularity ~= value
            // i.e. that value / granularity is (roughly) an integer.
            NOMAD::Double mult = (_value / granularity).roundd();
            NOMAD::Double verif_value = mult * granularity;
            if ((_value - verif_value).abs().todouble() < mult.abs().todouble() * NOMAD::DEFAULT_EPSILON)
            {
                isMult = true;
            }
        }
        else if (toBeDefined())
        {
            // Assume the granularity check will be done once the value is set.
            isMult = true;
        }
    }
    else if (granularity.isDefined() && (granularity < 0.0))
    {
        // Ignore negative granularity.
        isMult = false;
    }
    else
    {
        // Ignore granularity that is undefined or too small.
        isMult = true;
    }

    return isMult;
}


const NOMAD::Double NOMAD::Double::nextMult(const NOMAD::Double &granularity) const
{
    NOMAD::Double d;

    if (!granularity.isDefined() || !isDefined() || (granularity <= 0.0) || isMultipleOf(granularity))
    {
        d = _value;
    }
    else
    {
        // granularity > 0, and _value is not a multiple of granularity.
        // Adjust value with granularity
        int granMult = (int)(_value / granularity.todouble());
        if (_value > 0)
        {
            granMult++;
        }
        double bigGranExp = pow(10, granularity.nbDecimals());
        int bigGran = (int)(granularity.todouble() * bigGranExp);
        d = granMult * bigGran / bigGranExp;
    }

    return d;
}


const NOMAD::Double NOMAD::Double::previousMult(const NOMAD::Double &granularity) const
{
    NOMAD::Double d;

    if (!granularity.isDefined() || !isDefined() || (granularity <= 0.0) || isMultipleOf(granularity))
    {
        d = _value;
    }
    else
    {
        // granularity > 0, and _value is not a multiple of granularity.
        // Adjust value with granularity
        int granMult = (int)(_value / granularity.todouble());
        if (_value < 0)
        {
            granMult--;
        }
        double bigGranExp = pow(10, granularity.nbDecimals());
        int bigGran = (int)(granularity.todouble() * bigGranExp);
        d = granMult * bigGran / bigGranExp;
    }

    return d;
}
