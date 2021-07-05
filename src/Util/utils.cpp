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
 \file   utils.cpp
 \brief  Utility functions
 \author Sebastien Le Digabel, modified by Viviane Rochon Montplaisir
 \date   March 2017
 \see    utils.hpp
 */
#include <algorithm>    // for for_each
#include <cctype>       // for toupper
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"

/*-----------------------------------------------------------------*/
/*                         NOMAD::itos                             */
/*-----------------------------------------------------------------*/
std::string NOMAD::itos ( const int i )
{
    std::ostringstream oss;
    oss << i;
    return oss.str();
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::itos                             */
/*-----------------------------------------------------------------*/
std::string NOMAD::itos ( const size_t i )
{
    std::ostringstream oss;
    oss << i;
    return oss.str();
}


/*-----------------------------------------------------------------*/
/*                         NOMAD::toupper - 1/2                    */
/*-----------------------------------------------------------------*/
void NOMAD::toupper(std::string &s)
{
    // Warning: strings like "é" do not get converted.
    // Even when using locale.
    // This is a widechar issue, not investigating in it right now.
    // There is no easy way to show a clear warning, so do not do
    // anything about it - just file it is as a known issue.

// So 98
//    size_t ns = s.size();
//    for ( size_t i = 0 ; i < ns ; ++i )
//    {
//        s[i] = std::toupper(s[i]);
//    }

// modern obfuscated C++
    for_each(s.begin(), s.end(), [](char& in){ in = std::toupper(in); });
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::toupper - 2/2                    */
/*-----------------------------------------------------------------*/
void NOMAD::toupper(std::list<std::string> &ls)
{
    std::list<std::string>::iterator       it;
    std::list<std::string>::const_iterator end = ls.end();
    for ( it = ls.begin() ; it != end ; ++it )
    {
        NOMAD::toupper(*it);
    }
}


/*-----------------------------------------------------------------*/
/*                             NOMAD::trim                         */
/*-----------------------------------------------------------------*/
void NOMAD::trim(std::string &s)
{
    // Trim extra spaces at the beginning
    size_t space_index = s.find(' ');
    while (s.size() > 0 && 0 == space_index)
    {
        s.replace(0, 1, "");
        space_index = s.find(' ');
    }

    // Trim extra spaces at the end
    size_t space_rindex = s.rfind(' ');
    while (s.size() > 0 && s.size()-1 == space_rindex)
    {
        s.replace(space_rindex, 1, "");
        space_rindex = s.rfind(' ');
    }
}


/*-----------------------------------------------------------------*/
/*                             NOMAD::atoi                         */
/*-----------------------------------------------------------------*/
bool NOMAD::atoi ( const std::string & s , int & i )
{
    i = -1;
    if ( s.empty() )
    {
        return false;
    }

    size_t n = s.size();

    if ( s[0] == '-' )
    {
        if ( n > 1 && s[1] == '-' )
            return false;
        std::string ss = s;
        ss.erase(ss.begin());
        if ( NOMAD::atoi ( ss , i ) )
        {
            i = -i;
            return true;
        }
        return false;
    }

    std::string ts = s;
    NOMAD::toupper(ts);
    if ( ts == "INF" || ts == "+INF" )
    {
        i = NOMAD::P_INF_INT;
        return true;
    }
    if ( ts == "-INF" )
    {
        i = NOMAD::M_INF_INT;
        return true;
    }

    for ( size_t k = 0 ; k < n ; ++k )
    {
        if ( !isdigit(s[k]) )
        {
            return false;
        }
    }
    i = std::atoi(s.c_str());
    return true;
}

/*-----------------------------------------------------------------*/
/*                             NOMAD::atoi                         */
/*-----------------------------------------------------------------*/
bool NOMAD::atost ( const std::string & s , size_t & st )
{
    st = NOMAD::INF_SIZE_T;

    if ( s.empty() )
    {
        return false;
    }

    std::string ts = s;
    NOMAD::toupper(ts);
    if ( ts == "INF" || ts == "+INF" )
    {
        st = NOMAD::INF_SIZE_T;
        return true;
    }

    int i;
    bool success = NOMAD::atoi(s,i);
    if ( success )
    {
        if ( i < 0 )
            throw NOMAD::Exception(__FILE__, __LINE__, "Invalid value for size_t. Value must be >0");
        st = static_cast<size_t>(i);
    }

    return success;
}


bool NOMAD::atoi(const char c, int & i)
{
    std::string s = "-";
    s[0] = c;
    return NOMAD::atoi(s,i);
}





/*-----------------------------------------------------------*/
/* Convert a success type to a string (used for debugging).  */
/*-----------------------------------------------------------*/
std::string NOMAD::enumStr(NOMAD::SuccessType success)
{
    std::string str;
    switch (success)
    {
        case NOMAD::SuccessType::NOT_EVALUATED:
            str = "Not evaluated yet";
            break;
        case NOMAD::SuccessType::UNSUCCESSFUL:
            str = "Failure";
            break;
        case NOMAD::SuccessType::PARTIAL_SUCCESS:
            str = "Partial success (improving)";    // h is better, f is worse.
            break;
        case NOMAD::SuccessType::FULL_SUCCESS:
            str = "Full success (dominating)";
            break;
        default:
            str = "Error: Enum for success type is not recognized";
            throw NOMAD::Exception(__FILE__, __LINE__, str);
    }

    return str;
}

// Convert a string to index range
bool NOMAD::stringToIndexRange(const std::string & s           ,
                               int               & i           ,
                               int               & j           ,
                               bool               check_order   )
{
    if ( s.empty() )
        return false;

// For now we accept only range i-j and -j
//    if ( s == "*" )
//    {
//        if ( !n )
//            return false;
//        i = 0;
//        j = *n-1;
//        return true;
//    }

    if ( s[0] == '-' )
    {

        size_t ns = s.size();
        if ( ns > 1 && s[1] == '-' )
            return false;

        std::string ss = s;
        ss.erase ( ss.begin() );

        if ( NOMAD::stringToIndexRange ( ss , i , j , false ) )
        {
            i = -i;
            return true;
        }
        return false;
    }

    std::istringstream in (s);
    std::string        s1;

    getline ( in , s1 , '-' );

    if (in.fail())
        return false;

    size_t k , n1 = s1.size();

    if ( n1 >= s.size() - 1 )
    {
        for ( k = 0 ; k < n1 ; ++k )
            if (!isdigit(s1[k]))
                return false;
        size_t ist = (size_t)i;
        if ( ! atost ( s1 , ist ) )
            return false;
        i = (int)ist;
        if ( n1 == s.size() )
        {
            j = i;
            return true;
        }
        return false;
    }

    std::string s2;
    getline (in, s2);

    if (in.fail())
        return false;

    size_t n2 = s2.size();
    for ( k = 0 ; k < n2 ; ++k )
        if ( !isdigit(s2[k]) )
            return false;

    size_t ist = (size_t)i;
    size_t jst = (size_t)j;
    if ( ! atost ( s1, ist ) || ! atost ( s2 , jst ) )
        return false;
    i = (int)ist;
    j = (int)jst;

    return !check_order || i <= j;
}


// Convert a string to bool
bool NOMAD::stringToBool(const std::string &string)
{
    bool ret = false;
    std::string s = string;
    NOMAD::toupper(s);

    if ( s == "Y" || s == "YES" || s == "1" || s == "TRUE" )
    {
        ret = true;
    }
    else if ( s == "N" || s == "NO" || s == "0" || s == "FALSE" )
    {
        ret = false;
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized string for bool: " + s);
    }

    return ret;
}


// Convert a bool to string
std::string NOMAD::boolToString(bool boolean)
{
    return (boolean) ? "true" : "false";
}


// Return the number of decimals of a string representing a double.
// It is assumed that the number is properly formed.
// Ex. 123.4567 -> 4
//     123 -> 0
//     0.000 -> 3
std::size_t NOMAD::nbDecimals(const std::string& s)
{
    std::size_t nbDec;

    std::size_t ptPos = s.rfind(".");
    if (std::string::npos == ptPos)
    {
        nbDec = 0;
    }
    else
    {
        nbDec = s.size() - ptPos - 1;
    }
    return nbDec;
}


// s is a string representing a double.
// prec = number of digits to show after the point.
// spacePadding = recommended space padding after s.
// width = total width
// Ex. 5.123456789 with a prec of 6 shows as 5.123456.
// Ex. 5.12 with a prec of 6 will return as "5.12    " (with padding).
void NOMAD::getFormat(const std::string &s, const size_t prec, size_t &width, size_t &spacePadding)
{
    // These values non modifiable for now.
    const size_t nbDigitsBeforePoint = NOMAD::NB_DIGITS_BEFORE_POINT;
    const size_t intWidth = NOMAD::INT_DISPLAY_WIDTH;

    if (0 == prec)
    {
        // Integer.
        width = intWidth;
    }
    else
    {
        // Double
        width = nbDigitsBeforePoint + 1 + prec;
        size_t pointPos = s.find(".");
        spacePadding = 1 + prec; // Including one space for the decimal point
        if (pointPos != std::string::npos)
        {
            // Compute padding as needed.
            // Ex. if precision is 6 but the string number is "27.03" -> pad to "27.03    ".
            spacePadding -= (s.size() - pointPos);
            if (spacePadding >= width)
            {
                spacePadding = 0;
            }
        }
    }
}


bool NOMAD::separateFormat(const std::string &s, std::string &format, std::string &tag)
{
    // NOTE: format found here is used in StatsInfo, for display of OBJ, CONS_H
    // and H_MAX (or any Double). It is not used for values of type ArrayOfDouble
    // or its derived classes.

    format = "";
    tag = s;
    const std::string formatLetters = "eEfgGdi";
    const std::string allLetters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    bool doSeparate = false;

    if ('%' == s[0])
    {
        std::size_t firstLetter = s.find_first_of(allLetters, 1);
        std::size_t firstFormatLetter = s.find_first_of(formatLetters, 1);
        std::size_t cutIndex = firstLetter;
        if (std::string::npos != firstFormatLetter
            && firstLetter == firstFormatLetter)
        {
            cutIndex += 1;
        }
        if (std::string::npos != cutIndex)
        {
            std::string tempFormat = s.substr(0, cutIndex);
            doSeparate = validFormat(tempFormat);
            if (doSeparate)
            {
                format = tempFormat;
                tag = s.substr(cutIndex, s.length() - cutIndex);
            }
        }
    }

    return doSeparate;
}


// Pasted from NOMAD 3:
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

// %4.2 is also valid (no format letter). Append f.
bool NOMAD::validFormat(std::string &s)
{
    const std::string formatLetters = "eEfgGdi";
    bool isValid = true;

    // Verify size
    if (s.length() < 2)
    {
        isValid = false;
    }

    else
    {
        // Append f if last char is a digit.
        if (std::isdigit(s[s.length()-1]))
        {
            s = s + "f";
        }
        size_t indexFormatLetter = s.find_first_of(formatLetters, 1);

        // Verify we do have a format letter
        if (std::string::npos == indexFormatLetter)
        {
            isValid = false;
        }
        // Verify first character is '%'
        else if (s[0] != '%')
        {
            isValid = false;
        }
        // Verify last character is a formatting character
        else if (indexFormatLetter < s.length()-1)
        {
            isValid = false;
        }
        // Verify in between we have a decimal number
        else
        {
            bool pointEncountered = false;
            for (size_t i = 1; i < indexFormatLetter; i++)
            {
                if (std::isdigit(s[i]))
                {
                    // ok
                }
                else if (s[i] == '.')
                {
                    if (!pointEncountered)
                    {
                        pointEncountered = true;
                    }
                    else
                    {
                        // more than one point
                        isValid = false;
                    }
                }
                else
                {
                    isValid = false;
                }
            }
        }
    }

    return isValid;
}


int NOMAD::getThreadNum()
{
    int threadNum = 0;
#ifdef _OPENMP
    threadNum = omp_get_thread_num();
#endif // _OPENMP
    return threadNum;
}
