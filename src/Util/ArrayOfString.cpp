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
 \file   ArrayOfString.cpp
 \brief  Custom class for array of strings (implementation)
 \author Viviane Rochon Montplaisir
 \date   February 2018
 \see    ArrayOfString.hpp
 */
#include <algorithm>
#include <sstream>

#include "../Util/ArrayOfString.hpp"
#include "../Util/Exception.hpp"

std::ostream& NOMAD::operator<<(std::ostream& out, const NOMAD::ArrayOfString& arrayOfString)
{
    out << arrayOfString.display();

    return out;
}


/*-----------------------------------------------------------*/
/*                        Constructors                       */
/*-----------------------------------------------------------*/
NOMAD::ArrayOfString::ArrayOfString(size_t n, const std::string& initString)
  : _array()
{
    for (size_t i = 0; i < n; i++)
    {
        _array.push_back(initString);
    }
}


NOMAD::ArrayOfString::ArrayOfString(const std::string & inputString,
                                    const std::string & separators)
    : _array()
{
    _array = splitString(inputString, separators);
}


/*-----------------------------------------------*/
/*                    destructor                 */
/*-----------------------------------------------*/
NOMAD::ArrayOfString::~ArrayOfString()
{
}


/*-----------------------------------------------------------*/
/*                       '[]' operators                      */
/*-----------------------------------------------------------*/
const std::string& NOMAD::ArrayOfString::operator[](size_t index) const
{
    if ( index >= _array.size())
    {
        std::ostringstream oss;
        oss << "ArrayOfString: index = " << index << " is out of bounds [0, " << _array.size()-1 << "]";
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }

    return _array[index];
}


/*------------------------------------------------*/
/* Replace string s at position index.            */
/* Throw an error if the index is too large.      */
/* Note: operator[] returns const only, so to     */
/* modify an array, replace has to be used.       */
/*------------------------------------------------*/
void NOMAD::ArrayOfString::replace(const size_t index, const std::string& s)
{
    if ( index >= _array.size())
    {
        std::ostringstream oss;
        oss << "ArrayOfString: index = " << index << " is out of bounds [0, " << _array.size()-1 << "]";
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }
    _array[index] = s;
}


/*------------------------------------------------*/
/*  Erase element at position index               */
/*  Return true if string was removed from array. */
/*------------------------------------------------*/
bool NOMAD::ArrayOfString::erase(const size_t index)
{
    size_t k = 0;
    std::vector<std::string>::const_iterator it;

    for (it = _array.begin(); it != _array.end(); ++it)
    {
        if (index == k)
        {
            _array.erase(it);
            return true;
        }
        k++;
    }

    return false;
}


// Return first index of string s in array.
// Return -1 if s not found.
int NOMAD::ArrayOfString::find(const std::string& s) const
{
    int index = -1;

    for (size_t i = 0; i < _array.size(); i++)
    {
        if (s == _array[i])
        {
            index = (int)i;
            break;
        }
    }

    return index;
}


/*------------*/
/* Comparison */
/*------------*/
bool NOMAD::ArrayOfString::operator== (const NOMAD::ArrayOfString &array) const
{
    return (_array == array._array);
}


/*---------------*/
/*    Display    */
/*---------------*/
std::string NOMAD::ArrayOfString::display() const
{
    std::string s;

    if ( size() == 0 )
        s += " - ";

    for (size_t i = 0; i < size(); i++)
    {
        if (i > 0)
        {
            s += " ";
        }
        s += _array[i];
    }

    return s;
}


/*-----------------------------------*/
/* Split                             */
/* - Used as helper for constructor. */
/*-----------------------------------*/

std::vector<std::string> const
NOMAD::ArrayOfString::splitString(const std::string & inputString,
                                  const std::string & separators)
{
    std::vector<std::string> array;

    if (inputString.size() == 0)
    {
        return array;
    }

    size_t splitIndex = 0, index1 = 0, length;

    // Fill array with strings found in inputString.
    // inputString is a list of std::strings, separated by any char in separators.
    // Sample inputString:
    // OBJ PB EB

    while (splitIndex != std::string::npos)
    {
        // Find index of first non-separator.
        index1 = inputString.find_first_not_of(separators, index1);
        if (index1 == std::string::npos)
        {
            break;
        }
        splitIndex = inputString.find_first_of(separators, index1+1);
        if (splitIndex == std::string::npos)
        {
            length = inputString.size() - index1;
        }
        else
        {
            length = splitIndex - index1;
        }

        std::string si = inputString.substr(index1, length);
        array.push_back(si);

        index1 = splitIndex + 1;    // Start looking after separator
    }

    return array;
}


NOMAD::ArrayOfString NOMAD::ArrayOfString::combineAndAddPadding(const NOMAD::ArrayOfString & s1, const NOMAD::ArrayOfString & s2)
{
    // Add a padding on the first string of each pair of string. The padding must assure that the second strings are all aligned.
    // Also add a return at the end of the second string.

    size_t sizeS1 = s1.size();
    if (sizeS1 != s2.size() )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "s1 and s2 must have the same of number of elements.");
    }

    // Max length of first+second elements in s
    size_t maxL = 0;
    for (size_t e = 0 ; e < sizeS1 ; e++)
    {
        maxL = std::max(maxL,s1[e].length()+s2[e].length());
    }

    // Pad the first element to get the same overall length for all + combine first and second + add return
    NOMAD::ArrayOfString paddedString("\n");
    for (size_t e = 0 ; e < sizeS1 ; e++)
    {
        size_t padL = 1 + maxL - s1[e].length() - s2[e].length(); // Add one ' ' at least (for space after colon)
        std::string padS = s1[e];
        padS.insert(s1[e].length(), padL,' ');
        padS += s2[e]+'\n';
        paddedString.add(padS);
    }

    return paddedString;
}





