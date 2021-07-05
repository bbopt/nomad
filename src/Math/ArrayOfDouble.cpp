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
 \file   ArrayOfDouble.cpp
 \brief  Custom class for array of Doubles (implementation)
 \author Viviane Rochon Montplaisir
 \date   October 2017
 \see    ArrayOfDouble.hpp
 */
#include "../Math/ArrayOfDouble.hpp"
#include <iomanip>  // For std::setprecision, std::setw

// Initialize static variables
const std::string NOMAD::ArrayOfDouble::pStart = "(";
const std::string NOMAD::ArrayOfDouble::pEnd = ")";

std::ostream& NOMAD::operator<<(std::ostream& out, const NOMAD::ArrayOfDouble& arrayOfDouble)
{
    out << arrayOfDouble.display();

    return out;
}


std::istream& NOMAD::operator>>(std::istream& in, NOMAD::ArrayOfDouble& coords)
{
    size_t n = coords.size();
    for (size_t k = 0; k < n; ++k)
    {
        in >> coords[k];
    }
    if (in.fail() && !in.eof())
    {
        std::string err = "ArrayOfDouble: bad input for operator>>";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    return in;
}


/*-----------------------------------------------------------*/
/*                         constructor                       */
/*-----------------------------------------------------------*/
NOMAD::ArrayOfDouble::ArrayOfDouble(size_t n, const NOMAD::Double& d)
  : _n(n),
    _array(nullptr)
{
    if (_n > 0)
    {
        _array = new NOMAD::Double [_n];
        if (d.isDefined())
        {
            std::fill (_array, _array + _n, d);
        }
    }
    else
    {
        _n = 0;
    }
}

NOMAD::ArrayOfDouble::ArrayOfDouble(const std::vector<double> & v)
  : _n(v.size()),
    _array(nullptr)
{
    if (_n > 0)
    {
        _array = new NOMAD::Double[_n];
        for (size_t k = 0; k < _n; k++)
        {
            _array[k] = v[k];
        }
    }
    else
    {
        _n = 0;
    }
}


/*-----------------------------------------------------------*/
/*                        copy constructor                   */
/*-----------------------------------------------------------*/
NOMAD::ArrayOfDouble::ArrayOfDouble(const NOMAD::ArrayOfDouble &coord)
  : _n(coord._n),
    _array(nullptr)
{
    if (_n > 0)
    {
        NOMAD::Double       * array1 =  _array = new NOMAD::Double [_n];
        const NOMAD::Double * array2 = coord._array;
        for (size_t k = 0; k < _n; ++k, ++array1, ++array2)
        {
            *array1 = *array2;
        }
    }
}


/*-----------------------------------------------*/
/*                    destructor                 */
/*-----------------------------------------------*/
NOMAD::ArrayOfDouble::~ArrayOfDouble ()
{
    delete [] _array;
}


/*-----------------------------------------------*/
/*   This method changes the array's dimension   */
/*   and sets all values to d                    */
/*-----------------------------------------------*/
void NOMAD::ArrayOfDouble::reset (size_t n, const NOMAD::Double &d)
{
    if (n == 0)
    {
        _n = 0;
        delete [] _array;
        _array = nullptr;
    }
    else
    {
        _n = n;
        delete [] _array;
        _array = new NOMAD::Double [_n];

        if (d.isDefined())
        {
            std::fill(_array, _array + _n, d);
        }
    }
}


/*-----------------------------------------------*/
/*  This method changes the array's dimension.   */
/*  The values are kept.                         */
/*-----------------------------------------------*/
void NOMAD::ArrayOfDouble::resize(size_t n, const NOMAD::Double &d)
{
    if (n == _n)
    {
        return;
    }

    if (n == 0)
    {
        _n = 0;
        delete [] _array;
        _array = nullptr;
        return;
    }

    NOMAD::Double *newArray = new NOMAD::Double[n];
    if (_array)
    {
        size_t min = ( n < _n ) ? n : _n;

        NOMAD::Double       * array1 = newArray;
        const NOMAD::Double * array2 = _array;

        for (size_t i = 0; i < min; ++i, ++array1, ++array2)
        {
            *array1 = *array2;
        }
        if (d.isDefined())
        {
            std::fill(newArray + min, newArray + n, d);
        }

        delete [] _array;
    }
    _array  = newArray;
    _n      = n;
}


/*-----------------------------------------------------------*/
/*                       '[]' operators                      */
/*-----------------------------------------------------------*/

NOMAD::Double& NOMAD::ArrayOfDouble::operator[](size_t i) const
{
    if (!_array)
    {
        std::string err = "ArrayOfDouble: Array is not defined";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    if (i >= _n)
    {
        std::ostringstream oss;
        oss << "ArrayOfDouble: i = " << i << " is out of bounds [0, " << _n-1 << "]";
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }

    return _array[i];
}


/*-----------------------------------------------------------*/
/*                        snap to bounds                     */
/*-----------------------------------------------------------*/
void NOMAD::ArrayOfDouble::snapToBounds(const NOMAD::ArrayOfDouble &lowerBound,
                                        const NOMAD::ArrayOfDouble &upperBound)
{
    size_t n = size();

    if (!isComplete())
    {
        std::string err("snapToBounds: ");
        err += "ArrayOfDouble is not completely defined.";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    if (lowerBound.size() != n || upperBound.size() != n)
    {
        std::string err = "snapToBounds: ";
        err += "Inconsistent dimension for bounds. Expecting ";
        err += std::to_string(n);
        err += " but sizes are " + std::to_string(lowerBound.size());
        err += " and " + std::to_string(upperBound.size()) + ".";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (size_t i = 0; i < n; ++i)
    {
        if (lowerBound.isDefined() && lowerBound[i].isDefined() && _array[i] < lowerBound[i])
        {
            _array[i] = lowerBound[i]; // True snap, ignore mesh
        }
        if (upperBound.isDefined() && upperBound[i].isDefined() && upperBound[i] < _array[i])
        {
            _array[i] = upperBound[i]; // True snap, ignore mesh
        }
    }
}


bool NOMAD::ArrayOfDouble::inBounds(const NOMAD::ArrayOfDouble &lowerBound,
                                    const NOMAD::ArrayOfDouble &upperBound) const
{
    bool ret = true;
    for (size_t i = 0; i < _n && ret; ++i)
    {
        if (!_array[i].isDefined())
        {
            ret = false;
        }
        else if (lowerBound[i].isDefined() && _array[i] < lowerBound[i])
        {
            ret = false;
        }
        else if (upperBound[i].isDefined() && _array[i] > upperBound[i])
        {
            ret = false;
        }
    }

    return ret;
}


void NOMAD::ArrayOfDouble::readValuesAsArray(const NOMAD::ArrayOfString& strDouble)
{
    // Fill array with values found in strDouble.
    // strDouble is a list of n (DIMENSION) NOMAD::Doubles.
    // Can be starting with a pStart and ending with a pEnd.
    // Other possible format: * double.
    // Sample valueString:
    // ( -6.0 -6.0 -5.0 -6.0 -6.0 )
    // ( 5.0 6.0 7.0 NaN NaN )
    // ( 1 1 1 1 1 - - - - - - - - - - )
    // 1-4 0.5
    // * -6.0
    //

    std::string valueString = strDouble.display();  // For informations display

    NOMAD::Double d;

    if (_n < 1)
    {
        std::string err = "Cannot read values into an empty array";
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }

    if (strDouble[0] == pStart)
    {
        if (strDouble.size() != _n+2)
        {
            std::string err = "Cannot read values:" + valueString;
            err += ", size is not compatible with array size:" + std::to_string(_n);
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        if (strDouble[_n+1] != NOMAD::ArrayOfDouble::pEnd)
        {
            std::string err = "Incompatible format for array: " + valueString;
            err += ". Must be delimited by \"" + NOMAD::ArrayOfDouble::pStart;
            err += "\" and \"" + NOMAD::ArrayOfDouble::pEnd + ")\" parenthesis.";
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        for (size_t i = 1; i <= _n; i++)
        {
            if ("-" == strDouble[i])
            {
                // Undefined value
                continue;
            }
            else
            {
                // Convert string to double. If it does not work, throw an error.
                if (!d.atof(strDouble[i]))
                {
                    std::string err = "Error reading array " + valueString + ". ";
                    err += "Cannot convert string " + strDouble[i] + " to double.";
                    throw NOMAD::Exception(__FILE__, __LINE__, err);
                }
            }
            _array[i-1] = d;
        }
    }

    // Fill array with all the same value.
    // Sample valueString:
    // * -6.0
    else if (strDouble[0] == "*")
    {
        if (strDouble.size() != 2)
        {
            std::cerr << "Warning: Expecting an array of size 2 for this parameter entry: " << valueString << std::endl;
        }

        if (!d.atof(strDouble[1]))
        {
            std::string err = "Error reading array " + valueString + ". ";
            err += "Cannot convert string " + strDouble[1] + " to double.";
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        // Set all array values to d.
        this->set(d);
    }

    // Fill array indices with the given value, or set it as toBeDefined().
    // Sample valueStrings:
    // 2 -6.0
    // 2
    // 2-4 3.42
    // 2-4
    else if (strDouble.size() <= 2)
    {
        // Process first value and put its indices in a vector.
        // Ex. 2 -> { 2 }
        // Ex. 2-4 -> { 2, 3, 4 }
        std::vector<int> indexRange;
        bool firstValueProcessed = false;
        if (d.atof(strDouble[0]) && d.isInteger() && d < (double)_n)
        {
            // First value is a valid index.
            indexRange.push_back(d.round());
            firstValueProcessed = true;
        }
        else
        {
            size_t hyphenIndex = strDouble[0].find_first_of("-");
            if (hyphenIndex != std::string::npos)
            {
                // First value is an index range.
                std::string firstIndexStr = strDouble[0].substr(0, hyphenIndex);
                std::string lastIndexStr  = strDouble[0].substr(hyphenIndex+1, strDouble[0].size());
                if (d.atof(firstIndexStr) && d.isInteger() && d < (double)_n)
                {
                    NOMAD::Double dLast;
                    if (dLast.atof(lastIndexStr) && dLast.isInteger() && dLast < (double)_n)
                    {
                        // Push indices until dLast is reached.
                        for (int index = d.round(); index <= dLast; index++)
                        {
                            indexRange.push_back(index);
                        }
                        firstValueProcessed = true;
                    }
                }
            }
        }
        if (!firstValueProcessed)
        {
            std::string err = "Error: cannot use " + strDouble[0] + " as an index.";
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        // Process second value, if available, or set array[index] to be defined for all index in indexRange.
        if (1 == strDouble.size())
        {
            // Set array[index] to be defined.
            for (auto index : indexRange)
            {
                _array[index].setToBeDefined();
            }
        }
        else if (d.atof(strDouble[1]))
        {
            // Set array[index] to d
            for (auto index : indexRange)
            {
                _array[index] = d;
            }
        }
        else
        {
            std::string err = "Cannot read value as an array of doubles: " + valueString;
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }
    }

    else
    {
        std::string err = "Cannot read value as an array of doubles: " + valueString;
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


NOMAD::ArrayOfDouble NOMAD::ArrayOfDouble::abs() const
{
    NOMAD::ArrayOfDouble absArray(_n);
    for (size_t i = 0; i < size(); i++)
    {
        // Undefined values are ignored.
        if (_array[i].isDefined())
        {
            absArray[i] = _array[i].abs();
        }
    }

    return absArray;
}


NOMAD::Double NOMAD::ArrayOfDouble::max() const
{
    NOMAD::Double max;
    for (size_t j = 0; j < size(); j++)
    {
        if (!_array[j].isDefined())
        {
            continue;
        }
        if (!max.isDefined() || _array[j] > max)
        {
            max = _array[j];
        }
    }

    return max;
}


/*----------------------------------------------------------*/
/*                    scalar multiplication                 */
/*----------------------------------------------------------*/
const NOMAD::ArrayOfDouble & NOMAD::ArrayOfDouble::operator *= ( const NOMAD::Double & d )
{
    NOMAD::Double * p = _array;
    for (size_t k = 0 ; k < _n ; ++k , ++p)
    {
        *p *= d;
    }

    return *this;
}


/*----------------------------------------------------------*/
/*                    scalar division                       */
/*----------------------------------------------------------*/
const NOMAD::ArrayOfDouble & NOMAD::ArrayOfDouble::operator /= ( const NOMAD::Double & d )
{
    NOMAD::Double * p = _array;
    for (size_t k = 0 ; k < _n ; ++k , ++p)
    {
        *p /= d;
    }

    return *this;
}

/*----------------------------------------------------------*/
/*                   addition of two arrays                 */
/*----------------------------------------------------------*/
const NOMAD::ArrayOfDouble NOMAD::ArrayOfDouble::operator+(const NOMAD::ArrayOfDouble& p) const
{
    if (p._n != _n)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "x + y: x.size != y.size" );
    }
    NOMAD::ArrayOfDouble          tmp ( _n );
    NOMAD::Double       * p1 = tmp._array;
    const NOMAD::Double * p2 =     _array;
    const NOMAD::Double * p3 =   p._array;

    for (size_t k = 0 ; k < _n ; ++k , ++p1 , ++p2 , ++p3)
    {
        *p1 = *p2 + *p3;
    }

    return tmp;
}


/*----------------------------------------------------------*/
/*                 substraction of two arrays               */
/*----------------------------------------------------------*/
const NOMAD::ArrayOfDouble NOMAD::ArrayOfDouble::operator-(const NOMAD::ArrayOfDouble& p) const
{
    if (p._n != _n)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "x - y: x.size != y.size" );
    }
    NOMAD::ArrayOfDouble          tmp ( _n );
    NOMAD::Double       * p1 = tmp._array;
    const NOMAD::Double * p2 =     _array;
    const NOMAD::Double * p3 =   p._array;

    for (size_t k = 0 ; k < _n ; ++k , ++p1 , ++p2 , ++p3)
    {
        *p1 = *p2 - *p3;
    }

    return tmp;
}

bool NOMAD::ArrayOfDouble::roundToPrecision(const NOMAD::ArrayOfDouble & precision)
{
    bool modif = false;

    for (size_t i = 0; i < _n; i++)
    {
        if (_array[i].roundToPrecision(precision[i]))
        {
            modif = true;
        }
    }
    return modif;

}


bool NOMAD::ArrayOfDouble::isMultipleOf(const NOMAD::ArrayOfDouble &granularity, int &index) const
{
    bool allMult = true;
    index = -1;

    for (size_t i = 0; i < _n; i++)
    {
        if (0.0 == granularity[i])
        {
            continue;
        }
        NOMAD::Double xi = _array[i];
        // Verify that:
        //     round (xi / granularity) * granularity ~= xi
        // i.e. that xi / granularity is (roughly) an integer.
        if (!xi.isMultipleOf(granularity[i]))
        {
            index = static_cast<int>(i);
            allMult = false;
            break;
        }
    }

    return allMult;
}


/*-----------------------------------------------------------*/
/*                     affectation operator                  */
/*-----------------------------------------------------------*/
NOMAD::ArrayOfDouble& NOMAD::ArrayOfDouble::operator= (const NOMAD::ArrayOfDouble &arrayOfDouble)
{
    if (this == &arrayOfDouble)
    {
        return *this;
    }

    if (_n != arrayOfDouble._n)
    {
        delete [] _array;
        _n = arrayOfDouble._n;
        if (_n > 0)
        {
            _array = new NOMAD::Double [_n];
        }
        else
        {
            _array = nullptr;
        }
    }

    NOMAD::Double       * array1 = _array;
    const NOMAD::Double * array2 = arrayOfDouble._array;
    for (size_t k = 0; k < _n; ++k, ++array1, ++array2)
    {
        *array1 = *array2;
    }

    return *this;
}


/*----------------------------------------------------------------------*/
/*  Set the ArrayOfDouble's value given by index with the Double d      */
/*  If relative==true, set the value relative to the bounds lb and ub   */
/*----------------------------------------------------------------------*/
void NOMAD::ArrayOfDouble::set(size_t index,
                               const NOMAD::Double &d,
                               bool relative,
                               const NOMAD::Double &lb,
                               const NOMAD::Double &ub)
{
    if (index >= size())
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Set: invalid index");
    }

    if (relative)
    {
        if (!lb.isDefined() || !ub.isDefined())
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Set: invalid bounds");
        }

        if (!d.isDefined() || d < 0.0 || d > 1.0)
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Set: invalid value to set coordinate (0<=d<1) relative to bounds" );
        }

        _array[index] = d * (ub - lb);
    }
    else
    {
        _array[index] = d;
    }

}


/*----------------------------------------------------------------------*/
/*  Set the ArrayOfDouble's values with the Double array 'a' of size n  */
/*  The ArrayOfDouble's dimension is changed to n.                      */
/*----------------------------------------------------------------------*/
void NOMAD::ArrayOfDouble::set(size_t n, const NOMAD::Double *a)
{
    if (n == 0 || !a)
        return;

    if (_n != n)
    {
        delete [] _array;
        _n      = n;
        _array = new NOMAD::Double [_n];
    }

    NOMAD::Double* array = _array;
    for (size_t k = 0; k < _n; ++k, ++array, ++a)
    {
        *array = *a;
    }
}


/*-----------------------------------------------------------*/
/*    Check if all the values of the array are defined       */
/*-----------------------------------------------------------*/
bool NOMAD::ArrayOfDouble::isComplete() const
{
    if (_n == 0)
    {
        return false;
    }

    const NOMAD::Double* d = _array;
    for (size_t i = 0; i < _n; ++i, ++d)
    {
        if (!d->isDefined())
        {
            return false;
        }
    }
    return true;
}


/*---------------------------------------------------------------*/
/*  Check if at least one value is defined in the array _array   */
/*---------------------------------------------------------------*/
bool NOMAD::ArrayOfDouble::isDefined() const
{
    if (_n == 0)
        return false;
    const NOMAD::Double * array = _array;
    for (size_t i = 0; i < _n; ++i, ++array)
        if (array->isDefined())
            return true;
    return false;
}


/*---------------------------------------------------------------*/
/*          Count the number of values that are defined          */
/*---------------------------------------------------------------*/
size_t NOMAD::ArrayOfDouble::nbDefined() const
{
    const NOMAD::Double *array = _array;
    size_t k = 0;
    for (size_t i = 0; i < _n; ++i, ++array)
    {
        if (array->isDefined())
        {
            ++k;
        }
    }
    return k;
}


/*---------------------------------------------------------------*/
/*      Return true if at least one value is to be defined       */
/*---------------------------------------------------------------*/
bool NOMAD::ArrayOfDouble::toBeDefined() const
{
    const NOMAD::Double *array = _array;
    for (size_t i = 0; i < _n; ++i, ++array)
    {
        if (array->toBeDefined())
        {
            return true;
        }
    }

    return false;
}


/*---------------------------------------------------*/
/* Throw an exception if 2 array sizes do not match. */
/*---------------------------------------------------*/
void NOMAD::ArrayOfDouble::verifySizesMatch(size_t n1, size_t n2,
                                            std::string filename,
                                            size_t linenum) const
{
    std::ostringstream oss;
    if (n1 != n2 || 0 == n2)
    {
        oss << "ArrayOfDouble comparison operator: Cannot compare arrays of different sizes (";
        oss << n1;
        oss << " and ";
        oss << n2;
        oss << ")";
        throw NOMAD::Exception(filename, linenum, oss.str());
    }
    else if (0 == n1 || 0 == n2)
    {
        oss << "ArrayOfDouble comparison operator: Empty array";
        throw NOMAD::Exception(filename, linenum, oss.str());
    }
}


/*----------------------------------------------------------*/
/*                           operator==                     */
/*----------------------------------------------------------*/
bool NOMAD::ArrayOfDouble::operator== (const NOMAD::ArrayOfDouble &arrayOfDouble) const
{
    if (this == &arrayOfDouble)
    {
        return true;
    }
    if (arrayOfDouble._n != _n)
    {
        return false;
    }

    const NOMAD::Double * array1 = _array;
    const NOMAD::Double * array2 = arrayOfDouble._array;
    for (size_t k = 0; k < _n; ++k, ++array1, ++array2)
    {
        if (!array1->isDefined() || !array2->isDefined() || *array1 != *array2)
        {
            return false;
        }
    }

    return true;
}


/*---------------------------------------------------------------------------------*/
/* Comparison operator '<=': it is used to verify bounds.                          */
/* Verify all values of this are inferior or equal to all values of arrayOfDouble. */
/*---------------------------------------------------------------------------------*/
bool NOMAD::ArrayOfDouble::operator<= (const NOMAD::ArrayOfDouble& arrayOfDouble) const
{
    bool isInferior = true;
    bool isStrictlyInferior = false;

    this->compare(arrayOfDouble, isInferior, isStrictlyInferior);

    return isInferior;
}


/*---------------------------------------------------------------------------------*/
/* Comparison operator '<': it is used to verify bounds.                           */
/* Verify all values of this are inferior or equal to all values of arrayOfDouble, */
/* and that at least one value is strictly inferior to the value for the           */
/* same index in arrayOfDouble.                                                    */
/* Note:                                                                           */
/* An operator< is defined for Point, but it throws an exception.                  */
/*---------------------------------------------------------------------------------*/
bool NOMAD::ArrayOfDouble::operator< (const NOMAD::ArrayOfDouble& arrayOfDouble) const
{
    bool isInferior = true;
    bool isStrictlyInferior = false;

    this->compare(arrayOfDouble, isInferior, isStrictlyInferior);

    return isInferior && isStrictlyInferior;
}


// Helper function for operator< and operator<=
void NOMAD::ArrayOfDouble::compare(const NOMAD::ArrayOfDouble& arrayOfDouble,
                                   bool &isInferior,
                                   bool &isStrictlyInferior) const
{
    verifySizesMatch(_n, arrayOfDouble._n, __FILE__, __LINE__);

    isInferior = true;
    isStrictlyInferior = false;
    for (size_t i = 0; isInferior && i < _n; i++)
    {
        if (!_array[i].isDefined() || !arrayOfDouble[i].isDefined())
        {
            throw NOMAD::Exception(__FILE__, __LINE__,
                                   "ArrayOfDouble comparison operator: Undefined value in array");
        }
        if (_array[i] < arrayOfDouble[i])
        {
            isStrictlyInferior = true;
        }
        else if (_array[i] > arrayOfDouble[i])
        {
            isInferior = false;
        }
    }

}


/*---------------*/
/*    display    */
/*---------------*/
std::string NOMAD::ArrayOfDouble::display(const NOMAD::ArrayOfDouble &precision) const
{
    std::ostringstream oss;
    // Fixed
    oss.setf(std::ios::fixed, std::ios::floatfield);

    for (size_t i = 0; i < size(); i++)
    {
        if (i > 0)
        {
            oss << " ";
        }

        int prec = -1;
        if (precision.isDefined() && precision[i].isDefined())
        {
            prec = static_cast<int>(precision[i].round());
        }
        oss << ((*this)[i]).display(prec);

    }

    return oss.str();
}
