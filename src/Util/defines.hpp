/**
 \file   defines.hpp
 \brief  Definitions
 \author Sebastien Le Digabel, modified by Viviane Rochon Montplaisir
 \date   March 2017
 */
#ifndef __NOMAD400_DEFINES__
#define __NOMAD400_DEFINES__

#include <string>
#include <iostream>
#include <sstream>
#include <limits.h> // For INT_MAX
#include <limits>   // For numeric_limits
#include <cstdlib>
#include <memory>   // For shared_ptr, unique_ptr


// Define in order to display debug information
//#define DEBUG


// CASE Linux using gnu compiler
#ifdef __gnu_linux__
#define GCC_X
#endif

// CASE OSX using gnu compiler
#ifdef __APPLE__
#ifdef __GNUC__
#define GCC_X
#endif
#endif

// CASE Visual Studio C++ compiler
#ifdef _MSC_VER
#define WINDOWS
#pragma warning(disable:4996)
#endif

// For NOMAD random number generator
#if !defined(UINT32_MAX)
typedef unsigned int uint32_t;
#define UINT32_MAX 0xffffffff
#endif

#include "../nomad_nsbegin.hpp"

// Directory separator
#ifdef WINDOWS
const char        DIR_SEP = '\\';           ///< Directory separator
#else
const char        DIR_SEP = '/';            ///< Directory separator
#endif

/// Maximum number of variables.
const int MAX_DIMENSION = 1000;

/// Default epsilon used by Double
/** Use Parameters::set_EPSILON(), or parameter EPSILON,
 or Double::setEpsilon() to change it
 */
const double DEFAULT_EPSILON = 1e-13;

/// Default infinity string used by Double
/** Use Parameters::set_INF_STR(), or parameter INF_STR,
 or Double::setInfStr() to change it
 */
const std::string DEFAULT_INF_STR = "inf";

/// Default undefined value string used by Double
/** Use Parameters::set_UNDEF_STR(), or parameter UNDEF_STR,
 or Double::set_undefStr() to change it
 */
const std::string DEFAULT_UNDEF_STR = "NaN";
// Other strings recognized as NaN
const std::string DEFAULT_UNDEF_STR_HYPHEN = "-";
const std::string DEFAULT_UNDEF_STR_1 = "nan";

const double INF = std::numeric_limits<double>::max(); ///< Infinity
const double NaN = std::numeric_limits<double>::quiet_NaN(); ///< Quiet Not-A-Number
const int P_INF_INT = std::numeric_limits<int>::max(); ///< plus infinity for int
const int M_INF_INT = std::numeric_limits<int>::min(); ///< minus infinity for int
const size_t INF_SIZE_T = std::numeric_limits<size_t>::max();///< The infinity for \c size_t
const size_t INF_SHORT = std::numeric_limits<short>::max();///< The infinity for \c short

const double D_INT_MAX = UINT32_MAX; ///< The UINT32_MAX constant as a \c double


// Display precisions.
const int DISPLAY_PRECISION_STD = 6;  ///< Precision after decimal point (number of digits)
const int DISPLAY_PRECISION_FULL = 20;  ///< Display all decimals
const int NB_DIGITS_BEFORE_POINT = 3;   // "Precision" before decimal point
const int INT_DISPLAY_WIDTH = 3;        // Width for integers


// -------------------------
// Related to MADS algorithm
// -------------------------

/// Success type of an iteration.
//  Order is important.
enum class SuccessType
{
    NOT_EVALUATED,      ///< Not evaluated yet
    UNSUCCESSFUL,       ///< Failure
    PARTIAL_SUCCESS,    ///< Partial success (improving). Found an infeasible
    ///< solution with a better h. f is worse.
    FULL_SUCCESS        ///< Full success (dominating)
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_DEFINES__
