/**
 \file   utils.hpp
 \brief  Utility functions (headers)
 \author Sebastien Le Digabel, modified by Viviane Rochon Montplaisir
 \date   March 2017
 \see    utils.cpp
 */
#ifndef __NOMAD400_UTILS__
#define __NOMAD400_UTILS__

#include <list>

#include "../Util/defines.hpp"

#include "../nomad_nsbegin.hpp"

/// Transform an integer into a string.
/**
 \param i The integer -- \b IN.
 \return  The string.
 */
std::string itos ( const int i );


/// Transform a unsigned long (size_t) into a string.
/**
 \param i The unsigned long -- \b IN.
 \return  The string.
 */
std::string itos ( const size_t i );

/// Put a string into upper cases.
/**
 \param s The string -- \b IN/OUT.
 */
void toupper ( std::string & s );

/// Put a list of strings into upper cases.
/**
 \param ls The list of strings -- \b IN/OUT.
 */
void toupper  ( std::list<std::string> & ls );

/// Convert a string into an integer.
/**
 \param s The string  -- \b IN.
 \param i The integer -- \b OUT.
 \return  A boolean equal to \c true if the conversion was possible.
 */
bool atoi(const std::string &s, int &i);

/// Convert a character into an integer.
/**
 \param c The character -- \b IN.
 \param i The integer   -- \b OUT.
 \return  A boolean equal to \c true if the conversion was possible.
 */
bool atoi(const char c, int &i);

/// Convert a string into a size_t.
/**
 \param s   The string  -- \b IN.
 \param st  The size_t -- \b OUT.
 \return    A boolean equal to \c true if the conversion was possible.
 */
bool atost(const std::string &s, size_t &st);

/// Convert a string with format "i-j" into two size_t i and j.
/**

 \param  s The string              -- \b IN.
 \param  i The first index \c i  -- \b OUT.
 \param  j The second index \c j -- \b OUT.
 \param  check_order A boolean indicating if \c i and \c j are to be compared
 -- \b IN -- \b optional (default = \c true).
 \return A boolean equal to \c true if the conversion was possible.
 */
bool stringToIndexRange ( const std::string & s                  ,
                          size_t            & i                  ,
                          size_t            & j                  ,
                          bool                check_order = true   );


/// Convert a success type to a string.
std::string enumStr(SuccessType success);

/// Convert a string in {"YES","NO","Y","N","0","1","TRUE","FALSE"} to a boolean.
bool stringToBool(const std::string &string);

/// Convert a bool to "true" or "false"
std::string boolToString(bool boolean);

/// Return the number of decimals of a string representing a double.
std::size_t nbDecimals(const std::string& s);


/**
 Given a string s with precision prec after the decimal point,
return the suggested width and the string padding needed for alignment.
 */
void getFormat(const std::string &s,
               const size_t prec,
               size_t &width,
               size_t &spacePadding);

/**
 * Given a string starting with '%',
 * separate the actual formatting from the tag.
 * Ex. "%5.2fOBJ" -> "%5.2f", "OBJ"
 * Ex. "%dBBO"    -> "%d",    "BBO"
 * Ex. "%12BBE"   -> "%12f",  "BBE"
 * "TIME" returns false.
 * "%4.2" returns false.
 *
 \return \c true if a valid format was found, \c false otherwise.
 */
bool separateFormat(const std::string &s, std::string &format, std::string &tag);

bool validFormat(std::string &s);

// For OpenMP threads
int getThreadNum();


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_UTILS__
