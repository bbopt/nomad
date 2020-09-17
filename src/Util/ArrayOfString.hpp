/**
 \file   ArrayOfString.hpp
 \brief  Class to be used for an array of n strings
 \author Viviane Rochon Montplaisir
 \date   February 2018
 \see    ArrayOfString.cpp
 */

#ifndef __NOMAD400_ARRAYOFSTRING__
#define __NOMAD400_ARRAYOFSTRING__

#include <string>
#include <vector>

#include "../nomad_nsbegin.hpp"

/// Class for the representation of an array of n strings.
class ArrayOfString {
protected:
    /*---------*/
    /* Members */
    /*---------*/
    std::vector<std::string> _array;    ///< Each element of the array can be put in a string

public:

    /// Constructor #1
    /**
     \param n           Dimension of the array -- \b IN.
     \param initString  Initial string for all elements of the array -- \b IN.
     */
    explicit ArrayOfString(const size_t n = 0 ,
                           const std::string &initString = std::string());

    /// Constructor #2
    /**
     \param inputString     String to fill the array.
     \param separators      All characters that can be used as separator between the strings -- \b IN.
     */
    explicit ArrayOfString(const std::string & inputString,
                           const std::string & separators = " ");

    /// Destructor.
    virtual ~ArrayOfString();

    /*---------------*/
    /* Class methods */
    /*---------------*/
    /// Operator \c [].
    /**
     \param i   The index (0 for the first element) -- \b IN.
     \return    The \c (i+1)th string in the array.
     */
    const std::string& operator[](size_t i) const;

    /// Access to the dimension of the array.
    /**
     \return The dimension of the array.
     */
    size_t size() const { return _array.size(); }

    /// Check if empty
    /**
     \return \c true if the array is empty, \c false if not.
     */
    bool empty() const { return _array.empty(); }

    /// Add a string to the array
    /**
     \param s   The string to be add -- \b IN.
     */
    void add(const std::string& s) { _array.push_back(s); }

    /// Clear the array
    void clear() { _array.clear(); }

    /// Replace string at a given position.
    /**
     * Throw an error if the index is too large.
     \param index   Position -- \b IN.
     \param s       Replacement string -- \b IN.
     */
    void replace(const size_t index, const std::string& s);

    /// Erase the string which is at position index in the array.
    /**
     \param index   Position -- \b IN.
     \return        \c true if the string was removed from the array.
     */
    bool erase(const size_t index);

    /// Return first index of string s in array.
    /**
     \param s   String to find -- \b IN.
     \return    Position index if s is found and -1 if not found.
     */
    int find(const std::string& s) const;

    /*------------*/
    /* Comparison */
    /*------------*/
    /// Comparison operator \c ==.
    /**
     \param array   The right-hand side object -- \b IN.
     \return        \c true if  \c *this \c == \c array.
     */
    bool operator== (const ArrayOfString &array) const;

    /// Comparison operator \c !=.
    /**
     \param array    The right-hand side object -- \b IN.
     \return         \c true if  \c *this \c != \c array.
     */
    bool operator!= (const ArrayOfString &array) const { return !(*this == array); }

    /*-----------*/
    /*  Display  */
    /*-----------*/
    /// Display
    /**
    \return A string combining all elements of the array.
     */
    std::string display() const;

private:
    /// Helper method for splitting string
    std::vector<std::string> const splitString(const std::string & inputString,
                                               const std::string & separators = " ");


};

std::ostream& operator<< (std::ostream& out,
                          const ArrayOfString& arrayOfString);

#include "../nomad_nsend.hpp"
#endif // __NOMAD400_ARRAYOFSTRING__
