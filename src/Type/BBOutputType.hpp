/**
 \file   BBOutputType.hpp
 \brief  types for BBOutput
 \author Viviane Rochon Montplaisir
 \date   September 2018
 \see    BBOutput.hpp
 */
#ifndef __NOMAD400_BB_OUTPUT_TYPE__
#define __NOMAD400_BB_OUTPUT_TYPE__

#include <string>
#include <sstream>
#include <vector>

#include "../nomad_nsbegin.hpp"


/// Blackbox output types
enum class BBOutputType
{
    OBJ,        ///< Objective value
    EB,         ///< Extreme barrier constraint
    PB,         ///< Progressive barrier constraint
    CNT_EVAL,   ///< Output set to 0 or 1 to count the blackbox evaluation or not
    //STAT_AVG, ///< Stat (average)
    //STAT_SUM, ///< Stat (sum)
    BBO_UNDEFINED ///< Output ignored
};

/// Definition for the list of blackbox output types
typedef std::vector<BBOutputType> BBOutputTypeList;

typedef BBOutputTypeList::const_iterator BBOutputTypeListIt;


/// Utility for BBOutputType
/**
 Convert a string (ex "OBJ", "EB", "PB"...)
 to a BBOutputType.
 */
BBOutputType stringToBBOutputType(const std::string &s);

/// Utility for BBOutputType
/**
 Convert a string containing multiple BBOutputTypes (ex "OBJ EB PB PB")
 to a BBOutputTypeList.
 */
BBOutputTypeList stringToBBOutputTypeList(const std::string &s);

/// Utility for BBOutputType
/**
 Convert a BBOutputTypeList into a string
 */
std::string BBOutputTypeListToString ( const BBOutputTypeList & bbotList );

/// Helper to test if a BBOutputType is a constraint (PB, EB, ....)
bool BBOutputTypeIsConstraint(const BBOutputType & bbotType);

/// Count the number of constraints
size_t getNbConstraints(const BBOutputTypeList& bbotList);

/// Verify if the BBOutputType defines a constraint
bool isConstraint(const BBOutputType& bbot);

/// Count the number of objectives
size_t getNbObj(const BBOutputTypeList& bbotList);

/// Read and interpret BBOutputType
inline std::ostream& operator<<(std::ostream& os, const BBOutputType &bbot)
{
    switch (bbot)
    {
        case BBOutputType::OBJ:
            os << "OBJ";
            break;
        case BBOutputType::PB:
            os << "PB";
            break;
        case BBOutputType::EB:
            os << "EB";
            break;
        case BBOutputType::CNT_EVAL:
            os << "CNT_EVAL";
            break;
        case BBOutputType::BBO_UNDEFINED:
        default:
            return os << "BBO_UNDEFINED";
            break;
    }

    return os;
}

/// Display BBOutputType
inline std::ostream& operator<<(std::ostream& out, const BBOutputTypeList &bboutputtypelist)
{
    BBOutputTypeListIt it;
    bool first = true;
    for (it = bboutputtypelist.begin(); it != bboutputtypelist.end(); ++it)
    {
        if (!first)
        {
            out << " ";
        }
        out << *it;
        first = false;
    }
    return out;
}



#include "../nomad_nsend.hpp"

#endif // __NOMAD400_BB_OUTPUT_TYPE__
