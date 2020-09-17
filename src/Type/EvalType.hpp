/**
 \file   EvalType.hpp
 \brief  Types for Evaluation: Blackbox, Surrogate
 \author Viviane Rochon Montplaisir
 \date   November 2019
 \see    EvalType.cpp
 */

#ifndef __NOMAD400_EVAL_TYPE__
#define __NOMAD400_EVAL_TYPE__

#include <sstream>

#include "../nomad_nsbegin.hpp"

// Evaluator type
enum class EvalType
{
    BB,                 ///< The evaluator is a blackbox.
    SGTE,               ///< The evaluator is a surrogate function,
                        /// potentially much faster to run than a blackbox.
    UNDEFINED           ///< Undefined: This value may be used when the
                        ///< EvalType is not mandatory
};


// Convert a string (ex "BB", "SGTE")
// to an EvalType.
EvalType stringToEvalType(const std::string &s);

// Convert an EvalType to a string
std::string evalTypeToString (const EvalType& evalType);


inline std::ostream& operator<<(std::ostream& out, const EvalType &evalType)
{
    out << evalTypeToString(evalType);
    return out;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD400_EVAL_TYPE__
