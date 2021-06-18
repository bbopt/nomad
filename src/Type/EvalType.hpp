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
 \file   EvalType.hpp
 \brief  Types for Evaluation: Blackbox, Model
 \author Viviane Rochon Montplaisir
 \date   November 2019
 \see    EvalType.cpp
 */

#ifndef __NOMAD_4_0_EVAL_TYPE__
#define __NOMAD_4_0_EVAL_TYPE__

#include <sstream>

#include "../nomad_nsbegin.hpp"

// Evaluator type
enum class EvalType
{
    BB,                 ///< The evaluator is a blackbox.
    MODEL,              ///< The evaluator is a model function,
                        /// potentially much faster to run than a blackbox.
    SURROGATE,          ///< The evaluator is a static surrogate,
                        /// potentially much faster to run than a blackbox.
    LAST,               ///< For iterations; Note: UNDEFINED evals are ignored.
    UNDEFINED           ///< Undefined: This value may be used when the
                        ///< EvalType is not mandatory
};


// Convert a string (ex "BB", "MODEL")
// to an EvalType.
EvalType stringToEvalType(const std::string &s);

// Convert an EvalType to a string
std::string evalTypeToString(const EvalType& evalType);


inline std::ostream& operator<<(std::ostream& out, const EvalType &evalType)
{
    out << evalTypeToString(evalType);
    return out;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD_4_0_EVAL_TYPE__
