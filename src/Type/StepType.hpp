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
 \file   StepType.hpp
 \brief  Types for Steps that generate points: Search and Poll methods, and Algorithms
 \author Viviane Rochon Montplaisir
 \date   June 2021
 \see    StepType.cpp
 */

#ifndef __NOMAD_4_0_STEP_TYPE__
#define __NOMAD_4_0_STEP_TYPE__

#include <map>
#include <sstream>
#include <vector>

#include "../nomad_nsbegin.hpp"

enum class StepType
{
    ALGORITHM_LH,               ///< Algorithm Latin Hypercube
    ALGORITHM_MADS,             ///< Algorithm Mads
    ALGORITHM_NM,               ///< Algorithm Nelder-Mead
    ALGORITHM_PHASE_ONE,        ///< Phase One
    ALGORITHM_PSD_MADS_SUBPROBLEM, ///< Subproblem in PSD-Mads
    ALGORITHM_PSD_MADS,         ///< Algorithm PSD-Mads
    ALGORITHM_QUAD_MODEL,       ///< Algorithm Quad Model
    ALGORITHM_SGTELIB_MODEL,    ///< Algorithm Quad Model
    ALGORITHM_SSD_MADS,         ///< Algorithm SSD-Mads
    ALGORITHM_VNS_MADS,         ///< Algorithm VNS-Mads
    INITIALIZATION,             ///< Initialization step
    ITERATION,                  ///< Iteration step
    MAIN,                       ///< Main step
    MAIN_OBSERVE,               ///< Main step for Observe
    MAIN_SUGGEST,               ///< Main step for Suggest
    MEGA_ITERATION,             ///< MegaIteration step
    MEGA_SEARCH_POLL,           ///< MegaSearchPoll

    NM_CONTINUE,                ///< NM continue
    NM_EXPAND,                  ///< NM Expansion
    NM_INITIAL,                 ///< NM initial step type
    NM_INITIALIZE_SIMPLEX,      ///< NM initialize simplex
    NM_INSERT_IN_Y,             ///< NM insert in Y
    NM_INSIDE_CONTRACTION,      ///< NM Inside Contraction
    NM_OUTSIDE_CONTRACTION,     ///< NM Outside Contraction
    NM_REFLECT,                 ///< NM Reflect
    NM_SHRINK,                  ///< NM Shrink
    NM_UNSET,                   ///< NM step type not set

    OPTIMIZE,                   ///< Sub-optimization
    POLL,                       ///< Poll
    POLL_METHOD_DOUBLE,         ///< Double poll method
    POLL_METHOD_ORTHO_NPLUS1_NEG, ///< Ortho N+1 neg poll method
    POLL_METHOD_ORTHO_2N,       ///< Ortho 2N poll method
    POLL_METHOD_SINGLE,         ///< Single poll method
    POLL_METHOD_UNI_NPLUS1,     ///< Uniform N+1 poll method
    SEARCH,                     ///< Search
    SEARCH_METHOD_LH,           ///< Latin hypercube search method
    SEARCH_METHOD_NM,           ///< Nelder-Mead search method
    SEARCH_METHOD_QUAD_MODEL,   ///< Quadratic model search method
    SEARCH_METHOD_SGTELIB_MODEL,///< Sgtelib model search method
    SEARCH_METHOD_SPECULATIVE,  ///< Speculative search method
    SEARCH_METHOD_USER,         ///< User-defined search method
    SEARCH_METHOD_VNS_MADS,     ///< VNS Mads search method
    SURROGATE_EVALUATION,       ///< Evaluating trial points using static surrogate
    TERMINATION,                ///< Termination
    UNDEFINED,                  ///< Unknown value (default)
    UPDATE                      ///< Update step
};


/// Definition for a vector of StepTypes
typedef std::vector<StepType> StepTypeList;

/// Helper to test if a StepType represents an Algorithm (ALGORITHM_MADS, etc).
bool isAlgorithm(const StepType& stepType);

std::map<StepType, std::string>& dictStepType();

// Convert an StepType to a string
std::string stepTypeToString(const StepType& stepType);

// Convert a StepTypeList to a string; show only pertinent information.
std::string StepTypeListToString(const StepTypeList& stepTypeList);


inline std::ostream& operator<<(std::ostream& out, const StepType &stepType)
{
    out << stepTypeToString(stepType);
    return out;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD_4_0_STEP_TYPE__
