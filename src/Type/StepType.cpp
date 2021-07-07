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
 \file   StepType.cpp
 \brief  types for Types for Steps that generate points (implementation)
 \author Viviane Rochon Montplaisir
 \date   June 2021
 \see    StepType.hpp
 */

#include "../Type/StepType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"


bool NOMAD::isAlgorithm(const StepType& stepType)
{
    switch (stepType)
    {
        case NOMAD::StepType::ALGORITHM_LH:
        case NOMAD::StepType::ALGORITHM_MADS:
        case NOMAD::StepType::ALGORITHM_NM:
        case NOMAD::StepType::ALGORITHM_PHASE_ONE:
        case NOMAD::StepType::ALGORITHM_PSD_MADS_SUBPROBLEM:
        case NOMAD::StepType::ALGORITHM_PSD_MADS:
        case NOMAD::StepType::ALGORITHM_SGTELIB_MODEL:
        case NOMAD::StepType::ALGORITHM_SSD_MADS:
        case NOMAD::StepType::ALGORITHM_QUAD_MODEL:
        case NOMAD::StepType::ALGORITHM_VNS_MADS:
            return true;
        default:
            return false;
    }
    return false;   // Only here to avoid compilation warnings
}


std::map<NOMAD::StepType, std::string>& NOMAD::dictStepType()
{
    static std::map<NOMAD::StepType,std::string> dictionary = {
        {NOMAD::StepType::ALGORITHM_LH, "Latin Hypercube"},
        {NOMAD::StepType::ALGORITHM_MADS, "MADS"},
        {NOMAD::StepType::ALGORITHM_NM, "Nelder-Mead"},
        {NOMAD::StepType::ALGORITHM_PHASE_ONE, "Phase One"},
        {NOMAD::StepType::ALGORITHM_PSD_MADS_SUBPROBLEM, "PSD-Mads subproblem"},
        {NOMAD::StepType::ALGORITHM_PSD_MADS, "PSD-Mads"},
        {NOMAD::StepType::ALGORITHM_SGTELIB_MODEL, "Sgtelib Model"},
        {NOMAD::StepType::ALGORITHM_SSD_MADS, "SSD-Mads"},
        {NOMAD::StepType::ALGORITHM_QUAD_MODEL, "Quad Model"},
        {NOMAD::StepType::ALGORITHM_VNS_MADS, "VNS Mads"},
        {NOMAD::StepType::INITIALIZATION, "Initialization"},
        {NOMAD::StepType::ITERATION, "Iteration"},
        {NOMAD::StepType::MAIN, "Main"},
        {NOMAD::StepType::MAIN_OBSERVE, "Observe"},
        {NOMAD::StepType::MAIN_SUGGEST, "Suggest"},
        {NOMAD::StepType::MEGA_ITERATION, "MegaIteration"},
        {NOMAD::StepType::MEGA_SEARCH_POLL, "MegaSearchPoll"},

        {NOMAD::StepType::NM_CONTINUE, "NM Continue"},
        {NOMAD::StepType::NM_EXPAND, "NM Expansion"},
        {NOMAD::StepType::NM_INITIAL, "NM Initial step type"},
        {NOMAD::StepType::NM_INITIALIZE_SIMPLEX, "NM Initialize Simplex"},
        {NOMAD::StepType::NM_INSERT_IN_Y, "NM Insert in Y"},
        {NOMAD::StepType::NM_INSIDE_CONTRACTION, "NM Inside Contraction"},
        {NOMAD::StepType::NM_OUTSIDE_CONTRACTION, "NM Outside Contraction"},
        {NOMAD::StepType::NM_REFLECT, "NM Reflect"},
        {NOMAD::StepType::NM_SHRINK, "NM Shrink"},
        {NOMAD::StepType::NM_UNSET, "NM step type not set"},

        {NOMAD::StepType::OPTIMIZE, "Optimize"},
        {NOMAD::StepType::POLL, "Poll"},
        {NOMAD::StepType::POLL_METHOD_DOUBLE, "Double Poll Method"},
        {NOMAD::StepType::POLL_METHOD_ORTHO_NPLUS1_NEG, "Ortho N+1 Neg Poll Method"},
        {NOMAD::StepType::POLL_METHOD_ORTHO_2N, "Ortho 2N Poll Method"},
        {NOMAD::StepType::POLL_METHOD_SINGLE, "Single Poll Method"},
        {NOMAD::StepType::POLL_METHOD_UNI_NPLUS1, "Uniform N+1 Poll Method"},
        {NOMAD::StepType::SEARCH, "Search"},
        {NOMAD::StepType::SEARCH_METHOD_LH, "Latin Hypercube Search Method"},
        {NOMAD::StepType::SEARCH_METHOD_NM, "Nelder-Mead Search Method"},
        {NOMAD::StepType::SEARCH_METHOD_QUAD_MODEL, "Quadratic Model Search Method"},
        {NOMAD::StepType::SEARCH_METHOD_SGTELIB_MODEL, "Sgtelib Model Search Method"},
        {NOMAD::StepType::SEARCH_METHOD_SPECULATIVE, "Speculative Search Method"},
        {NOMAD::StepType::SEARCH_METHOD_USER, "User-Defined Search Method"},
        {NOMAD::StepType::SEARCH_METHOD_VNS_MADS, "VNS Mads Search Method"},
        {NOMAD::StepType::SURROGATE_EVALUATION, "Points evaluated using static surrogate"},
        {NOMAD::StepType::TERMINATION, "Termination"},
        {NOMAD::StepType::UNDEFINED, "Undefined"},
        {NOMAD::StepType::UPDATE, "Update"}
    };

    return dictionary;
}


// Convert a NOMAD::StepType to a string for display.
std::string NOMAD::stepTypeToString(const NOMAD::StepType& stepType)
{
    std::map<NOMAD::StepType, std::string>::iterator it = dictStepType().find(stepType);
    return it->second;
}


// Convert a vector of StepTypes to a string for display.
std::string NOMAD::StepTypeListToString(const StepTypeList& stepTypeList)
{
    std::string s;
    bool showAlgorithm = false; // Too much information in general.
    if (stepTypeList.size() > 1 && NOMAD::StepType::INITIALIZATION == (*stepTypeList.begin()))
    {
        showAlgorithm = true;
    }

    bool first = true;  // First Step Type to be shown
    for (auto it = stepTypeList.rbegin(); it < stepTypeList.rend(); ++it)
    {
        if (!showAlgorithm && isAlgorithm(*it))
        {
            continue;
        }
        if (!first)
        {
            s += " - ";
        }
        s += NOMAD::stepTypeToString(*it);
        first = false;
    }
    return s;
}
