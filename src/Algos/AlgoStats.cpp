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
\file   AlgoStats.cpp
\brief  Manage stats for Algos
\author Christophe Tribes
\date   August 2021
*/
#include "../Algos/AlgoStats.hpp"

void NOMAD::AlgoStats::startCounting()
{
    _nbEvalsDoneBB = _nbEvalsDoneModels = _nbEvalsDoneSurrogate = 0;

    _nbTrialPointsGeneratedBB = _nbTrialPointsGeneratedModels = _nbTrialPointsGeneratedSurrogate = 0;
    
    _nbGenerateCalls = 0;
    
#ifdef TIME_STATS
    _totalRealAlgoTime = 0.0;
    _startTime = NOMAD::Clock::getCPUTime();
    _totalCPUAlgoTime = 0.0;
#endif // TIME_STATS
    
}

void NOMAD::AlgoStats::endCounting()
{
#ifdef TIME_STATS
    _totalRealAlgoTime = NOMAD::Clock::getTimeSinceStart();
    _totalCPUAlgoTime += NOMAD::Clock::getCPUTime() - _startTime;
#endif // TIME_STATS
}


void NOMAD::AlgoStats::incrementEvalsDone(size_t nb, EvalType evalType)
{
    switch (evalType)
    {
        case NOMAD::EvalType::BB:
            _nbEvalsDoneBB += nb;
            break;
        case NOMAD::EvalType::SURROGATE:
            _nbEvalsDoneSurrogate += nb;
            break;
        case NOMAD::EvalType::MODEL:
            _nbEvalsDoneModels += nb;
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__,"Method stats increment evals done not implemented for " + evalTypeToString(evalType));
            break;
    }
        
}

void NOMAD::AlgoStats::incrementTrialPointsGenerated(size_t nb, EvalType evalType)
{
    switch (evalType)
    {
        case NOMAD::EvalType::BB:
            _nbTrialPointsGeneratedBB += nb;
            break;
        case NOMAD::EvalType::SURROGATE:
            _nbTrialPointsGeneratedSurrogate += nb;
            break;
        case NOMAD::EvalType::MODEL:
            _nbTrialPointsGeneratedModels += nb;
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__,"Algo stats increment trial points generated not implemented for " + evalTypeToString(evalType));
            break;
    }
        
}

size_t NOMAD::AlgoStats::getNbTrialPointsGenerated(EvalType evalType) const
{
    switch (evalType)
    {
        case NOMAD::EvalType::BB:
            return _nbTrialPointsGeneratedBB;
            break;
        case NOMAD::EvalType::SURROGATE:
            return _nbTrialPointsGeneratedSurrogate;
            break;
        case NOMAD::EvalType::MODEL:
            return _nbTrialPointsGeneratedModels;
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__,"Method get nb trial points generated not implemented for " + evalTypeToString(evalType));
            break;
    }
        
}

size_t NOMAD::AlgoStats::getNbEvalsDone(EvalType evalType) const
{
    switch (evalType)
    {
        case NOMAD::EvalType::BB:
            return _nbEvalsDoneBB;
            break;
        case NOMAD::EvalType::SURROGATE:
            return _nbEvalsDoneSurrogate;
            break;
        case NOMAD::EvalType::MODEL:
            return _nbEvalsDoneModels;
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__,"Method get nb evals done not implemented for " + evalTypeToString(evalType));
            break;
    }
        
}


/*---------------*/
/*    display    */
/*---------------*/
std::string NOMAD::AlgoStats::display() const
{
    std::string s;
    throw NOMAD::Exception(__FILE__, __LINE__,"Not yet implemented ");
    return s;
}

std::ostream& operator<<(std::ostream& os, NOMAD::AlgoStats& stats)
{
    
    std::ostringstream oss;
    oss << stats.display();
    return os;
}
