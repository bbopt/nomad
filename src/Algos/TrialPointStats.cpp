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
\file   TrialPointStats.cpp
\brief  Manage stats for point generation methods (ex. Search and Poll)
\author Christophe Tribes
\date   August 2021
*/

#include "../Algos/Algorithm.hpp"
#include "../Algos/IterationUtils.hpp"
#include "../Algos/TrialPointStats.hpp"

#ifdef _OPENMP
omp_lock_t NOMAD::TrialPointStats::_updateLock  ;
#endif

void NOMAD::TrialPointStats::init()
{    
    _nbCalls = 0;
    
    initializeMap(_nbTotalEvalsDone);
    initializeMap(_nbCurrentEvalsDone);
    
    initializeMap(_nbTotalTrialPointsGenerated);
    initializeMap(_nbCurrentTrialPointsGenerated);
    
}

void NOMAD::TrialPointStats::initializeMap(std::map<EvalType, size_t> & counter)
{
    counter.clear();
    for (auto et: _allEvalType )
    {
        counter.insert(std::make_pair(et, 0));
    }
        
}

void NOMAD::TrialPointStats::incrementEvalsDone(size_t nb, EvalType evalType)
{
    _nbTotalEvalsDone.at(evalType) += nb;
    _nbCurrentEvalsDone.at(evalType) += nb;
}

void NOMAD::TrialPointStats::incrementTrialPointsGenerated(size_t nb, EvalType evalType)
{
    _nbTotalTrialPointsGenerated.at(evalType) += nb;
    _nbCurrentTrialPointsGenerated.at(evalType) += nb;
}

size_t NOMAD::TrialPointStats::getNbTrialPointsGenerated(EvalType evalType, bool totalCount) const
{
    if (totalCount)
        return size_t(_nbTotalTrialPointsGenerated.at(evalType));
    else
        return size_t(_nbCurrentTrialPointsGenerated.at(evalType));
}

size_t NOMAD::TrialPointStats::getNbEvalsDone(EvalType evalType, bool totalCount) const
{
    if (totalCount)
        return size_t(_nbTotalEvalsDone.at(evalType));
    else
        return size_t(_nbCurrentEvalsDone.at(evalType));
}

void NOMAD::TrialPointStats::updateWithCurrentStats(const TrialPointStats &trialPointStats)
{
    // Use the CURRENT stats of the given trialPointStats to update this current trialPointStats (CURRENT and TOTAL)
    for (auto et: _allEvalType )
    {
        _nbTotalEvalsDone.at(et) += trialPointStats.getNbEvalsDone(et, false);
        _nbCurrentEvalsDone.at(et) += trialPointStats.getNbEvalsDone(et, false);
        
        _nbTotalTrialPointsGenerated.at(et) += trialPointStats.getNbTrialPointsGenerated(et, false);
        _nbCurrentTrialPointsGenerated.at(et) += trialPointStats.getNbTrialPointsGenerated(et, false);
    }
}

void NOMAD::TrialPointStats::resetCurrentStats()
{
    for (auto et: _allEvalType )
    {
        _nbCurrentEvalsDone[et] = 0 ;
        _nbCurrentTrialPointsGenerated[et] = 0;
    }
    
}

void NOMAD::TrialPointStats::updateParentStats()
{
    // First try to update iteration utils parent
    // The parent can be an IterationUtils using an Algorithm to generate and evaluate trial point.
    // For example, VNS Search Method (IU) runs a VNS (Algo) which runs a Mads (Algo), etc. We need to pass the stats from Mads to VNS and from VNS to VNS Search Method.

    Step* step = const_cast<Step*>(_parentStep);
    while (nullptr != step)
    {
        if (nullptr != dynamic_cast<NOMAD::IterationUtils*>(step))
        {
            auto iu = dynamic_cast<NOMAD::IterationUtils*>(step);
            #ifdef _OPENMP
			    omp_init_lock(&_updateLock);
                omp_set_lock(&_updateLock);
            #endif // _OPENMP
                iu->updateStats(*this);
            #ifdef _OPENMP
                omp_unset_lock(&_updateLock);
            #endif // _OPENMP
            break;
        }
        else if (nullptr != dynamic_cast<NOMAD::Algorithm*>(step))
        {
            auto algo = dynamic_cast<NOMAD::Algorithm*>(step);
            #ifdef _OPENMP
			    omp_init_lock(&_updateLock);
                omp_set_lock(&_updateLock);
            #endif // _OPENMP
                algo->updateStats(*this);
            #ifdef _OPENMP
                omp_unset_lock(&_updateLock);
            #endif // _OPENMP
            break;
        }
        step = const_cast<Step*>(step->getParentStep());
    }
}

/*---------------*/
/*    display    */
/*---------------*/
std::string NOMAD::TrialPointStats::display() const
{
    std::string s;
    throw NOMAD::Exception(__FILE__, __LINE__,"Not yet implemented ");
    return s;
}

std::ostream& operator<<(std::ostream& os, NOMAD::TrialPointStats& stats)
{
    
    std::ostringstream oss;
    oss << stats.display();
    return os;
}
