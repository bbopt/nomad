/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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

#include "../../Cache/CacheBase.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/COOPMads/CacheSearchMethod.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Algos/EvcInterface.hpp"

void NOMAD::CacheSearchMethod::init()
{
    
    setStepType(NOMAD::StepType::SEARCH_METHOD_CACHE);
    
    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    if (evc->getCurrentEvalType() == NOMAD::EvalType::MODEL)
    {
        setEnabled(false);
        return;
    }
    if (!evc->getUseCache())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For COOP-Mads cache search, we need a cache.");
    }
    
    const bool isEnabled = getRunParams()->getAttributeValue<bool>("COOP_MADS_OPTIMIZATION_CACHE_SEARCH");
    setEnabled(isEnabled);
    
}


void NOMAD::CacheSearchMethod::generateTrialPointsFinal()
{
    
    // Get the barrier.
    if (nullptr == _megaIterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For COOP-Mads cache search, we need a MegaIteration among the parents.");
    }
    auto barrier = _megaIterAncestor->getBarrier();
    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For COOP-Mads cache search, we need a MegaIteration among the parents with a barrier.");
    }
    auto computeType = barrier->getFHComputeType();
    
    // Look in the cache for the best feasible points
    NOMAD::CacheInterface cacheInterface(this);
    std::vector<NOMAD::EvalPoint> evalPointList;
    cacheInterface.findBestFeas(evalPointList, computeType);
    for (auto & ep: evalPointList)
    {
        // Test best points on barrier
        if ( NOMAD::SuccessType::FULL_SUCCESS == barrier->getSuccessTypeOfPoints(std::make_shared<EvalPoint>(ep), nullptr))
        {
            OUTPUT_INFO_START
            std::string s = "Cache search found a point in cache that dominates the best feasible point in barrier: ";
            s += ep.display(computeType.Short()) ;
            AddOutputInfo(s);
            OUTPUT_INFO_END
            
            // Update point from using current feasible incumbent
            NOMAD::Point fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(NOMAD::IterationUtils::_parentStep);
            
            ep.setPointFrom(barrier->getCurrentIncumbentFeas(), fixedVariable);
            
            insertTrialPoint(ep);
        }
    }
    // Look in the cache for the best infeasible points
    cacheInterface.findBestInf(evalPointList, barrier->getHMax(), computeType);
    for (auto & ep: evalPointList)
    {
        if ( NOMAD::SuccessType::FULL_SUCCESS == barrier->getSuccessTypeOfPoints(std::make_shared<EvalPoint>(ep), nullptr))
        {
            OUTPUT_INFO_START
            std::string s = "Cache search found a point in cache that dominates one of the best infeasible point in barrier: ";
            s += ep.display(computeType.Short()) ;
            AddOutputInfo(s);
            OUTPUT_INFO_END
            
            // Update point from using current feasible incumbent
            NOMAD::Point fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(NOMAD::IterationUtils::_parentStep);
            
            ep.setPointFrom(barrier->getCurrentIncumbentInf(), fixedVariable);
            
            insertTrialPoint(ep);
        }
    }
    
}

bool NOMAD::CacheSearchMethod::evalTrialPoints(const NOMAD::Step *step,
                                            const size_t keepN,
                                            NOMAD::StepType removeStepType)
{
    if (_trialPoints.size()==0)
    {
        return false;
    }
    
    // Trial points should not be re-evaluated because
    // they are already in cache. We consider it is a success.
    NOMAD::EvcInterface evcInterface(step);
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    if (!evc->getUseCache())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"Must use cache to determine if trial points are success");
    }
    
    if ( evcInterface.countPointsThatNeedEval(_trialPoints) != 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"No trial point should need evaluation.");
    }
    
    // After generation, the trial points have no eval avail.
    // Let's update the trial points with evaluation from cache.
    NOMAD::CacheInterface cacheInterface(this);
    NOMAD::EvalPointSet evalPointSet;
    for (const auto & tp: _trialPoints)
    {
        NOMAD::EvalPoint ep;
        cacheInterface.find(*tp.getX(), ep);
        evalPointSet.insert(ep);
    }
    _trialPoints.clear();
    _trialPoints = evalPointSet;
    
    
    _trialPointsSuccess = NOMAD::SuccessType::FULL_SUCCESS;
    
    // Propagate trial points success type to generating method step (for example, poll and search)
    setSuccessType(_trialPointsSuccess);
    
    // Update step success stats from evc success stats
    updateStepSuccessStats(this);
    
    return true /* trial points are always considered as better */;
}

