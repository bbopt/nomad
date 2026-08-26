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


#include <algorithm>    // For std::merge and std::unique

#include "../../Algos/Ads/AdsIteration.hpp"
#include "../../Algos/Ads/AdsProgressiveBarrier.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Cache/CacheSet.hpp"
#include "../../Math/Direction.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/utils.hpp"

bool NOMAD::AdsIteration::runImp()
{
    bool iterationSuccess = false;
    bool searchSuccessOutsidePunctured = false;

    // Parameter Update is handled at the upper level - MegaIteration.
    if ( nullptr != _megasearchpoll
        && !_stopReasons->checkTerminate())
    {
        iterationSuccess = doMegaSearchPoll();
    }
    else
    {
        // 0. Get the tag of the last eval point created -> used below
        // for the punctured space at iteration k
        const size_t iterationTag = NOMAD::EvalPoint::getCurrentTag();
        
        // 1. Search
        if ( nullptr != _search && ! _stopReasons->checkTerminate() )
        {
            iterationSuccess = doSearch();
        }

        if ( nullptr != _search && ! _stopReasons->checkTerminate() )
        {

            if (iterationSuccess)
            {
                
                // Check if trial point that gave a success is in the punctured space.
                auto barrier = getMegaIterationBarrier();
                auto newBestFeas = barrier->getCurrentIncumbentFeas();
                auto newBestInf  = barrier->getCurrentIncumbentInf();
                
                auto deltaMeshSize = getMesh()->getdeltaMeshSize();
                
                // Compute minimum of deltaMeshSize once
                NOMAD::Double minDelta;
                for (size_t i = 0; i < deltaMeshSize.size(); i++)
                {
                    if (!deltaMeshSize[i].isDefined())
                        continue;
                    if (!minDelta.isDefined() || deltaMeshSize[i] < minDelta)
                    {
                        minDelta = deltaMeshSize[i];
                    }
                }
                
                bool newBestIsInPuncturedSpace = false;
                
                // Try to cast barrier to AdsProgressiveBarrier
                auto adsBarrier = std::dynamic_pointer_cast<NOMAD::AdsProgressiveBarrier>(barrier);
                
                if (adsBarrier != nullptr)
                {
                    // Use AdsProgressiveBarrier's isInPuncturedSpace method
                    // Separate if-s to manage a new feasible best and a new infeasible best
                    if (nullptr != newBestFeas)
                    {
                        newBestIsInPuncturedSpace = adsBarrier->isInPuncturedSpace(*newBestFeas->getX(),
                                                                                    minDelta,
                                                                                    iterationTag);
                    }
                    if (!newBestIsInPuncturedSpace && nullptr != newBestInf)
                    {
                        newBestIsInPuncturedSpace = adsBarrier->isInPuncturedSpace(*newBestInf->getX(),
                                                                                    minDelta,
                                                                                    iterationTag);
                    }
                    
                    // Add the tag of the evaluated point to the barrier's tag list
                    if (newBestIsInPuncturedSpace)
                    {
                        if (nullptr != newBestFeas)
                        {
                            adsBarrier->addCachePointTag(newBestFeas->getTag());
                        }
                        else if (nullptr != newBestInf)
                        {
                            adsBarrier->addCachePointTag(newBestInf->getTag());
                        }
                    }
                }
                else
                {
                    // Fallback to original method if barrier is not AdsProgressiveBarrier
                    std::vector<NOMAD::EvalPoint> evalPointList;
                    auto evalType = barrier->getEvalType();
                    
                    auto crit = [&](const NOMAD::Point& newBest, const NOMAD::EvalPoint& cacheEvalPoint)
                    {
                        if (cacheEvalPoint.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK
                            && cacheEvalPoint.getTag() <= iterationTag)
                            return false;
                        
                        // Difference between cache point and new best
                        NOMAD::ArrayOfDouble diff = (*cacheEvalPoint.getX() - newBest);
                        
                        // Use Direction to compute the norm
                        NOMAD::Direction dir(diff);
                        NOMAD::Double dist = dir.norm(NOMAD::NormType::L2);
                        
                        // Compare to precomputed minimum delta mesh size
                        return dist <= minDelta;
                    };
                    
                    // Separate if-s to manage a new feasible best and a new infeasible best
                    if (nullptr != newBestFeas)
                    {
                        // If most than one point is found, the newBest is too close to
                        // an existing point up to iteration k (evaluated successfully)
                        // Note: the new best is already in cache and it counts for 1!
                        // If 1 point found -> pt is in punctured space at iteration k
                        newBestIsInPuncturedSpace = (1 == NOMAD::CacheBase::getInstance()->find(*newBestFeas->getX(),
                                                                                                crit,evalPointList,
                                                                                                2)) ;
                    }
                    if (!newBestIsInPuncturedSpace && nullptr != newBestInf)
                    {
                        newBestIsInPuncturedSpace = (1 == NOMAD::CacheBase::getInstance()->find(*newBestInf->getX(),
                                                                                                crit,evalPointList,
                                                                                                2)) ;
                    }
                }
                
                if (newBestIsInPuncturedSpace)
                {
                    OUTPUT_INFO_START
                    AddOutputInfo("Search Successful. Enlarge Delta frame size.");
                    OUTPUT_INFO_END
                }
                else
                {
                    // ADS rule:
                    // Search can return success, but if this point is outside the punctured
                    // space we reframe this Search success as iteration failure and force Poll.
                    // If Poll also fails afterwards, this iteration is classified as REFRAMING.
                    iterationSuccess = false;
                    searchSuccessOutsidePunctured = true;
                    
                    // The new best is used as frame center for poll
                    // Update the barrier ref best
                    barrier->updateRefBests();
                    
                    OUTPUT_INFO_START
                    AddOutputInfo("Search step produced a dominating point that is not in punctured space (too close to a tagged point). Let's do a poll around this point.");
                    OUTPUT_INFO_END
                }
            }
            
            // Iteration success may have been updated after check punctured space inclusion
            if (!iterationSuccess)
            {
                iterationSuccess = doPoll();
                
                if (!iterationSuccess && nullptr != _extendedPoll && _extendedPoll->isEnabled())
                {
                    OUTPUT_INFO_START
                    AddOutputInfo("Poll unsuccessful. Let's do an extended poll.");
                    OUTPUT_INFO_END
                    _extendedPoll->start();
                    iterationSuccess = _extendedPoll->run();
                    _extendedPoll->end();                    
                }
                
                if (!iterationSuccess && searchSuccessOutsidePunctured)
                {
                    OUTPUT_INFO_START
                    AddOutputInfo("Iteration classified as REFRAMING: Search success outside punctured space, then Poll failed to find an improving point.");
                    OUTPUT_INFO_END
                }
                
                // If poll was successful, check if the new best point is in the punctured space
                if (iterationSuccess)
                {
                    // Check if trial point that gave a success is in the punctured space.
                    auto barrier = getMegaIterationBarrier();
                    auto newBestFeas = barrier->getCurrentIncumbentFeas();
                    auto newBestInf  = barrier->getCurrentIncumbentInf();
                    
                    auto deltaMeshSize = getMesh()->getdeltaMeshSize();
                    
                    // Compute minimum of deltaMeshSize once
                    NOMAD::Double minDelta;
                    for (size_t i = 0; i < deltaMeshSize.size(); i++)
                    {
                        if (!deltaMeshSize[i].isDefined())
                            continue;
                        if (!minDelta.isDefined() || deltaMeshSize[i] < minDelta)
                        {
                            minDelta = deltaMeshSize[i];
                        }
                    }
                    
                    bool newBestIsInPuncturedSpace = false;
                    
                    // Try to cast barrier to AdsProgressiveBarrier
                    auto adsBarrier = std::dynamic_pointer_cast<NOMAD::AdsProgressiveBarrier>(barrier);
                    
                    if (adsBarrier != nullptr)
                    {
                        // Use AdsProgressiveBarrier's isInPuncturedSpace method
                        // Separate if-s to manage a new feasible best and a new infeasible best
                        if (nullptr != newBestFeas)
                        {
                            newBestIsInPuncturedSpace = adsBarrier->isInPuncturedSpace(*newBestFeas->getX(),
                                                                                        minDelta,
                                                                                        iterationTag);
                        }
                        if (!newBestIsInPuncturedSpace && nullptr != newBestInf)
                        {
                            newBestIsInPuncturedSpace = adsBarrier->isInPuncturedSpace(*newBestInf->getX(),
                                                                                        minDelta,
                                                                                        iterationTag);
                        }
                        
                        // Add the tag of the evaluated point to the barrier's tag list
                        if (newBestIsInPuncturedSpace)
                        {
                            if (nullptr != newBestFeas)
                            {
                                adsBarrier->addCachePointTag(newBestFeas->getTag());
                            }
                            else if (nullptr != newBestInf)
                            {
                                adsBarrier->addCachePointTag(newBestInf->getTag());
                            }
                        }
                    }
                }
            }
        }
    }


    // End of the iteration: iterationSuccess is true iff we have a full success.
    return iterationSuccess;
}
