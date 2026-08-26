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

#include <algorithm>    // For std::find, std::remove, std::sort

#include "../../Cache/CacheBase.hpp"
#include "../../Algos/Ads/AdsProgressiveBarrier.hpp"
#include "../../Math/Direction.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::AdsProgressiveBarrier::addCachePointTag(size_t tag)
{
    // Add tag only if not already present
    if (std::find(_cachePointTags.begin(), _cachePointTags.end(), tag) == _cachePointTags.end())
    {
        _cachePointTags.push_back(tag);
    }
}

void NOMAD::AdsProgressiveBarrier::addCachePointTags(const std::vector<size_t>& tags)
{
    for (size_t tag : tags)
    {
        addCachePointTag(tag);
    }
}

void NOMAD::AdsProgressiveBarrier::removeCachePointTag(size_t tag)
{
    _cachePointTags.erase(
        std::remove(_cachePointTags.begin(), _cachePointTags.end(), tag),
        _cachePointTags.end());
}

void NOMAD::AdsProgressiveBarrier::clearCachePointTags()
{
    _cachePointTags.clear();
}

bool NOMAD::AdsProgressiveBarrier::isInPuncturedSpace(const Point& trialPoint,
                                                      const Double& minDelta,
                                                      size_t iterationTag) const
{
    if (_cachePointTags.empty())
    {
        // No tags stored: treat point as in the punctured space
        // (default behavior when the list is empty)
        return true;
    }

    auto evalType = getEvalType();
    std::vector<EvalPoint> relevantPoints;
    
    // Collect all cache points whose tags are stored and that were successfully evaluated
    NOMAD::CacheBase::getInstance()->browse([&](const NOMAD::EvalPoint& cacheEvalPoint)
    {
        // Check whether the point's tag is in our list
        size_t pointTag = cacheEvalPoint.getTag();
        if (std::find(_cachePointTags.begin(), _cachePointTags.end(), pointTag) != _cachePointTags.end())
        {
            // Require successful evaluation
            if (cacheEvalPoint.getEvalStatus(evalType) == NOMAD::EvalStatusType::EVAL_OK
                && cacheEvalPoint.getTag() <= iterationTag)
            {
                relevantPoints.push_back(cacheEvalPoint);
            }
        }
    });

    if (relevantPoints.empty())
    {
        // No relevant points: treat as in the punctured space
        return true;
    }

    // Sort points by barrier ordering (best points first by f and h)
    auto computeType = getFHComputeType();
    std::sort(relevantPoints.begin(), relevantPoints.end(),
        [&](const EvalPoint& ep1, const EvalPoint& ep2)
        {
            auto eval1 = ep1.getEval(evalType);
            auto eval2 = ep2.getEval(evalType);
            
            if (nullptr == eval1 || nullptr == eval2)
            {
                return false;
            }
            
            bool feas1 = eval1->isFeasible(computeType.fhComputeTypeS);
            bool feas2 = eval2->isFeasible(computeType.fhComputeTypeS);
            
            // Feasible points always rank above infeasible ones
            if (feas1 != feas2)
            {
                return feas1;
            }
            
            // Both feasible: compare by f
            if (feas1)
            {
                return ep1.getF(computeType) < ep2.getF(computeType);
            }
            
            // Both infeasible: compare by h then f
            Double h1 = ep1.getH(computeType);
            Double h2 = ep2.getH(computeType);
            if (h1 != h2)
            {
                return h1 < h2;
            }
            return ep1.getF(computeType) < ep2.getF(computeType);
        });

    // Count how many points lie within distance minDelta of trialPoint
    size_t count = 0;
    for (const auto& cacheEvalPoint : relevantPoints)
    {
        // Distance between cache point and trial point
        NOMAD::ArrayOfDouble diff = (*cacheEvalPoint.getX() - trialPoint);
        
        // Use Direction to compute the norm
        NOMAD::Direction dir(diff);
        NOMAD::Double dist = dir.norm(NOMAD::NormType::L2);
        
        // Compare to precomputed minimum delta mesh size
        if (dist <= minDelta)
        {
            count++;
        }
    }
    
    // In punctured space if no close point found (count == 0).
    // A point is added to _cachePointTags only after passing the punctured check,
    // so the trial point itself is not in the list during this check.
    return count == 0;
}

