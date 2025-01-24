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
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/DMultiMads/DMultiMadsIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMegaIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMiddlePointSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"

void NOMAD::DMultiMadsMiddlePointSearchMethod::init()
{
    if (nullptr == getParentOfType<DMultiMadsIteration*>())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMadsMiddlePointSearchMethod only works for DMultiMads");
    }

    setStepType(NOMAD::StepType::SEARCH_METHOD_DMULTIMADS_MIDDLEPOINT);

    const bool isEnabled = _runParams->getAttributeValue<bool>("DMULTIMADS_MIDDLEPOINT_SEARCH");
    setEnabled(isEnabled);
}

void NOMAD::DMultiMadsMiddlePointSearchMethod::preRunValidations()
{
    std::shared_ptr<NOMAD::BarrierBase> barrier = nullptr;

    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    const auto refComputeType = evc->getComputeType();

    if (NOMAD::ComputeType::STANDARD != refComputeType)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do Middle Point Search for DMultiMads on ComputeType other than STANDARD.");
    }

    // Get barrier from upper MegaIteration, if available.
    const auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    if (megaIter != nullptr)
    {
        barrier = megaIter->getBarrier();
    }

    if (nullptr == barrier || nullptr == std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using Middle Point search, we need a DMultiMadsBarrier.");
    }

    if (evc->getCurrentEvalType() == NOMAD::EvalType::MODEL)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do Middle Point search for DMultiMads on EvalType::MODEL.");
    }

    // Keep a reference to the DMultiMads barrier.
    _ref_dmads_barrier = barrier;
}

void NOMAD::DMultiMadsMiddlePointSearchMethod::generateTrialPointsFinal()
{
    preRunValidations();

    const auto dMadsBarrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(_ref_dmads_barrier);

    // Get all the points available in the Pareto front approximation.
    auto getCurrentParetoFrontApproximation = [](
            const NOMAD::DMultiMadsBarrier &barrier) -> std::vector<NOMAD::EvalPointPtr>
    {
        if (barrier.getCurrentIncumbentFeas() != nullptr)
        {
            return barrier.getAllXFeas();
        }
        else
        {
            return barrier.getAllXInf();
        }
    };
    std::vector<NOMAD::EvalPointPtr> evalPointList = getCurrentParetoFrontApproximation(*dMadsBarrier);

    if (evalPointList.size() == 1)
    {
        OUTPUT_INFO_START
        AddOutputInfo("Only one point in the set of solutions. Stop Middle Point search step.");
        OUTPUT_INFO_END

        return;
    }

    const auto currentFrameCenter = dMadsBarrier->getCurrentIncumbentFeas() != nullptr ?
            dMadsBarrier->getCurrentIncumbentFeas() : dMadsBarrier->getCurrentIncumbentInf();

    if (evalPointList.size() == 2)
    {
        OUTPUT_INFO_START
        AddOutputInfo("Two points in the set of solutions. Take the middle point.");
        OUTPUT_INFO_END

        NOMAD::Point xcandidate(*evalPointList[0]->getX() + *evalPointList[1]->getX());
        xcandidate /= 2.0;

        // Generate candidate: no need to check if it exists in cache. We are ready to pay
        // the potential cache hit.
        auto evalPoint = NOMAD::EvalPoint(xcandidate);
        evalPoint.setPointFrom(std::make_shared<NOMAD::EvalPoint>(*currentFrameCenter),
                               NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        evalPoint.addGenStep(getStepType());
        insertTrialPoint(evalPoint);

        return;
    }

    // Obtain cache.
    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    const bool useCache = evc->getUseCache();
    const NOMAD::ArrayOfDouble lb = getPbParams()->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    const NOMAD::ArrayOfDouble ub = getPbParams()->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");

    const size_t nbMaxCacheSearchPerObj = getRunParams()->getAttributeValue<size_t>("DMULTIMADS_MIDDLEPOINT_SEARCH_CACHE_MAX");

    const size_t nbObj = dMadsBarrier->getNbObj();
    const NOMAD::FHComputeType initFHComputeType = dMadsBarrier->getFHComputeType();
    for (size_t obj = 0; obj < nbObj; ++obj)
    {
        std::sort(evalPointList.begin(), evalPointList.end(),
                  [&initFHComputeType, obj](const NOMAD::EvalPointPtr& ev1, const EvalPointPtr& ev2) -> bool
                  {
                      return ev1->getFs(initFHComputeType)[obj] < ev2->getFs(initFHComputeType)[obj];
                  });

        // Order points from largest to smallest gap according to the current objective.
        std::vector<std::tuple<NOMAD::Double, NOMAD::EvalPointPtr, NOMAD::EvalPointPtr>> gapValues(
                evalPointList.size() - 1);
        for (size_t i = 0; i < evalPointList.size() - 1; ++i)
        {
            const auto& evPrev = evalPointList[i];
            const auto& evNext = evalPointList[i + 1];
            const NOMAD::Double curGap = evNext->getFs(initFHComputeType)[obj] -
                                         evPrev->getFs(initFHComputeType)[obj];
            gapValues[i] = std::make_tuple(curGap, evPrev, evNext);
        }
        std::sort(gapValues.begin(), gapValues.end(),
                  [](const std::tuple<NOMAD::Double, NOMAD::EvalPointPtr, NOMAD::EvalPointPtr> &evt1,
                     const std::tuple<NOMAD::Double, NOMAD::EvalPointPtr, NOMAD::EvalPointPtr> &evt2) -> bool
                     {
                         return get<0>(evt1) > get<0>(evt2);
                     });

        if (get<0>(gapValues[0]) == 0)
        {
            OUTPUT_INFO_START
            AddOutputInfo("Maximal gap value for " + std::to_string(obj + 1) + " is zero: continue");
            OUTPUT_INFO_END

            continue;
        }

        // Starting from the largest gap, try to compute a middle point
        for (size_t i = 0; i < gapValues.size(); ++i)
        {
            if (i == nbMaxCacheSearchPerObj)
            {
                OUTPUT_INFO_START
                std::string s = "Middle Point search has reached the maximum number of cache search trials allowed for objective f";
                s += std::to_string(obj + 1);
                AddOutputInfo(s);
                OUTPUT_INFO_END
                break;
            }

            const auto& gapElt = gapValues[i];
            const auto gapValue = get<0>(gapElt);
            const auto evPrev = get<1>(gapElt);
            const auto evNext = get<2>(gapElt);

            NOMAD::Point xcandidate(*evPrev->getX() + *evNext->getX());
            xcandidate /= 2.0;

            OUTPUT_INFO_START
            AddOutputInfo("Middle Point search has generated the following point for objective f" + std::to_string(obj+1) + ":");
            std::string s = xcandidate.display();
            s += " from " + evPrev->getX()->display();
            s += " and " + evNext->getX()->display();
            s += " with gap in objective space " + gapValue.display();
            AddOutputInfo(s);
            OUTPUT_INFO_END

            // Generate a middle point candidate
            auto evalPoint = NOMAD::EvalPoint(xcandidate);
            evalPoint.setPointFrom(std::make_shared<NOMAD::EvalPoint>(*currentFrameCenter),
                                   NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
            evalPoint.addGenStep(getStepType());
        
            snapPointToBoundsAndProjectOnMesh(evalPoint, lb, ub);

            // When there is no cache, no need to continue the procedure.
            if (!useCache)
            {
                insertTrialPoint(evalPoint);
                break;
            }

            NOMAD::CacheInterface cacheInterface(this);
            NOMAD::EvalPoint foundEvalPoint;
            if (cacheInterface.find(*evalPoint.getX(), foundEvalPoint))
            {
                OUTPUT_INFO_START
                std::string s = "Middle point search candidate ";
                s += evalPoint.getX()->display();
                s += " already in cache: continue";
                AddOutputInfo(s);
                OUTPUT_INFO_END
                continue;
            }

            insertTrialPoint(evalPoint);
            break;
        }
    }
}
