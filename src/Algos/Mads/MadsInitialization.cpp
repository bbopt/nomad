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

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/MicroSleep.hpp"

void NOMAD::MadsInitialization::init()
{
    setStepType(NOMAD::StepType::INITIALIZATION);
}


bool NOMAD::MadsInitialization::runImp()
{
    _initialMesh = std::make_shared<NOMAD::GMesh>(_pbParams);

    bool doContinue = ! _stopReasons->checkTerminate();

    if (doContinue)
    {
        eval_x0s();
        doContinue = ! _stopReasons->checkTerminate();
    }
    return doContinue;
}


void NOMAD::MadsInitialization::validateX0s() const
{
    auto x0s = _pbParams->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    bool validX0available = false;
    std::string err;

    for (size_t x0index = 0; x0index < x0s.size(); x0index++)
    {
        auto x0 = x0s[x0index];
        if (!x0.isComplete() || x0.size() != n)
        {
            err += "Initialization: eval_x0s: Invalid X0 " + x0.display() + ".";
        }
        else
        {
            validX0available = true;
        }
    }

    if (validX0available)
    {
        if (!err.empty())
        {
            // Show invalid X0s
            AddOutputWarning(err);
        }
    }
    else
    {
        // No valid X0 available. Throw exception.
        size_t cacheSize = NOMAD::CacheBase::getInstance()->size();
        if (cacheSize > 0)
        {
            err += " Hint: Try not setting X0 so that the cache is used (";
            err += std::to_string(cacheSize) + " points).";
        }
        else
        {
            err += " Cache is empty. Hint: Try setting LH_SEARCH so that the Latin Hypercube Search is used to find initial points.";
        }
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

}


bool NOMAD::MadsInitialization::eval_x0s()
{
    bool evalOk = false;
    std::string s;

    auto x0s = _pbParams->getAttributeValue<NOMAD::ArrayOfPoint>("X0");

    validateX0s();

    // Add X0s that need evaluation to eval queue
    NOMAD::CacheInterface cacheInterface(this);
    NOMAD::EvcInterface evcInterface(this);
    auto evc = evcInterface.getEvaluatorControl();
    NOMAD::EvalType evalType = NOMAD::EvalType::BB;
    NOMAD::ComputeType computeType = NOMAD::ComputeType::STANDARD;
    if (nullptr != evc)
    {
        evalType = evc->getEvalType();
        computeType = evc->getComputeType();
        evc->lockQueue();
    }

    NOMAD::EvalPointSet evalPointSet;
    for (size_t x0index = 0; x0index < x0s.size(); x0index++)
    {
        auto x0 = x0s[x0index];
        NOMAD::EvalPoint evalPointX0(x0);
        evalPointSet.insert(evalPointX0);
    }

    // Add points to the eval queue.
    // Convert to full dimension if needed.
    // Note: Queue is already locked - it needs to be locked to add points.
    evcInterface.keepPointsThatNeedEval(evalPointSet, false);   // false: no mesh

    if (nullptr != evc)
    {
        // Enforce no opportunism.
        auto previousOpportunism = evc->getOpportunisticEval();
        evc->setOpportunisticEval(false);
        evc->unlockQueue(false); // false: do not sort eval queue

        // Evaluate all x0s. Ignore returned success type.
        // Note: EvaluatorControl would not be able to compare/compute success since there is no barrier.
        evcInterface.startEvaluation();

        // Reset opportunism to previous values.
        evc->setOpportunisticEval(previousOpportunism);
    }

    bool x0Failed = true;

    // Construct barrier using points evaluated by this step.
    // The points are cleared from the EvaluatorControl.
    auto evaluatedPoints = evcInterface.retrieveAllEvaluatedPoints();
    std::vector<NOMAD::EvalPoint> evalPointX0s;

    if (_stopReasons->checkTerminate())
    {
        return false;
    }

    for (auto x0 : x0s)
    {
        NOMAD::EvalPoint evalPointX0(x0);

        // Look for x0 in freshly evaluated points
        bool x0Found = findInList(x0, evaluatedPoints, evalPointX0);

        if (!x0Found)
        {
            // Look for x0 in EvaluatorControl barrier
            x0Found = evcInterface.findInBarrier(x0, evalPointX0);
            if (!x0Found)
            {
                // Look for x0 in cache
                // Note: Even if we are not currently using cache in this sub-algorithm,
                // we may have interesting points in the global cache.
                if (nullptr != evc && evc->getUseCache())
                {
                    // If status of point in cache is IN_PROGRESS, wait for evaluation to be completed.
                    x0Found = (cacheInterface.find(x0, evalPointX0, evalType) > 0);
                }
                else
                {
                    // Look for X0 in cache, but do not wait for evaluation.
                    x0Found = (cacheInterface.find(x0, evalPointX0) > 0);
                }
            }
        }

        if (x0Found && evalPointX0.isEvalOk(evalType))
        {
            // evalOk is true if at least one evaluation is Ok
            evalOk = true;
            evalPointX0s.push_back(evalPointX0);

            x0Failed = false;   // At least one good X0.
        }
    }

    if (x0Failed)
    {
        // All x0s failed. Show an error.
        auto madsStopReason = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get(_stopReasons);
        madsStopReason->set(NOMAD::MadsStopType::X0_FAIL);

        for (auto x0 : x0s)
        {
            auto x0Full = x0.makeFullSpacePointFromFixed(NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
            AddOutputError("X0 evaluation failed for X0 = " + x0Full.display());
        }
    }
    else
    {
        OUTPUT_INFO_START
        for (auto evalPointX0 : evalPointX0s)
        {
            s = "Using X0: ";
            // BB: Simple display. MODEL: Full display.
            s += (NOMAD::EvalType::BB == evalType) ? evalPointX0.display() : evalPointX0.displayAll();
        }
        AddOutputInfo(s);
        OUTPUT_INFO_END

        // Construct barrier using x0s
        auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
        _barrier = std::make_shared<NOMAD::Barrier>(hMax,
                                                    NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                                    evalType,
                                                    computeType,
                                                    evalPointX0s,
                                                    _barrierInitializedFromCache);
    }

    NOMAD::OutputQueue::Flush();

    return evalOk;
}
