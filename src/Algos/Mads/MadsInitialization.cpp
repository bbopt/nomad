
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::MadsInitialization::init()
{
    _name = NOMAD::Initialization::getName();
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
            err += std::to_string(cacheSize) + " points)";
        }
        else
        {
            err += ". Cache is empty.";
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
    auto evalType = evc->getEvalType();
    evc->lockQueue();

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

    // Enforce no opportunism.
    auto previousOpportunism = evc->getOpportunisticEval();
    evc->setOpportunisticEval(false);
    evc->unlockQueue(false); // false: do not sort eval queue

    // Evaluate all x0s. Ignore returned success type.
    // Note: EvaluatorControl would not be able to compare/compute success since there is no barrier.
    evcInterface.startEvaluation();

    // Reset opportunism to previous values.
    evc->setOpportunisticEval(previousOpportunism);

    bool x0Failed = true;

    // Construct barrier using points evaluated by this step.
    // The points are cleared from the EvaluatorControl.
    auto evaluatedPoints = evcInterface.retrieveAllEvaluatedPoints();
    std::vector<NOMAD::EvalPoint> evalPointX0s;
    for (auto x0 : x0s)
    {
        NOMAD::EvalPoint evalPointX0(x0);

        // Look for x0 in freshly evaluated points
        bool x0Found = findInList(x0, evaluatedPoints, evalPointX0);

        if (!x0Found)
        {
            auto barrier = evc->getBarrier();
            if (nullptr != barrier)
            {
                // Look for x0 in EvaluatorControl barrier
                x0Found = findInList(x0, barrier->getAllPoints(), evalPointX0);
            }
            if (!x0Found && evc->getUseCache())
            {
                // Look for x0 in cache
                x0Found = (cacheInterface.find(x0, evalPointX0) > 0);
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
            AddOutputError("X0 evaluation failed for X0 = " + x0.display());
        }
    }
    else
    {
        OUTPUT_INFO_START
        for (auto evalPointX0 : evalPointX0s)
        {
            s = "Using X0: ";
            // BB: Simple display. SGTE: Full display.
            s += (NOMAD::EvalType::BB == evalType) ? evalPointX0.display() : evalPointX0.displayAll();
        }
        AddOutputInfo(s);
        OUTPUT_INFO_END

        // Construct barrier using x0s
        auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
        _barrier = std::make_shared<NOMAD::Barrier>(hMax, NOMAD::SubproblemManager::getSubFixedVariable(this), evalType, evalPointX0s);
    }

    NOMAD::OutputQueue::Flush();

    return evalOk;
}
