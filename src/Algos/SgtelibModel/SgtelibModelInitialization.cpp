
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelInitialization.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::SgtelibModelInitialization::init()
{
    _name = NOMAD::Initialization::getName();
    verifyParentNotNull();
}


/*-------------------------*/
/*       Destructor        */
/*-------------------------*/
NOMAD::SgtelibModelInitialization::~SgtelibModelInitialization()
{
}


void NOMAD::SgtelibModelInitialization::startImp()
{
}


bool NOMAD::SgtelibModelInitialization::runImp()
{
    bool doContinue = ! _stopReasons->checkTerminate();

    if (doContinue)
    {
        eval_x0s();
        doContinue = ! _stopReasons->checkTerminate();
    }

    return doContinue;
}


void NOMAD::SgtelibModelInitialization::endImp()
{
}


void NOMAD::SgtelibModelInitialization::validateX0s() const
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


// Eval X0s, using blackbox.
// If we are here, it means we are in standalone mode. Either X0s were provided,
// or the best points were found in the cache, in MainStep, and put in parameter
// X0.
// Method is copied from MadsInitialization.
bool NOMAD::SgtelibModelInitialization::eval_x0s()
{
    bool evalOk = false;

    auto x0s = _pbParams->getAttributeValue<NOMAD::ArrayOfPoint>("X0");

    validateX0s();

    // Add X0s that need evaluation to eval queue
    NOMAD::CacheInterface cacheInterface(this);
    NOMAD::EvcInterface evcInterface(this);
    auto evc = evcInterface.getEvaluatorControl();
    evc->lockQueue();

    NOMAD::EvalPointSet evalPointSet;
    for (size_t x0index = 0; x0index < x0s.size(); x0index++)
    {
        auto x0 = x0s[x0index];
        NOMAD::EvalPoint evalPoint_x0(x0);
        evalPointSet.insert(evalPoint_x0);
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

    auto evalPointList = evcInterface.retrieveAllEvaluatedPoints();
    for (auto x0 : x0s)
    {
        NOMAD::EvalPoint evalPoint_x0(x0);
        cacheInterface.find(x0, evalPoint_x0);
        // To evaluate X0, use blackbox, not sgte.
        if (evalPoint_x0.isEvalOk(NOMAD::EvalType::BB))
        {
            // evalOk is true if at least one evaluation is Ok
            evalOk = true;
            AddOutputInfo("Using X0: " + evalPoint_x0.displayAll());
        }
        else
        {
            AddOutputError("X0 evaluation failed for X0 = " + x0.display());
        }
    }

    if (evalOk)
    {
        // Construct barrier using x0s
        auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
        _barrier = std::make_shared<NOMAD::Barrier>(hMax, NOMAD::SubproblemManager::getSubFixedVariable(this), evc->getEvalType(), evalPointList);
    }
    else
    {
        auto sgteStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        sgteStopReason->set(NOMAD::ModelStopType::X0_FAIL);
    }

    NOMAD::OutputQueue::Flush();

    return evalOk;
}
