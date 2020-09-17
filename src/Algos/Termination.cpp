
#include "../Algos/Algorithm.hpp"
#include "../Algos/EvcInterface.hpp"
#include "../Algos/Termination.hpp"
#include "../Cache/CacheBase.hpp"
#include "../Util/Clock.hpp"

void NOMAD::Termination::init()
{
    // Usually, we do not have the Algorithm name yet, so we cannot use it here
    // for _name
    _name = "Termination";
    verifyParentNotNull();
}


void NOMAD::Termination::startImp()
{
    _name = getAlgoName() + "Termination";
}


bool NOMAD::Termination::runImp()
{
    return _stopReasons->checkTerminate() ;
}


bool NOMAD::Termination::terminate(size_t iteration)
{
    bool stop = _stopReasons->checkTerminate();
    if (stop)
    {
        // A stop condition was already reached.
        return stop;
    }

    // Set stopReason due to criterions other than AlgoStopReasons<>
    auto maxIterations = _runParams->getAttributeValue<size_t>("MAX_ITERATIONS");
    auto maxTime = _runParams->getAttributeValue<size_t>("MAX_TIME");

    // Termination conditions go here.
    // This is also tested in EvaluatorControl
    if (NOMAD::Step::getUserTerminate())
    {
        // Force quit (by pressing CTRL-C):
        _stopReasons->set(NOMAD::BaseStopType::CTRL_C);
    }
    else if (maxIterations < NOMAD::INF_SIZE_T && iteration >= maxIterations)
    {
        // Max iterations reached
        _stopReasons->set(NOMAD::IterStopType::MAX_ITER_REACHED);
    }
    else if (maxTime < NOMAD::INF_SIZE_T && NOMAD::Clock::getTimeSinceStart() >= maxTime)
    {
        // Max time reached
        _stopReasons->set(NOMAD::BaseStopType::MAX_TIME_REACHED);
    }
    else if (_pbParams->getAttributeValue<bool>("STOP_IF_FEASIBLE") && solHasFeas())
    {
        _stopReasons->set(NOMAD::IterStopType::STOP_ON_FEAS );
    }
    else
    {
        // Need to check on MaxEval and MaxBBEval a last time because in evaluatorControl
        // the stop reason may have been set due to all queue points evaluated.
        auto evc = NOMAD::EvcInterface::getEvaluatorControl();
        stop = evc->reachedMaxEval() || evc->reachedMaxStepEval();
    }

    stop = stop || _stopReasons->checkTerminate() ;
    return stop;
}


bool NOMAD::Termination::solHasFeas() const
{
    bool hasFeas = NOMAD::CacheBase::getInstance()->hasFeas();
    if (!hasFeas)
    {
        // No feasible point in cache, but possibly in parent step's barrier.
        if (nullptr != _parentStep)
        {
            auto barrier = _parentStep->getMegaIterationBarrier();
            hasFeas = (nullptr != barrier && nullptr != barrier->getFirstXFeas());
        }
    }

    return hasFeas;
}


void NOMAD::Termination::endImp()
{
    const NOMAD::Algorithm* currentAlgo = getParentOfType<NOMAD::Algorithm*>();
    NOMAD::OutputLevel outputLevel = currentAlgo->isSubAlgo() ? NOMAD::OutputLevel::LEVEL_INFO
                                                              : NOMAD::OutputLevel::LEVEL_HIGH;

    if (_stopReasons->checkTerminate())
    {
        std::string terminationInfo = "A termination criterion is reached: ";
        terminationInfo += _stopReasons->getStopReasonAsString();
        auto evc = NOMAD::EvcInterface::getEvaluatorControl();
        if (_stopReasons->testIf(NOMAD::EvalStopType::MAX_BB_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getBbEval());
        }
        else if (_stopReasons->testIf(NOMAD::EvalStopType::MAX_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getNbEval());
        }
        else if (_stopReasons->testIf(NOMAD::EvalStopType::MAX_BLOCK_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getBlockEval());
        }
        else if (_stopReasons->testIf(NOMAD::EvalStopType::MAX_SGTE_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getTotalSgteEval());
        }
        else if (_stopReasons->testIf(NOMAD::EvalStopType::LAP_MAX_BB_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getLapBbEval());
        }
        AddOutputInfo(terminationInfo, outputLevel);
    }
    else
    {
        std::string terminationInfo = "No termination criterion reached";
        AddOutputInfo(terminationInfo, outputLevel);
    }

}
