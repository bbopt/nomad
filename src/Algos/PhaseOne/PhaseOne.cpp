
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/PhaseOne/PhaseOne.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Eval/ComputeSuccessType.hpp"

NOMAD::BBOutputTypeList NOMAD::PhaseOne::_bboutputtypes;

void NOMAD::PhaseOne::init()
{
    _name = "Phase One";
    verifyParentNotNull();

}

void NOMAD::PhaseOne::startImp()
{

    // Setup EvalPoint success computation to be based on h rather than f.
    NOMAD::EvcInterface::getEvaluatorControl()->setComputeSuccessTypeFunction(NOMAD::ComputeSuccessType::computeSuccessTypePhaseOne);
    NOMAD::Eval::setComputeSuccessTypeFunction(NOMAD::Eval::computeSuccessTypePhaseOne);
    NOMAD::Eval::setComputeHFunction(NOMAD::Eval::computeHPB);

    // The cache may not be empty.
    // Recompute the h for cache points that were read from cache file.
    NOMAD::CacheBase::getInstance()->processOnAllPoints(NOMAD::PhaseOne::recomputeHPB);

    // Comment to appear at the end of stats lines
    setAlgoComment("(Phase One)", true); // true: force comment

    // Setup the pb parameters to stop once a feasible point is obtained
    _pbParams->setAttributeValue("STOP_IF_FEASIBLE", true);
    _pbParams->checkAndComply();


    // Setup Mads
    _madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
    _mads = std::make_shared<NOMAD::Mads>(this, _madsStopReasons, _runParams, _pbParams);

}

void NOMAD::PhaseOne::readInformationForHotRestart()
{
}



bool NOMAD::PhaseOne::runImp()
{
    bool ret = false;

    // Run Mads on Phase One.
    _mads->start();
    ret = _mads->run();
    _mads->end();

    return ret;
}


void NOMAD::PhaseOne::endImp()
{
    // Remove any remaining points from eval queue.
    EvcInterface::getEvaluatorControl()->clearQueue();
    // Ensure evaluation of queue will continue
    NOMAD::EvcInterface::getEvaluatorControl()->restart();

    // reset to the previous stats comment
    resetPreviousAlgoComment(true); // true: release lock on comment

    // Reset success computation function
    NOMAD::EvcInterface::getEvaluatorControl()->setComputeSuccessTypeFunction(NOMAD::ComputeSuccessType::defaultComputeSuccessType);
    // VRM THIS ALSO HAS TO BE FIXED
    NOMAD::Eval::setComputeSuccessTypeFunction(NOMAD::Eval::defaultComputeSuccessType);
    NOMAD::Eval::setComputeHFunction(NOMAD::Eval::defaultComputeH);

    // All points in the cache must be recomputed for their h.
    // Note: Cache is ordered on the Point part only, and we recompute the Eval
    // part, so the cache remains coherent.
    NOMAD::CacheBase::getInstance()->processOnAllPoints(NOMAD::PhaseOne::recomputeH);

    bool hasFeas = NOMAD::CacheBase::getInstance()->hasFeas();
    if (!hasFeas)
    {
        // If cache is not used, feasible points remain in the barrier
        auto barrier = _mads->getMegaIterationBarrier();
        if (nullptr != barrier)
        {
            hasFeas = (nullptr != barrier->getFirstXFeas());
        }
    }

    // Update PhaseOne stop reasons
    auto PhaseOneStopReasons = NOMAD::AlgoStopReasons<NOMAD::PhaseOneStopType>::get( _stopReasons );
    if (!hasFeas)
    {
        if ( _madsStopReasons->checkTerminate() )
        {
            PhaseOneStopReasons->set ( NOMAD::PhaseOneStopType::MADS_FAIL );
        }
        else
        {
            PhaseOneStopReasons->set ( NOMAD::PhaseOneStopType::NO_FEAS_PT );
        }
    }
}


void NOMAD::PhaseOne::recomputeH(NOMAD::EvalPoint& evalPoint)
{
    // EvalType BB: Never use Sgte in Phase One
    auto eval = evalPoint.getEval(NOMAD::EvalType::BB);
    if (nullptr != eval && !eval->getBBO().empty())
    {
        eval->setH(NOMAD::Eval::defaultComputeH(*eval, _bboutputtypes));
    }
}


void NOMAD::PhaseOne::recomputeHPB(NOMAD::EvalPoint& evalPoint)
{
    auto eval = evalPoint.getEval(NOMAD::EvalType::BB);
    if (nullptr != eval && !eval->getBBO().empty())
    {
        eval->setH(NOMAD::Eval::computeHPB(*eval, _bboutputtypes));
    }
}


