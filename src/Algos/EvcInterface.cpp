
#include "../Algos/EvcInterface.hpp"
#include "../Algos/Iteration.hpp"
#include "../Algos/Mads/MegaSearchPoll.hpp"
#include "../Algos/SubproblemManager.hpp"
#include "../Cache/CacheBase.hpp"
#include "../Output/OutputQueue.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
std::shared_ptr<NOMAD::EvaluatorControl> NOMAD::EvcInterface::_evaluatorControl = nullptr;

void NOMAD::EvcInterface::init()
{
    verifyStepNotNull();
    verifyEvaluatorControlNotNull();

    _fixedVariable = NOMAD::SubproblemManager::getSubFixedVariable(_step);
}


void NOMAD::EvcInterface::verifyStepNotNull()
{
    if (nullptr == _step)
    {
        std::string err = "Step for EvcInterface should not be NULL";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::EvcInterface::verifyEvaluatorControlNotNull()
{
    if (nullptr == _evaluatorControl)
    {
        std::string err = "EvaluatorControl for EvcInterface should not be NULL";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::EvcInterface::setEvaluatorControl(const std::shared_ptr<NOMAD::EvaluatorControl>& evaluatorControl)
{
    _evaluatorControl = evaluatorControl;
    verifyEvaluatorControlNotNull();
}


// For each point, look if it is in the cache.
// If it is, remove it from the EvalPointSet.
// If not, add it to EvaluatorControl's Queue.
void NOMAD::EvcInterface::keepPointsThatNeedEval(const NOMAD::EvalPointSet &trialPoints, bool useMesh)
{
    NOMAD::EvalType evalType = _evaluatorControl->getEvalType();

    // Create EvalPoints and send them to EvaluatorControl
    if (nullptr == _evaluatorControl)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, _step->getName() + ": EvaluatorControl not found");
    }

    // Currently, this method may be used inside an Iteration (Search or Poll, NM, ...),
    // or inside a MegaSearchPoll.
    auto iteration = _step->getParentOfType<NOMAD::Iteration*>();
    auto megaSearchPoll = dynamic_cast<const NOMAD::MegaSearchPoll*>(_step);

    if (useMesh && nullptr == iteration && nullptr == megaSearchPoll)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, _step->getName() + ": In keepPointsThatNeedEval: need a parent of type Iteration or MegaSearchPoll");
    }

    if (trialPoints.size() > 0)
    {
        OUTPUT_INFO_START
        _step->AddOutputInfo("Add points to eval queue for step " + _step->getName(), true, false);
        OUTPUT_INFO_END
        OUTPUT_DEBUG_START
        _step->AddOutputDebug(NOMAD::itos(trialPoints.size()) + " points to add to eval queue");
        OUTPUT_DEBUG_END
    }

    for (auto trialPoint : trialPoints)
    {
        // Try to insert the trialPoint in the cache.
        // If it is not in the cache, it will be added.
        // If the point is already in the cache, depending on its
        // EvalStatus we might want to evaluate it.

        // First, convert trial point to full dimension, since we are
        // now only working with the cache and the EvaluatorControl.
        auto trialPointSub = trialPoint;    // Used to get iteration
        trialPoint = trialPoint.makeFullSpacePointFromFixed(_fixedVariable);

        // Compute if we should evaluate, maybe re-evaluate, this point
        bool doEval = true;
        // maxNumberEval is used to compute if we should re-evaluate this point.
        // Default value is 1.
        // This will be a parameter in the future. Currently not implemented.
        const int maxNumberEval = 1;
        if (_evaluatorControl->getUseCache())
        {
            doEval = NOMAD::CacheBase::getInstance()->smartInsert(trialPoint, maxNumberEval, evalType);
        }
        else
        {
            // smartInsert would have taken care of updating tag, but since
            // cache is not used, update tag here.
            trialPoint.updateTag();
            // Look in EvaluatorControl's Barrier if the point is already evaluated.
            auto barrier = _evaluatorControl->getBarrier();
            if (nullptr != barrier)
            {
                NOMAD::EvalPoint foundEvalPoint;
                // Either point is not in barrier, or
                // point was evaluated, but not with the current eval type.
                doEval = !findInList(*trialPoint.getX(), barrier->getAllPoints(), foundEvalPoint)
                         || (nullptr == foundEvalPoint.getEval(evalType));
                {
                    doEval = true;
                }
            }
        }

        if (doEval)
        {
            NOMAD::EvalQueuePointPtr evalQueuePoint(new NOMAD::EvalQueuePoint(trialPoint, evalType));
            if (useMesh && nullptr == iteration)
            {
                iteration = dynamic_cast<NOMAD::Iteration*>(megaSearchPoll->getIterForPoint(trialPointSub).get());
                if (nullptr == iteration)
                {
                    std::string s = _step->getName();
                    s += ": In keepPointsThatNeedEval: Could not determine iteration for point ";
                    s += trialPoint.display();
                    throw NOMAD::Exception(__FILE__,__LINE__, s);
                }
            }
            if ( useMesh )
            {
                auto mesh = iteration->getMesh();
                if ( mesh != nullptr )
                {
                    evalQueuePoint->setMeshSize(mesh->getdeltaMeshSize());
                    evalQueuePoint->setFrameSize(mesh->getDeltaFrameSize());
                    evalQueuePoint->setK(iteration->getK());
                }
            }

            evalQueuePoint->setComment(_step->getAlgoComment());
            evalQueuePoint->setGenStep(_step->getName());
#ifdef _OPENMP
            evalQueuePoint->setThreadAlgo(omp_get_thread_num());
#endif // _OPENMP

            if (_evaluatorControl->addToQueue(evalQueuePoint))
            {
                OUTPUT_DEBUG_START
                _step->AddOutputDebug("New point added to eval queue: " + trialPoint.display());
                OUTPUT_DEBUG_END
            }
            else
            {
                OUTPUT_DEBUG_START
                _step->AddOutputDebug("Point not added to eval queue: " + trialPoint.display());
                OUTPUT_DEBUG_END
            }
        }
        else
        {
            // Point already evaluated
            OUTPUT_DEBUG_START
            _step->AddOutputDebug("Point not re-evaluated: " + trialPoint.display());
            OUTPUT_DEBUG_END
        }
    }

    OUTPUT_DEBUG_START
    size_t evcNbPoints = _evaluatorControl->getQueueSize();
    if (evcNbPoints > 0)
    {
        // Not conditional to trialPoints size:
        // There could be leftover points to evaluate from previous evaluations.
        _step->AddOutputDebug("Eval queue now has " + NOMAD::itos(evcNbPoints) + " points.");
    }
    OUTPUT_DEBUG_END

    OUTPUT_INFO_START
    if (trialPoints.size() > 0)
    {
        _step->AddOutputInfo("Add points to eval queue for step " + _step->getName(), false, true);
    }

    NOMAD::OutputQueue::Flush();
    OUTPUT_INFO_END
}


void NOMAD::EvcInterface::setBarrier(const std::shared_ptr<NOMAD::Barrier>& subBarrier)
{
    if (subBarrier == nullptr)
    {
        return;
    }

    // Input is the barrier from MegaIteration, which may belong to a subspace.
    // EvaluatorControl's barrier must be in full dimension.
    auto fullBarrier = std::make_shared<NOMAD::Barrier>(*subBarrier);

    // Clear xFeas and xInf lists and recompute them
    fullBarrier->clearXFeas();
    fullBarrier->clearXInf();
    for (auto xFeas : subBarrier->getAllXFeas())
    {
        auto xFeasFull = xFeas.makeFullSpacePointFromFixed(_fixedVariable);
        fullBarrier->addXFeas(xFeasFull, _evaluatorControl->getEvalType());
    }
    for (auto xInf : subBarrier->getAllXInf())
    {
        auto xInfFull = xInf.makeFullSpacePointFromFixed(_fixedVariable);
        fullBarrier->addXInf(xInfFull);
    }
    auto refBestFeas = subBarrier->getRefBestFeas();
    auto refBestInf  = subBarrier->getRefBestInf();
    if (nullptr != refBestFeas)
    {
        fullBarrier->setRefBestFeas(std::make_shared<NOMAD::EvalPoint>(refBestFeas->makeFullSpacePointFromFixed(_fixedVariable)));
    }
    if (nullptr != refBestInf)
    {
        fullBarrier->setRefBestInf(std::make_shared<NOMAD::EvalPoint>(refBestInf->makeFullSpacePointFromFixed(_fixedVariable)));
    }

    _evaluatorControl->setBarrier(fullBarrier);
}


std::vector<NOMAD::EvalPoint> NOMAD::EvcInterface::retrieveAllEvaluatedPoints()
{
    std::vector<NOMAD::EvalPoint> evaluatedPoints;

    for (auto evalPoint : _evaluatorControl->retrieveAllEvaluatedPoints())
    {
        // Convert from full to subspace dimension
        evalPoint = evalPoint.makeSubSpacePointFromFixed(_fixedVariable);
        evaluatedPoints.push_back(evalPoint);
    }

    return evaluatedPoints;
}


// When points are generated and added to queue,
// we can start evaluation.
NOMAD::SuccessType NOMAD::EvcInterface::startEvaluation()
{
    OUTPUT_INFO_START
    _step->AddOutputInfo("Evaluate points for " + _step->getName(), true, false);
    OUTPUT_INFO_END

    NOMAD::SuccessType success = NOMAD::SuccessType::UNSUCCESSFUL;
    std::shared_ptr<NOMAD::AllStopReasons> stopReasons = _step->getAllStopReasons();

    if ( ! stopReasons->checkTerminate() )
    {
        // Evaluate points
        success = _evaluatorControl->run();
    }

    OUTPUT_DEBUG_START
    std::string s = _step->getName() + ": " + NOMAD::enumStr(success);
    s += ". Stop reasons: " + stopReasons->getStopReasonAsString() ;
    _step->AddOutputDebug(s);
    OUTPUT_DEBUG_END

    OUTPUT_INFO_START
    NOMAD::OutputQueue::Flush();

    _step->AddOutputInfo("Evaluate points for " + _step->getName(), false, true);
    OUTPUT_INFO_END

    return success;
}


bool NOMAD::EvcInterface::evalSinglePoint(NOMAD::EvalPoint &evalPoint,
                                          const NOMAD::Double &hMax)
{
    // Convert to full dimension before calling EvaluatorControl
    evalPoint = evalPoint.makeFullSpacePointFromFixed(_fixedVariable);
    bool ret = _evaluatorControl->evalSinglePoint(evalPoint, NOMAD::getThreadNum(), hMax);
    // Convert back to subspace dimension
    evalPoint = evalPoint.makeSubSpacePointFromFixed(_fixedVariable);

    return ret;
}
