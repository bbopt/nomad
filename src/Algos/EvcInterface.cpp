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

#include "../Algos/EvcInterface.hpp"
#include "../Algos/Iteration.hpp"
#include "../Algos/Mads/MegaSearchPoll.hpp"
#include "../Algos/PhaseOne/PhaseOne.hpp"
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

    _fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(_step);
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
// If not, add it to EvaluatorControl's points for sort.
std::vector<NOMAD::EvalPoint> NOMAD::EvcInterface::getSortedTrialPoints(const NOMAD::EvalPointSet &trialPoints, bool forceRandom, bool flagTrimIfNotDoEval)
{
    // Create EvalPoints and send them to EvaluatorControl
    if (nullptr == _evaluatorControl)
    {
        throw NOMAD::StepException(__FILE__,__LINE__, _step->getName() + ": EvaluatorControl not found", _step);
    }
    
    NOMAD::EvalType evalType = _evaluatorControl->getCurrentEvalType();
    if (NOMAD::EvalType::BB != evalType && NOMAD::EvalType::SURROGATE != evalType)
    {
        throw NOMAD::StepException(__FILE__,__LINE__, _step->getName() + ": Not suppose to sort points that are not BB", _step);
    }
    if (trialPoints.empty())
        return {};
    
    
    std::vector<EvalPoint> sortedTrialPoints;
    std::vector<EvalQueuePointPtr> evalPointsPtrToSort;
    
    OUTPUT_INFO_START
    _step->AddOutputInfo("Sort trial points for step " + _step->getName(), true, false);
    OUTPUT_INFO_END
    
    
    OUTPUT_INFO_START
    _step->AddOutputInfo(NOMAD::itos(trialPoints.size()) + " points added");
    OUTPUT_INFO_END
    
    for (auto trialPoint : trialPoints)
    {
        // Try to insert the trialPoint in the cache.
        // If it is not in the cache, it will be added.
        // If the point is already in the cache, depending on its
        // EvalStatus we might want to evaluate it.
        
        // First, convert trial point to full dimension, since we are
        // now only working with the cache and the EvaluatorControl.
        trialPoint = trialPoint.makeFullSpacePointFromFixed(_fixedVariable);

        bool doEval = true;
        if (flagTrimIfNotDoEval && _evaluatorControl->getUseCache())
        {
            const int maxNumberEval = 1;
            doEval = NOMAD::CacheBase::getInstance()->smartInsert(trialPoint, maxNumberEval, evalType);
            
        }
        if (doEval)
        {
            evalPointsPtrToSort.push_back(std::make_shared<EvalQueuePoint>(trialPoint,evalType));
            
            OUTPUT_DEBUG_START
            _step->AddOutputDebug("New point added for sorting: " + trialPoint.display());
            OUTPUT_DEBUG_END
            
        }
        else
        {
            // Point already evaluated
            OUTPUT_DEBUG_START
            _step->AddOutputDebug("Point already evaluated, trimmed from sorted eval points: " + trialPoint.display());
            OUTPUT_DEBUG_END
        }
    }
    
    OUTPUT_INFO_START
    _step->AddOutputDebug("Trial points sorted, before trim: " + std::to_string(trialPoints.size()) + " after trim: " + std::to_string(evalPointsPtrToSort.size()));
    OUTPUT_INFO_END
    
    if (!evalPointsPtrToSort.empty())
    {
        // Sort un-trimmed trial points
        _evaluatorControl->sort(evalPointsPtrToSort, forceRandom);
        
        for (const auto& evalPointPtr: evalPointsPtrToSort )
        {
            sortedTrialPoints.insert(sortedTrialPoints.begin(),evalPointPtr->makeSubSpacePointFromFixed(_fixedVariable));
        }
    }
    
    OUTPUT_INFO_START
    _step->AddOutputInfo("Sort trial points for step " + _step->getName(), false, true);
    NOMAD::OutputQueue::Flush();
    OUTPUT_INFO_END
        
    
    return sortedTrialPoints;
}


// For each point, look if it is in the cache.
// If it is, remove it from the EvalPointSet.
// If not, add it to EvaluatorControl's Queue.
void NOMAD::EvcInterface::keepPointsThatNeedEval(const NOMAD::EvalPointSet &trialPoints, bool useMesh)
{
    // Create EvalPoints and send them to EvaluatorControl
    if (nullptr == _evaluatorControl)
    {
        throw NOMAD::StepException(__FILE__,__LINE__, _step->getName() + ": EvaluatorControl not found", _step);
    }

    NOMAD::EvalType evalType = _evaluatorControl->getCurrentEvalType();

    // Currently, this method may be used inside an Iteration (Search or Poll, NM, ...),
    // or inside a MegaSearchPoll.
    auto iteration = _step->getParentOfType<NOMAD::Iteration*>();
    auto megaSearchPoll = dynamic_cast<const NOMAD::MegaSearchPoll*>(_step);

    if (useMesh && nullptr == iteration && nullptr == megaSearchPoll)
    {
        throw NOMAD::StepException(__FILE__,__LINE__, _step->getName() + ": In keepPointsThatNeedEval: need a parent of type Iteration or MegaSearchPoll", _step);
    }

    if (!trialPoints.empty())
    {
        OUTPUT_INFO_START
        _step->AddOutputInfo("Add points (full space) to eval queue for step " + _step->getName(), true, false);
        OUTPUT_INFO_END
        OUTPUT_INFO_START
        _step->AddOutputInfo(NOMAD::itos(trialPoints.size()) + " points added");
        OUTPUT_INFO_END
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
            // Cache is in full space.
            doEval = NOMAD::CacheBase::getInstance()->smartInsert(trialPoint, maxNumberEval, evalType);
        }
        else
        {

            // Look in EvaluatorControl's Barrier if the point is already evaluated.
            // Only do this when EvalType is BB. If it is MODEL, always reevaluate.
            if (NOMAD::EvalType::BB == evalType)
            {
                auto barrier = _evaluatorControl->getBarrier();
                if (nullptr != barrier)
                {
                    NOMAD::EvalPoint foundEvalPoint;
                    // Either point is not in barrier, or
                    // point was evaluated, but not with the current eval type.
                    // Barrier is in full space.
                    doEval = !barrier->findPoint(*trialPoint.getX(), foundEvalPoint)
                            || (nullptr == foundEvalPoint.getEval(evalType));
                }
            }
        }

        if (doEval)
        {
            NOMAD::EvalQueuePointPtr evalQueuePoint(new NOMAD::EvalQueuePoint(trialPoint, evalType));
            if (useMesh && nullptr == iteration)
            {
                std::string s = _step->getName();
                s += ": In keepPointsThatNeedEval: Could not determine iteration for point ";
                s += trialPoint.display();
                throw NOMAD::StepException(__FILE__,__LINE__, s, _step);
            }
            if ( useMesh )
            {
                auto mesh = iteration->getMesh();
                if ( mesh != nullptr )
                {
                    if (NOMAD::EvalType::BB == evalType || NOMAD::EvalType::SURROGATE == evalType )
                    {
                        evalQueuePoint->setMesh(mesh);
                    }
                    evalQueuePoint->setK(iteration->getK());
                }
            }

            evalQueuePoint->addGenStep(_step->getStepType());
            // Additional infos
            auto algo = _step->getParentOfType<NOMAD::Algorithm*>();
            while ( nullptr != algo)
            {
                evalQueuePoint->addGenStep(algo->getStepType());
                algo = algo->getParentOfType<NOMAD::Algorithm*>();
            }

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

    OUTPUT_INFO_START
    size_t evcNbPoints = _evaluatorControl->getQueueSize();
    if (evcNbPoints > 0)
    {
        // Not conditional to trialPoints size:
        // There could be leftover points to evaluate from previous evaluations.
        _step->AddOutputDebug("Eval queue now has " + NOMAD::itos(evcNbPoints) + " points.");
    }
    OUTPUT_INFO_END

    OUTPUT_INFO_START
    if (!trialPoints.empty())
    {
        _step->AddOutputInfo("Add points (full space) to eval queue for step " + _step->getName(), false, true);
    }

    NOMAD::OutputQueue::Flush();
    OUTPUT_INFO_END
}


// For each point, look if it is in the cache.
// If not, count it but DO NOT add it to EvaluatorControl's Queue.
size_t NOMAD::EvcInterface::countPointsThatNeedEval(const NOMAD::EvalPointSet &trialPoints)
{
    size_t nbPointsThatNeedEval = 0;
    // Create EvalPoints and send them to EvaluatorControl
    if (nullptr == _evaluatorControl)
    {
        throw NOMAD::StepException(__FILE__,__LINE__, _step->getName() + ": EvaluatorControl not found", _step);
    }

    NOMAD::EvalType evalType = _evaluatorControl->getCurrentEvalType();
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
            // Cache is in full space.
            doEval = NOMAD::CacheBase::getInstance()->smartInsert(trialPoint, maxNumberEval, evalType);
        }
        else
        {

            // Look in EvaluatorControl's Barrier if the point is already evaluated.
            // Only do this when EvalType is BB. If it is MODEL, always reevaluate.
            if (NOMAD::EvalType::BB == evalType)
            {
                auto barrier = _evaluatorControl->getBarrier();
                if (nullptr != barrier)
                {
                    NOMAD::EvalPoint foundEvalPoint;
                    // Either point is not in barrier, or
                    // point was evaluated, but not with the current eval type.
                    // Barrier is in full space.
                    doEval = !barrier->findPoint(*trialPoint.getX(), foundEvalPoint)
                            || (nullptr == foundEvalPoint.getEval(evalType));
                }
            }
        }

        if (doEval)
        {
            nbPointsThatNeedEval++;
        }
    }


    OUTPUT_INFO_START
    _step->AddOutputInfo("Number of points for step " + _step->getName() + " of eval type " + NOMAD::evalTypeToString(evalType) + " that would need eval: " + std::to_string(nbPointsThatNeedEval));

    NOMAD::OutputQueue::Flush();
    OUTPUT_INFO_END
    
    return nbPointsThatNeedEval;
}


std::vector<NOMAD::EvalPoint> NOMAD::EvcInterface::retrieveEvaluatedPointsFromCache(const NOMAD::EvalPointSet &trialPoints)
{
    std::vector<NOMAD::EvalPoint> evaluatedPoints;
    
    // Create EvalPoints and send them to EvaluatorControl
    if (nullptr == _evaluatorControl)
    {
        throw NOMAD::StepException(__FILE__,__LINE__, _step->getName() + ": EvaluatorControl not found", _step);
    }

    NOMAD::EvalType evalType = _evaluatorControl->getCurrentEvalType();

    if (_evaluatorControl->getUseCache())
    {
        
        for (auto trialPoint : trialPoints)
        {
            // Find the trialPoint in the cache.
            // If the point is in the cache, depending on its
            // EvalStatus, we don't need to evaluate it and we return it.
            
            // First, convert trial point to full dimension, since we are
            // now only working with the cache and the EvaluatorControl.
            
            trialPoint = trialPoint.makeFullSpacePointFromFixed(_fixedVariable);
            
            NOMAD::EvalPoint evalPoint;
            
            // Cache is in full space.
            NOMAD::CacheBase::getInstance()->find(trialPoint, evalPoint, evalType, false /* false: do not wait if point is not available */);
            
            if (evalPoint.isComplete() && evalPoint.isEvalOk(evalType) )
            {
                evalPoint = evalPoint.makeSubSpacePointFromFixed(_fixedVariable);
                evaluatedPoints.push_back(evalPoint);
            }
        }
    }

    return evaluatedPoints;
    
}


void NOMAD::EvcInterface::setBarrier(const std::shared_ptr<NOMAD::BarrierBase>& subBarrier)
{
    // The barrier may belong to a subspace but EvaluatorControl cares only of outputs for detecting success. No need create a barrier from scratch.
    _evaluatorControl->setBarrier(subBarrier);
}


std::vector<NOMAD::EvalPoint> NOMAD::EvcInterface::retrieveAllEvaluatedPoints()
{
    std::vector<NOMAD::EvalPoint> evaluatedPoints;

#ifdef _OPENMP
#pragma omp critical
#endif
    {
        for (auto evalPoint : _evaluatorControl->retrieveAllEvaluatedPoints())
        {
            // Convert from full to subspace dimension
            try
            {
                evalPoint = evalPoint.makeSubSpacePointFromFixed(_fixedVariable);
            }
            catch(...)
            {
                OUTPUT_INFO_START
                _step->AddOutputInfo("Fail to convert from full to subspace. Point is not retrieved. Let's continue in " + _step->getName(), true, false);
                OUTPUT_INFO_END
                continue;
            }
            evaluatedPoints.push_back(evalPoint);
        }
    }
    return evaluatedPoints;
}


// When points are generated and added to queue,
// we can start evaluation.
NOMAD::SuccessType NOMAD::EvcInterface::startEvaluation()
{
    OUTPUT_INFO_START
    _step->AddOutputInfo("Evaluate points (full space) for " + _step->getName(), true, false);
    OUTPUT_INFO_END

    std::shared_ptr<NOMAD::AllStopReasons> stopReasons = _step->getAllStopReasons();

    // Evaluate points
    // Note: do not use checkTerminate() here. If it is time to terminate, EvaluatorControl will take
    // care of clearing the queue.
    NOMAD::SuccessType success = _evaluatorControl->run();

    OUTPUT_DEBUG_START
    std::string s = _step->getName() + ": " + NOMAD::enumStr(success);
    s += ". Stop reasons: " + stopReasons->getStopReasonAsString() ;
    _step->AddOutputDebug(s);
    OUTPUT_DEBUG_END

    OUTPUT_INFO_START
    _step->AddOutputInfo("Evaluate points (full space) for " + _step->getName(), false, true);
    NOMAD::OutputQueue::Flush();
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
