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
#include "../Algos/IterationUtils.hpp"
#include "../Algos/Mads/Search.hpp"
#include "../Algos/SubproblemManager.hpp"
#include "../Algos/SurrogateEvaluation.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Type/EvalSortType.hpp"
#include "../Algos/MainStep.hpp"

void NOMAD::IterationUtils::init()
{
    // Set values of _fromAlgo and _iterAncestor.
    
    // Is the direct parent an Algorithm?
    auto algoParent = dynamic_cast<const NOMAD::Algorithm*>(_parentStep);
    _fromAlgo = (nullptr != algoParent);
    
    
    // Is the direct parent an Iteration?
    auto iterParent = dynamic_cast<const NOMAD::Iteration*>(_parentStep);
    
    
    // Check if there is an Iteration among ancestors
    _iterAncestor = const_cast<NOMAD::Iteration*>(iterParent);
    if (nullptr == _iterAncestor)
    {
        auto iterAncestorConst = _parentStep->getParentOfType<NOMAD::Iteration*>();
        _iterAncestor = const_cast<NOMAD::Iteration*>(iterAncestorConst);
    }
    
    // Find the MegaIteration ancestor.
    // Is this IterationUtils a MegaIteration?
    auto megaIter = dynamic_cast<const NOMAD::MegaIteration*>(this);
    if (nullptr == megaIter)
    {
        // Is the direct parent a MegaIteration?
        megaIter = dynamic_cast<const NOMAD::MegaIteration*>(_parentStep);
        if (nullptr == megaIter)
        {
            // Check if there is a MegaIteration among ancestors
            megaIter = _parentStep->getParentOfType<NOMAD::MegaIteration*>();
        }
    }
    _megaIterAncestor = const_cast<NOMAD::MegaIteration*>(megaIter);
    
    if (!_fromAlgo && nullptr == _iterAncestor && nullptr == _megaIterAncestor)
    {
        throw NOMAD::StepException(__FILE__, __LINE__,
                                   "An instance of class IterationUtils must have either an Iteration or a MegaIteration as ancestor or an Algorithm as direct parent",
                                   _parentStep);
    }
    
    const NOMAD::Search * search = dynamic_cast<const NOMAD::Search*>(_parentStep);
    if ( nullptr == search)
    {
        if ( nullptr != _iterAncestor)
        {
            search = _iterAncestor->getParentOfType<NOMAD::Search*>(false /*do not stop at algo*/);
        }
        else if ( nullptr != _megaIterAncestor)
        {
            search = _megaIterAncestor->getParentOfType<NOMAD::Search*>(false /*do not stop at algo*/);
        }
    }
    if (nullptr != search && nullptr != search->getRunParams())
    {
        _projectOnMesh = search->getRunParams()->getAttributeValue<bool>("SEARCH_METHOD_MESH_PROJECTION");
    }
    
    auto runParams = _parentStep->getRunParams();
    _frameCenterUseCache = false;
    if (nullptr != runParams )
    {
        _frameCenterUseCache = _parentStep->getRunParams()->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE");
    
        _pointPrecisionFull = _parentStep->getPbParams()->getAttributeValue<NOMAD::ArrayOfDouble>("POINT_FORMAT");
    }

}


bool NOMAD::IterationUtils::snapPointToBoundsAndProjectOnMesh(
                                NOMAD::EvalPoint& evalPoint,
                                const NOMAD::ArrayOfDouble& lowerBound,
                                const NOMAD::ArrayOfDouble& upperBound)
{
    const NOMAD::EvalPoint evalPoint0 = evalPoint; // Remember first value in case snap does not work.
    NOMAD::Point point = *evalPoint.getX(); // Working locally on point only.

    // Compute fixedVariable
    NOMAD::Point fixedVariable(evalPoint.size());
    // Try/catch ensures that method getSubFixedVariable
    // does not throw an exception.
  
      
    try
    {
        fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(_parentStep);
    }
    catch (NOMAD::Exception&)
    {
        if (nullptr != evalPoint.getPointFrom())
        {
            fixedVariable.resize(evalPoint.getPointFrom()->size());
        }
    }

    if (nullptr == _iterAncestor)
    {
        // If no iterAncestor. Snap the points and the corresponding direction to the bounds
        point.snapToBounds(lowerBound, upperBound);
    }
    else
    {
        // Case with iterAncestor.
        auto mesh = _iterAncestor->getMesh();

        // No mesh --> just snap to bounds
        if (nullptr == mesh)
        {
            point.snapToBounds(lowerBound, upperBound);
        }
        else
        {
            auto center = evalPoint.getPointFrom(fixedVariable);
            if (nullptr == center)
            {
                throw NOMAD::StepException(__FILE__, __LINE__, "snapPointToBoundsAndProjectOnMesh needs a frame center", _parentStep);
            }
            
            #ifdef USE_IBEX
            point = projectWithIbex(point);
            #endif
            
            // First, project on mesh.
            if (_projectOnMesh)
            {
                point = mesh->projectOnMesh(point, *center);
            }
            
            // Second, snap to bounds.
            point.snapToBounds(lowerBound, upperBound);
            
        }
    }

    // Round to POINT_FORMAT number of decimals
    NOMAD::Point pointPrecision = _pointPrecisionFull.projectPointToSubspace(fixedVariable);
    bool modif = point.roundToPrecision(pointPrecision);
    
    if (modif || *evalPoint0.getX() != point )
    {
        // Point is not the same.
        // Update evalPoint
        evalPoint = NOMAD::EvalPoint(point);
        evalPoint.setPointFrom(evalPoint0.getPointFrom(), fixedVariable);
        evalPoint.setGenSteps(evalPoint0.getGenSteps());
        evalPoint.setTag(-1);
    }

    OUTPUT_DEBUG_START
    std::string s = "Point before projection: " + evalPoint0.getX()->display();
    _parentStep->AddOutputDebug(s);
    s = "Point after projection:  " + point.display();
    _parentStep->AddOutputDebug(s);
    OUTPUT_DEBUG_END

    return true;
}


// Sanity check.
// Verify that all points in trialPoints are on the current mesh.
// If a point is not on the mesh -> exception.
bool NOMAD::IterationUtils::verifyPointsAreOnMesh(const std::string& name) const
{
    auto mesh = _iterAncestor->getMesh();
    std::string err;

    if (nullptr == mesh)
    {
        err = "No mesh on iteration (point generated by " + name + ")";
        throw NOMAD::StepException(__FILE__,__LINE__,err, _parentStep);
    }

    for (auto point : _trialPoints)
    {
        auto meshCenter = *point.getPointFrom();
        if (point.size() < meshCenter.size())
        {
            auto fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(_parentStep);
            meshCenter = meshCenter.makeSubSpacePointFromFixed(fixedVariable);
        }
        
        if (!mesh->verifyPointIsOnMesh(point, meshCenter))
        {
            // Let the algorithm decide what to do with that. For example, for NM it is harder to project points on mesh and verify that they are on mesh.
            return false;
        }
    }
    return true;
}

bool NOMAD::IterationUtils::evalTrialPoints(const NOMAD::Step *step,
                                            const size_t keepN,
                                            NOMAD::StepType removeStepType)
{
    bool foundBetter = false;

    // Put the trial points into the evaluation queue
    keepTrialPointsThatNeedEval(step, keepN, removeStepType);
    
    // Send trial EvalPoints to EvaluatorControl
    NOMAD::EvcInterface evcInterface(step);
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    if (_nbEvalPointsThatNeedEval > 0)
    {
                
        _trialPointsSuccess = evcInterface.startEvaluation();

        if (_trialPointsSuccess >= NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            foundBetter = true;
        }

        // Update trial points with evaluated trial points.
        // Note: If cache is not used, Points that are not evaluated yet
        // will be forgotten.
        NOMAD::EvalPointSet evalPointSet;
        for (const auto & evalPoint : evcInterface.retrieveAllEvaluatedPoints())
        {
            evalPointSet.insert(evalPoint);
        }
        OUTPUT_DEBUG_START
        std::string s;
        s = "Number of trial points: " + std::to_string(_trialPoints.size());
        _parentStep->AddOutputDebug(s);
        s = "Number of trial points that needed eval: " + std::to_string(_nbEvalPointsThatNeedEval);
        _parentStep->AddOutputDebug(s);
        s = "Number of evaluated points: " + std::to_string(evalPointSet.size());
        _parentStep->AddOutputDebug(s);
        OUTPUT_DEBUG_END

        _trialPoints.clear();
        _trialPoints = evalPointSet;
    }
    else
    {
        // No new evaluation, clear trial point list.
        _trialPoints.clear();
        
        // No trial point produced
        _trialPointsSuccess = NOMAD::SuccessType::NO_TRIALS;
    }
    
    // Propagate trial points success type to generating method step (for example, poll and search)
    Step* genMethod = const_cast<Step*>(step);
    genMethod->setSuccessType(_trialPointsSuccess);
    
    // Update step success stats from evc success stats
    updateStepSuccessStats(step);
    
    return foundBetter;
}

void NOMAD::IterationUtils::keepTrialPointsThatNeedEval(const Step *step,
                                                        const size_t keepN,
                                                        NOMAD::StepType removeStepType)
{
    // Send trial EvalPoints to EvaluatorControl
    NOMAD::EvcInterface evcInterface(step);
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    // Pass the barrier to the evaluator control for detecting success.
    evcInterface.setBarrier(step->getMegaIterationBarrier());
    
    // Queue will be unlocked one points that need eval are put in the queue
    // Note: lock without unlocking first jams the thread
    evc->lockQueue();

    // If we have a mesh, the evc interface can add some additional information
    bool useMesh = ! _fromAlgo;
    
    // Put trial points that need eval in the evaluation queue
    evcInterface.keepPointsThatNeedEval(_trialPoints, useMesh);

    // The number of points in the queue for evaluation
    // Watch out, this number must be fetched before unlocking the queue,
    // otherwise there is a risk that the evaluations already restarted
    // and the queue may be empty.
    
    // Return the number of eval points in the evaluation queue for the current main thread.
    // These trial points should be of the same type.
    // NOTE: Added trial points are registered by thread number in main thread info.
    _nbEvalPointsThatNeedEval = evc->getQueueSize(NOMAD::getThreadNum());
    
    // Arguments to unlockQueue:
    // true: do sort
    // keepN: keep a maximum of N points
    // removeStepType: only remove points of this StepType.
    evc->unlockQueue(true, keepN, removeStepType);

}

void NOMAD::IterationUtils::countTrialPointsThatNeedEval(const Step *step)
{
    // Send trial EvalPoints to EvaluatorControl
    NOMAD::EvcInterface evcInterface(step);
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
        
    evcInterface.setBarrier(step->getMegaIterationBarrier());
    
    _nbEvalPointsThatNeedEval = evcInterface.countPointsThatNeedEval(_trialPoints);
}


// Post-processing of the points after evaluation.
// Update the barrier: If generating method is a success, update hMax and incumbents of the Barrier.
// Add some stats
bool NOMAD::IterationUtils::postProcessing()
{
    
    const NOMAD::Step * step  = dynamic_cast<NOMAD::Step*>(this);
    if (nullptr != step)
    {
        bool stop=false;     // should be initialized to false or may lead to strange behaviour if empty callback
        step->runCallback(NOMAD::CallbackType::POSTPROCESSING_CHECK, *step, stop);
        
        // Convert CUSTOM_OPPORTUNISTIC_ITER_STOP (evc) into USER_ITER_STOP (iter)
        updateStopReasonForIterStop(step);
        
        // Do we have a global stop?
        // This can only be from a custom callback because default callback returns stop=false;
        if (!step->getAllStopReasons()->checkTerminate() && stop)
        {
            step->getAllStopReasons()->set(NOMAD::BaseStopType::USER_GLOBAL_STOP);
        }
    }
    else
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"An instance of class IterationUtils must also be a step");
    }
    
    // No post processing required when no trial points available
    if ( _trialPoints.size() == 0 )
    {
        return false;
    }
    
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    auto evalType = NOMAD::EvalType::BB;
    auto computeType = NOMAD::ComputeType::STANDARD;
    if (nullptr != evc)
    {
        evalType = evc->getCurrentEvalType();
        computeType = evc->getComputeType();
    }

    bool changeOccured = false;
    auto megaIterBarrier = _megaIterAncestor->getBarrier();

    // The post processing is done when a MegaIteration is used
    if ( megaIterBarrier == nullptr )
    {
        return false;
    }

    // Current hMax of the barrier.
    auto hMax = megaIterBarrier->getHMax();
    const auto hMaxRef = hMax;
    
    // Update incumbents and hMax only according to flag
    bool barrierModified = false;
    if (nullptr != _megaIterAncestor)
    {
        // Make a vector from the set _trialPoints
        std::vector<NOMAD::EvalPoint> evalPointList;
        std::copy(_trialPoints.begin(), _trialPoints.end(),
                  std::back_inserter(evalPointList));
        
         barrierModified = megaIterBarrier->updateWithPoints(
                                evalPointList,
                                evalType,
                                computeType,
                                _frameCenterUseCache /* not used by progressive barrier */,
                                _updateIncumbentsAndHMax /* set by trial point generating method */);
        
        hMax = _megaIterAncestor->getBarrier()->getHMax();
    }
    

    
    // Update hMax
    if (hMax < hMaxRef)
    {
        OUTPUT_DEBUG_START
        _parentStep->AddOutputDebug("hMax went from " + hMaxRef.tostring() + " to " + hMax.tostring());
        OUTPUT_DEBUG_END
        _parentStep->getMegaIterationBarrier()->setHMax(hMax);
        changeOccured = true;
    }
    changeOccured = changeOccured || barrierModified;

    NOMAD::OutputQueue::Flush();
    
    // Update local stats with evaluations done
    size_t nbTrialPointsEvaluated = 0;
    for (const auto & trialPoint : _trialPoints)
    {
        // We are looking for trial points that have been evaluated
        if (trialPoint.isEvalOk(evalType))
        {
            nbTrialPointsEvaluated ++;
        }
    }

    _trialPointStats.incrementEvalsDone(nbTrialPointsEvaluated, evalType);
    _trialPointStats.updateParentStats();

    return changeOccured;
}
// End postProcessing

void NOMAD::IterationUtils::updateStats(TrialPointStats &trialPointStats)
{
    _trialPointStats.updateWithCurrentStats(trialPointStats);
}

void NOMAD::IterationUtils::updateStepSuccessStats(const Step* step)
{
    
    
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    NOMAD::EvalType evalType;
    if (nullptr != evc)
    {
        evalType = evc->getCurrentEvalType();
    }
    else
    {
        return;
    }
    
    // Update step success stats with evc success stats
    // Important: For steps directly generating/evaluating trial points, each evaluated trial point is counted (UNSUCCESSFUL, PARTIAL_SUCCESS, FULL_SUCCESS). Unevaluated trial points are not counted.
    if (NOMAD::EvalType::BB == evalType)
    {
        const SuccessStats & evcSuccessStats= evc->getSuccessStats();
        
        if (evcSuccessStats.hasStatsForPropagation())
        {
            Step* stepToUpdate = const_cast<Step*>(step);
            NOMAD::SuccessStats & stepStats = stepToUpdate->getSuccessStats();
            stepStats.updateStats(evcSuccessStats); // Update the stats of current step
        }
        
        // Each evaluated trial point should only be counted once. Reset evaluator control stats after transfer.
        evc->resetSuccessStats();
    }
}

bool NOMAD::IterationUtils::insertTrialPoint(const NOMAD::EvalPoint &evalPoint)
{
    // We do not need the Eval part of EvalPoint right now,
    // but it will be used soon. Could be refactored, but
    // not high priority. Note that an EvalPointSet compares
    // the Point part of the EvalPoints only.

    //Set the eval point tag and increment for the next point.
    evalPoint.updateTag();

    std::pair<NOMAD::EvalPointSet::iterator,bool> ret = _trialPoints.insert(evalPoint);

    OUTPUT_INFO_START
    std::string s = "Point:";
    s += (ret.second) ? " " : " not inserted: ";
    s += evalPoint.display();
    NOMAD::OutputInfo("",s,NOMAD::OutputLevel::LEVEL_INFO);
    OUTPUT_INFO_END

    return ret.second ;
}

void NOMAD::IterationUtils::generateTrialPoints()
{
    clearTrialPoints();
     
    // Reset the trial point stats (this may have been done before)
    _trialPointStats.resetCurrentStats();
    
    // Call implementation to generate trial points
    generateTrialPointsImp();
    
    // Update success type if no points are generated
    if (_trialPoints.size() == 0)
    {
        _trialPointsSuccess = NOMAD::SuccessType::NO_TRIALS;
    }
    
    // Update counters of generated trial points
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    auto evalType = NOMAD::EvalType::BB;
    if (nullptr != evc)
    {
        evalType = evc->getCurrentEvalType();
    }
    _trialPointStats.incrementTrialPointsGenerated(_trialPoints.size(), evalType);
}

void NOMAD::IterationUtils::generateTrialPointsSecondPass()
{
    
    // Call implementation to generate trial points for second pass
    generateTrialPointsSecondPassImp();
    
    // Update counters of generated trial points
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    auto evalType = NOMAD::EvalType::BB;
    if (nullptr != evc)
    {
        evalType = evc->getCurrentEvalType();
    }
    _trialPointStats.incrementTrialPointsGenerated(_trialPoints.size(), evalType);
}



bool NOMAD::IterationUtils::meshIsFinest() const
{
    if (nullptr == _iterAncestor)
    {
        return false;
    }
    else
    {
        // Case with iterAncestor.
        auto mesh = _iterAncestor->getMesh();

        // No mesh --> just snap to bounds
        if (nullptr == mesh)
        {
            throw NOMAD::StepException(__FILE__, __LINE__,
                            "An instance of class IterationUtils call meshIsFinest must have a mesh",
                            _parentStep);
        }
        return mesh->isFinest();
    }
}


// Complete trial points information for sorting before evaluation
void NOMAD::IterationUtils::completeTrialPointsInformation()
{
    
    // Send trial EvalPoints to EvaluatorControl
    NOMAD::EvcInterface evcInterface(_parentStep);
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    
    std::unique_ptr<NOMAD::SurrogateEvaluation> surrogateEvaluation = nullptr;

    // If sort type is MODEL, but Evaluator type is not MODEL,
    // start by evaluating points using the surrogate (MODEL) Evaluator.
    // No need to sort if not opportunistic
    if ( NOMAD::EvalSortType::QUADRATIC_MODEL == evc->getEvalSortType()
        && NOMAD::EvalType::MODEL != evc->getCurrentEvalType()
        && _trialPoints.size() > 1
        && evc->getOpportunisticEval())
    {
        // Reset the counter (otherwise the cumulative model evals for sorting may exceed the limit MODEL_MAX_EVAL)
        evc->resetModelEval();
        
        // Construction of quadratic model
        surrogateEvaluation = std::make_unique<NOMAD::SurrogateEvaluation>(_parentStep, _trialPoints, NOMAD::EvalType::MODEL);
    }
    // If sort type is SURROGATE, but Evaluator type is not SURROGATE,
    // start by evaluating points using the surrogate Evaluator.
    else if ( NOMAD::EvalSortType::SURROGATE == evc->getEvalSortType()
        && NOMAD::EvalType::SURROGATE != evc->getCurrentEvalType()
        && _trialPoints.size() > 1
        && evc->getOpportunisticEval())
    {
        surrogateEvaluation = std::make_unique<NOMAD::SurrogateEvaluation>(_parentStep,_trialPoints, NOMAD::EvalType::SURROGATE);
    }
    
    if (nullptr != surrogateEvaluation)
    {
        surrogateEvaluation->start(); // start sets the eval type to MODEL or SURROGATE, perform MODEL construction if it is its eval type
        surrogateEvaluation->run(); // Perform MODEL or SURROGATE evaluations on the trial points
        surrogateEvaluation->end();
    }
}

#ifdef USE_IBEX
NOMAD::Point NOMAD::IterationUtils::projectWithIbex(NOMAD::Point point)
{
	auto mainStepAncestorConst = _parentStep->getParentOfType<NOMAD::MainStep*>(false);
	auto mainStepAncestor = const_cast<NOMAD::MainStep*>(mainStepAncestorConst);
	auto set = mainStepAncestor->getIbexSet();
	
    size_t n = point.size();
        
    ibex::Vector v(n);
    for (size_t i = 0; i < n; i++)
    {
    	v[i] = point[i].trunk();
    }
    
    ibex::Vector v_projected = (*set).move_inside(v);
    
    for (size_t i = 0; i < n; i++)
    {
    	point[i] = v_projected[i];
    }

    return point;
}
#endif

void NOMAD::IterationUtils::updateStopReasonForIterStop(const Step* step)
{
    // Test for NOMAD::EvalMainThreadStopType::CUSTOM_OPPORTUNISTIC_ITER_STOP)
    // and transform into NOMAD::IterStopType::USER_ITER_STOP
    
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    // This is postprocessing for BB only
    if (NOMAD::EvalType::BB != evc->getCurrentEvalType())
    {
        return;
    }
    auto evcStopReason = evc->getStopReason(-1);
    
    if (evcStopReason.checkStopType(NOMAD::EvalMainThreadStopType::CUSTOM_OPPORTUNISTIC_ITER_STOP))
    {
        // Reset evcStopReason
        evc->setStopReason(-1, NOMAD::EvalMainThreadStopType::STARTED);
        
        // Only replace the default iter stop reason (no stop)
        if (step->getAllStopReasons()->testIf(NOMAD::IterStopType::STARTED))
        {
            step->getAllStopReasons()->set(NOMAD::IterStopType::USER_ITER_STOP);
            
            OUTPUT_INFO_START
            NOMAD::OutputQueue::Add("User iter stop in "+step->getName(), NOMAD::OutputLevel::LEVEL_INFO);
            NOMAD::OutputQueue::Flush();
            OUTPUT_INFO_END
        }
        
    }
}

