/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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
#include "../Algos/SubproblemManager.hpp"
#include "../Output/OutputQueue.hpp"

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
    catch (NOMAD::Exception &/*e*/)
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

            // First, project on mesh.
            point = mesh->projectOnMesh(point, *center);
            // Second, snap to bounds.
            point.snapToBounds(lowerBound, upperBound);
        }
    }

    if (*evalPoint0.getX() != point)
    {
        // Point is not the same.
        // Update evalPoint
        evalPoint = NOMAD::EvalPoint(point);
        evalPoint.setPointFrom(evalPoint0.getPointFrom(), fixedVariable);
        evalPoint.setGenStep(evalPoint0.getGenStep());
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
void NOMAD::IterationUtils::verifyPointsAreOnMesh(const std::string& name) const
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
            err = "This point (generated by " + name + ")";
            err += " is not on the mesh: " + point.display() + ".";
            throw NOMAD::StepException(__FILE__,__LINE__,err, _parentStep);
        }
    }
}


bool NOMAD::IterationUtils::evalTrialPoints(NOMAD::Step *step)
{
    bool foundBetter = false;

    // Send trial EvalPoints to EvaluatorControl
    NOMAD::EvcInterface evcInterface(step);

    NOMAD::EvcInterface::getEvaluatorControl()->lockQueue();

    // If we have a mesh, the evc interface can add some additional information
    bool useMesh = ! _fromAlgo;
    evcInterface.keepPointsThatNeedEval(_trialPoints, useMesh);

    evcInterface.setBarrier(step->getMegaIterationBarrier());

    // The number of points in the queue for evaluation
    // Watch out, this number must be fetched before unlocking the queue,
    // otherwise there is a risk that the evaluations already restarted
    // and the queue may be empty.
    _nbEvalPointsThatNeedEval = NOMAD::EvcInterface::getEvaluatorControl()->getQueueSize(NOMAD::getThreadNum());

    NOMAD::EvcInterface::getEvaluatorControl()->unlockQueue(true);  // true: do sort

    if (_nbEvalPointsThatNeedEval > 0)
    {
        _success = evcInterface.startEvaluation();

        if (_success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            foundBetter = true;
        }

        // Update trial points with evaluated trial points.
        // Note: If cache is not used, Points that are not evaluated yet
        // will be forgotten.
        NOMAD::EvalPointSet evalPointSet;
        for (auto evalPoint : evcInterface.retrieveAllEvaluatedPoints())
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
    }
    return foundBetter;
}


// Post-processing of the points after evaluation.
// For instance, computation of a new hMax and update of the Barrier.
bool NOMAD::IterationUtils::postProcessing()
{
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    auto evalType = NOMAD::EvalType::BB;
    auto computeType = NOMAD::ComputeType::STANDARD;
    if (nullptr != evc)
    {
        evalType = evc->getEvalType();
        computeType = evc->getComputeType();
    }

    bool changeOccured = false;
    auto megaIterBarrier = _megaIterAncestor->getBarrier();

    // The post processing is done when a MegaIteration is used
    if ( megaIterBarrier == nullptr )
    {
        return false;
    }

    auto xInf = megaIterBarrier->getFirstXInf();
    NOMAD::Double fxInf = -NOMAD::INF;
    NOMAD::Double hxInf = megaIterBarrier->getHMax();
    if (nullptr != xInf)
    {
        fxInf = xInf->getF(evalType, computeType);
        hxInf = xInf->getH(evalType, computeType);
    }

    // Current hMax is hMax of the barrier.
    auto hMax = megaIterBarrier->getHMax();
    const auto hMaxRef = hMax;

    // Compute hMax in case of PARTIAL_SUCCESS.
    if (NOMAD::SuccessType::PARTIAL_SUCCESS == _success)
    {
        NOMAD::Double tempHMax;
        // Trial point is already updated with its Eval.
        for (auto trialPoint : _trialPoints)
        {
            // We are looking for improving, non-dominating trial points.
            // I.e. h is better, but f is less good.

            // Note: Searching for updated trial points in the cache.
            if (trialPoint.isFeasible(evalType, computeType))
            {
                // We are only interested in infeasible points, i.e., h > 0.
                continue;
            }

            NOMAD::Double ftrialPoint = trialPoint.getF(evalType, computeType);
            NOMAD::Double htrialPoint = trialPoint.getH(evalType, computeType);

            bool evalOk = (NOMAD::EvalStatusType::EVAL_OK == trialPoint.getEvalStatus(evalType));

            if (evalOk && (htrialPoint < hxInf)
                && (ftrialPoint > fxInf))
            {
                // improving
                if (!tempHMax.isDefined() || tempHMax < htrialPoint)
                {
                    // Keep highest hMax from all improving points.
                    tempHMax = htrialPoint;
                }
            }
        }
        // Failsafe.
        if (!tempHMax.isDefined())
        {
            tempHMax = hxInf;
        }
        hMax = tempHMax;
    }
    // In the case of full success, or unsuccessful, hMax becomes the h of xInf.
    else if (hxInf.isDefined())
    {
        hMax = hxInf;
    }



    // Update Barrier right away.
    bool barrierModified = false;
    if (nullptr != _megaIterAncestor)
    {
        // Make a vector from the set _trialPoints
        std::vector<NOMAD::EvalPoint> evalPointList;
        std::copy(_trialPoints.begin(), _trialPoints.end(),
                  std::back_inserter(evalPointList));
         barrierModified = _megaIterAncestor->getBarrier()->updateWithPoints(evalPointList,
                                evalType,
                                computeType,
                                _parentStep->getRunParams()->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"));
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

    return changeOccured;
}
// End postProcessing


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
    std::string s = "xt:";
    s += (ret.second) ? " " : " not inserted: ";
    s += evalPoint.display();
    NOMAD::OutputInfo("",s,NOMAD::OutputLevel::LEVEL_INFO);
    OUTPUT_INFO_END

    return ret.second ;
}
