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
            //auto megaIterAncestorConst = _parentStep->getParentOfType<NOMAD::MegaIteration*>();
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
                                NOMAD::Point& point,
                                const NOMAD::ArrayOfDouble& lowerBound,
                                const NOMAD::ArrayOfDouble& upperBound)
{
    bool snapWorked = true;
    const NOMAD::Point point0 = point; // Remember first value in case snap does not work.

    if ( nullptr == _iterAncestor )
    {
        // If no iterAncestor. Snap the points and the corresponding direction to the bounds
        if (!point.inBounds(lowerBound, upperBound))
        {
            point.snapToBounds(lowerBound, upperBound, NOMAD::ArrayOfDouble(), NOMAD::ArrayOfDouble());
            if (!point.inBounds(lowerBound, upperBound))
            {
                snapWorked = false;
            }
        }
    }
    else
    {
        // Case with iterAncestor.

        auto center = _iterAncestor->getFrameCenter();
        auto mesh = _iterAncestor->getMesh();

        // No mesh --> just snap to bounds
        if (nullptr == mesh)
        {
            if (!point.inBounds(lowerBound, upperBound))
            {
                point.snapToBounds(lowerBound, upperBound, *center, NOMAD::ArrayOfDouble());
                if (!point.inBounds(lowerBound, upperBound))
                {
                    snapWorked = false;
                }
            }
        }
        else
        {
            bool firstProjectionWorked = false;

            // These points are for debug info only.
            NOMAD::Point pointBefore = point;
            NOMAD::Point pointFirstProj = point;
            NOMAD::Point pointSnap = point;

            // First, project on mesh.
            if (!mesh->verifyPointIsOnMesh(point, *center))
            {
                point = mesh->projectOnMesh(point, *center);
                if (mesh->verifyPointIsOnMesh(point, *center))
                {
                    firstProjectionWorked = true;
                }
            }
            pointFirstProj = point;
            // Second, snap to bounds.
            if (!point.inBounds(lowerBound, upperBound))
            {
                point.snapToBounds(lowerBound, upperBound, *center, mesh->getdeltaMeshSize());
            }
            pointSnap = point;
            // Third, if needed, project on mesh again.
            if (!mesh->verifyPointIsOnMesh(point, *center))
            {
                point = mesh->projectOnMesh(point, *center);
            }

            if (!point.inBounds(lowerBound, upperBound) || !mesh->verifyPointIsOnMesh(point, *center))
            {
                snapWorked = false;

                // Debug info
                OUTPUT_DEBUG_START
                if (firstProjectionWorked)
                {
                    NOMAD::OutputInfo outputInfo("Snap", "Warning: point was not snapped properly on mesh:", NOMAD::OutputLevel::LEVEL_DEBUG);
                    // First projection worked, but then the snapToBounds offset it from mesh.
                    outputInfo.addMsg("Point before projection:");
                    NOMAD::ArrayOfDouble debugPrecision(point.size(), 20);
                    outputInfo.addMsg(pointBefore.display(debugPrecision));
                    outputInfo.addMsg("Point after first projection:");
                    outputInfo.addMsg(pointFirstProj.display(debugPrecision));
                    outputInfo.addMsg("Point after snapping to bounds:");
                    outputInfo.addMsg(pointSnap.display(debugPrecision));
                    outputInfo.addMsg("Point after second projection:");
                    outputInfo.addMsg(point.display(debugPrecision));
                    outputInfo.addMsg("Lower bound: " + lowerBound.display(debugPrecision));
                    outputInfo.addMsg("Upper bound: " + upperBound.display(debugPrecision));
                    outputInfo.addMsg("Center: " + center->display(debugPrecision));
                    outputInfo.addMsg("Mesh size: " + mesh->getdeltaMeshSize().display(debugPrecision));
                    NOMAD::OutputQueue::Add(std::move(outputInfo));
                }
                OUTPUT_DEBUG_END
            }
        }
    }

    if (!snapWorked)
    {
        // Revert point to first value
        point = point0;
    }
    else if (point0 != point)
    {
        // In the case we are working on an EvalPoint,
        // Point is not the same. Tag will be updated.
        auto evalPoint = dynamic_cast<NOMAD::EvalPoint*>(&point);
        if (nullptr != evalPoint)
        {
            evalPoint->setTag(0);
        }
    }

    return snapWorked;
}


// Sanity check.
// Verify that all points in trialPoints are on the current mesh.
// If a point is not on the mesh -> exception.
void NOMAD::IterationUtils::verifyPointsAreOnMesh(const std::string& name) const
{
    auto mesh = _iterAncestor->getMesh();
    auto frameCenter = _iterAncestor->getFrameCenter();
    std::string err;

    if (nullptr == frameCenter)
    {
        err = "No frame center on iteration (point generated by " + name + ")";
        throw NOMAD::StepException(__FILE__,__LINE__,err, _parentStep);
    }
    if (nullptr == mesh)
    {
        err = "No mesh on iteration (point generated by " + name + ")";
        throw NOMAD::StepException(__FILE__,__LINE__,err, _parentStep);
    }

    for (auto point : _trialPoints)
    {
        if (!mesh->verifyPointIsOnMesh(point, *frameCenter))
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

    // When doSort is false, lexicographical order is used.
    auto disable = step->getRunParams()->getAttributeValue<NOMAD::ArrayOfString>("DISABLE");
    bool doSort = (-1 == disable.find("EVAL_SORT"));
    NOMAD::EvcInterface::getEvaluatorControl()->unlockQueue(doSort);

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
bool NOMAD::IterationUtils::postProcessing(const NOMAD::EvalType& evalType)
{
    bool changeOccured = false;
    auto megaIterBarrier = _megaIterAncestor->getBarrier();

    // The post processing is done when a MegaIteration is used
    if ( megaIterBarrier == nullptr )
    {
        return false;
    }

    auto xInf = megaIterBarrier->getFirstXInf();
    NOMAD::Double fxInf, hxInf;
    if (nullptr != xInf)
    {
        fxInf = xInf->getF(evalType);
        hxInf = xInf->getH(evalType);
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
            if (trialPoint.isFeasible(evalType))
            {
                // We are only interested in infeasible points, i.e., h > 0.
                continue;
            }

            NOMAD::Double ftrialPoint = trialPoint.getF(evalType);
            NOMAD::Double htrialPoint = trialPoint.getH(evalType);

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

    // Update hMax
    if (hMax < hMaxRef)
    {
        OUTPUT_DEBUG_START
        _parentStep->AddOutputDebug("hMax went from " + hMaxRef.tostring() + " to " + hMax.tostring());
        OUTPUT_DEBUG_END
        _parentStep->getMegaIterationBarrier()->setHMax(hMax);
        changeOccured = true;
    }

    // Update Barrier right away.
    if (nullptr != _megaIterAncestor)
    {
        // Make a vector from the set _trialPoints
        std::vector<NOMAD::EvalPoint> evalPointList;
        std::copy(_trialPoints.begin(), _trialPoints.end(),
                  std::back_inserter(evalPointList));
        bool barrierModified = _megaIterAncestor->getBarrier()->updateWithPoints(evalPointList,
                                                              evalType,
                                                              _parentStep->getRunParams()->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"));
        changeOccured = changeOccured || barrierModified;
    }

    NOMAD::OutputQueue::Flush();

    return changeOccured;
}
// End postProcessing


void NOMAD::IterationUtils::updatePointsWithFrameCenter()
{
    if ( nullptr == _iterAncestor )
    {
        OUTPUT_DEBUG_START
        _parentStep->AddOutputDebug("No ancestor, no frame center");
        OUTPUT_DEBUG_END
        return;
    }

    auto frameCenter = _iterAncestor->getFrameCenter();
    if (nullptr == frameCenter)
    {
        OUTPUT_DEBUG_START
        _parentStep->AddOutputDebug("Cannot update point with NULL frame center from iteration.");
        OUTPUT_DEBUG_END
        return ;
    }


    // frameCenter has to be converted to full dimension, to be able to refer
    // to it consistently later.
    auto fixedVariable = NOMAD::SubproblemManager::getSubFixedVariable(_parentStep);
    std::shared_ptr<NOMAD::Point> frameCenterFull = std::make_shared<NOMAD::Point>(frameCenter->getX()->makeFullSpacePointFromFixed(fixedVariable));

    for (auto it = _trialPoints.begin(); it != _trialPoints.end(); it++)
    {
        // Update EvalPoint directly.
        // Since we are not changing the Point part, which is the only part
        // used for sorting, the EvalPointSet remains coherent. This
        // is why we can use the const_cast.
        auto evalPoint = const_cast<NOMAD::EvalPoint*>(&*it);

        evalPoint->setPointFrom(frameCenterFull);

        // Debug info
        OUTPUT_DEBUG_START
        std::string s = "Set pointFrom of point ";
        s += evalPoint->getX()->display();
        s += " to ";
        s += (nullptr == frameCenterFull) ? "NULL" : frameCenterFull->display();
        _parentStep->AddOutputDebug(s);
        OUTPUT_DEBUG_END
    }
}


bool NOMAD::IterationUtils::insertTrialPoint(const NOMAD::EvalPoint &evalPoint)
{
    // We do not need the Eval part of EvalPoint right now,
    // but it will be used soon. Could be refactored, but
    // not high priority. Note that an EvalPointSet compares
    // the Point part of the EvalPoints only.

    std::pair<NOMAD::EvalPointSet::iterator,bool> ret = _trialPoints.insert(evalPoint);

    OUTPUT_INFO_START
    std::string s = "xt:";
    s += (ret.second) ? " " : " not inserted: ";
    s += evalPoint.display();
    NOMAD::OutputInfo("",s,NOMAD::OutputLevel::LEVEL_INFO);
    OUTPUT_INFO_END

    return ret.second ;
}
