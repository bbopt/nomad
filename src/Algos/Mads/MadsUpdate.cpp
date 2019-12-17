/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Output/OutputInfo.hpp"

void NOMAD::MadsUpdate::init()
{
    _name = getAlgoName() + "Update";
    verifyParentNotNull();

    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>(getParentOfType<NOMAD::MadsMegaIteration*>());
    if (nullptr == megaIter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"An instance of class MadsUpdate must have a MegaIteration among its ancestors");
    }

}


bool NOMAD::MadsUpdate::runImp()
{
    // megaIter barrier is already in subproblem.
    // So no need to convert refBestFeas and refBestInf
    // from full dimension to subproblem.
    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>(getParentOfType<NOMAD::MadsMegaIteration*>());
    auto barrier = megaIter->getBarrier();
    auto mesh = megaIter->getMesh();

    if (NOMAD::CacheBase::getInstance()->size() > 1)
    {
        // Remember previous values
        auto refBestFeas = barrier->getFirstXFeas();
        auto refBestInf  = barrier->getFirstXInf();

        // Update Best feasible and infeasible with cache.
        updateBestInBarrier(true);
        updateBestInBarrier(false);

        // Compute success
        NOMAD::EvalPointPtr newBestFeas = barrier->getFirstXFeas();
        NOMAD::EvalPointPtr newBestInf  = barrier->getFirstXInf();

        std::shared_ptr<NOMAD::EvalPoint> newBest;

        // Get which of newBestFeas and newBestInf is improving
        // the solution. Check newBestFeas first.
        NOMAD::ComputeSuccessType computeSuccess;
        NOMAD::SuccessType success = computeSuccess(newBestFeas, refBestFeas);
        if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            // newBestFeas is the improving point.
            newBest = newBestFeas;
            // Output Warning: When using '\n', the computed indentation for the
            // Step will be ignored. Leaving it like this for now. Using an
            // OutputInfo with AddMsg() would resolve the output layout.
            std::string s = "Update: improving feasible point";
            if (refBestFeas)
            {
                s += " from\n    " + refBestFeas->display() + "\n";
            }
            s += " to " + newBestFeas->display();
            AddOutputDebug(s);
        }
        else
        {
            // Check newBestInf
            NOMAD::SuccessType success2 = computeSuccess(newBestInf, refBestInf);
            if (success2 > success)
            {
                success = success2;
            }
            if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                // newBestInf is the improving point.
                newBest = newBestInf;
                std::string s = "Update: improving infeasible point";
                if (refBestInf)
                {
                    s+= " from\n    " + refBestInf->display() + "\n";
                }
                s += " to " + newBestInf->display();
                AddOutputDebug(s);
            }
        }
        if (success == NOMAD::SuccessType::UNSUCCESSFUL)
        {
            std::string s = "Update: no success found";
            AddOutputDebug(s);
        }

        // NOTE enlarge or refine might have to be done multiple times
        // in a row, if we are working on multiple meshes at the same
        // time.
        // If we found a better refBest but only for a mesh that is
        // already refined a few times, we must refine mainMesh accordingly
        // before enlarging it.
        // If we found a better refBest for a mesh that is enlarged
        // a few times, we must enlarge mainMesh this many times.
        // If we did not find a better refBest, we refine the mesh
        // once, we might fall on points that are already evaluated - which
        // is fine, or on new points - which will be evaluated.
        //
        // For now, leave it as is, since multiple mesh sizes are not
        // currently evaluated.
        //
        // TODO Analyze, write tests, and implement.


        // Debug verification
        // Compare computed success with value from MegaIteration.
        // This is the value from the previous MegaIteration. If it
        // was not evaluated, ignore the test.
        // It is possible that the MegaIteration found a partial success,
        // and then the EvaluatorControl found a full success before 
        // Update is run.
        // For this reason, only test the boolean value success vs. failure.
        const bool megaIterSuccessful = (megaIter->getSuccessType() >= NOMAD::SuccessType::PARTIAL_SUCCESS);
        const bool successful = (success >= NOMAD::SuccessType::PARTIAL_SUCCESS);
        if (   (NOMAD::SuccessType::NOT_EVALUATED != megaIter->getSuccessType())
            && (   (successful != megaIterSuccessful)
                || (NOMAD::SuccessType::NOT_EVALUATED == success)) )
        {
            std::string s = "Warning: MegaIteration success type: ";
            s += NOMAD::enumStr(megaIter->getSuccessType());
            s += ". Is different than computed success type: " + NOMAD::enumStr(success);
            if (refBestFeas)
            {
                s += "\nRef best feasible:   " + refBestFeas->display();
            }
            if (newBestFeas)
            {
                s += "\nNew best feasible:   " + newBestFeas->display();
            }
            if (refBestInf)
            {
                s += "\nRef best infeasible: " + refBestInf->display();
            }
            if (newBestInf)
            {
                s += "\nNew best infeasible: " + newBestInf->display();
            }
            AddOutputWarning(s);
        }

        if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            // Compute new direction for main mesh.
            // The direction is related to the frame center which generated
            // newBest.
            auto pointFromPtr = newBest->getPointFrom(); 
            auto pointNewPtr = newBest->getX();
            if (nullptr == pointFromPtr)
            {
                std::string s = "Update cannot compute new direction for successful point: pointFromPtr is NULL ";
                s += newBest->display();
                throw NOMAD::Exception(__FILE__,__LINE__, s);
            }
            // PointFrom is in full dimension. Convert it to subproblem
            // to compute direction.
            auto fixedVariable = _parentStep->getSubFixedVariable();
            auto pointFromSub = std::make_shared<NOMAD::Point>(pointFromPtr->makeSubSpacePointFromFixed(fixedVariable));

            NOMAD::Direction dir = NOMAD::Point::vectorize(*pointFromSub, *pointNewPtr);
            std::string dirStr = "New direction " + dir.display();
            AddOutputInfo(dirStr);

            // Update frame size for main mesh
            AddOutputInfo("Last Iteration Successful.");
            auto anisotropyFactor = _runParams->getAttributeValue<NOMAD::Double>("ANISOTROPY_FACTOR");
            bool anistropicMesh = _runParams->getAttributeValue<bool>("ANISOTROPIC_MESH");

            if (mesh->enlargeDeltaFrameSize(dir, anisotropyFactor, anistropicMesh))
            {
                AddOutputInfo("Delta is enlarged.");
            }
            else
            {
                AddOutputInfo("Delta is not enlarged.");
            }
        }
        else
        {
            AddOutputInfo("Last Iteration Unsuccessful. Refine Delta.");
            mesh->refineDeltaFrameSize();
        }
    }

    mesh->updatedeltaMeshSize();
    AddOutputInfo("delta mesh  size = " + mesh->getdeltaMeshSize().display());
    AddOutputInfo("Delta frame size = " + mesh->getDeltaFrameSize().display());


    return true;
}


void NOMAD::MadsUpdate::updateBestInBarrier(bool feasible)
{
    std::vector<NOMAD::EvalPoint> evalPointList;
    auto barrier = getMegaIterationBarrier();
    auto evalType = getEvalType();
    std::string s;  // For output info

    NOMAD::Double hMax = 0.0;
    if (!feasible)
    {
        hMax = barrier->getHMax();
    }
    // Find best points in cache.
    NOMAD::CacheInterface cacheInterface(this);
    if (feasible)
    {
        cacheInterface.findBestFeas(evalPointList, getSubFixedVariable(), evalType);
    }
    else
    {
        cacheInterface.findBestInf(evalPointList, hMax, getSubFixedVariable(), evalType);
    }

    // Debug info
    s = "Found " + std::to_string(evalPointList.size());
    s += (feasible) ? " feasible " : " infeasible ";
    s += "best points";
    NOMAD::OutputInfo outputInfo(_name, s, NOMAD::OutputLevel::LEVEL_DEBUG);
    size_t nbPoints = 0;
    for (auto evalPoint : evalPointList)
    {
        outputInfo.addMsg(evalPoint.display());
        nbPoints++;
        if (nbPoints >=4)
        {
            outputInfo.addMsg("...");
            break;
        }
    }
    NOMAD::OutputQueue::Add(std::move(outputInfo));

    if (_runParams->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"))
    {
        // Clear barrier and add newly found points.
        // Note: evalPoint for MegaIteration Barrier are in subspace - not full dimension.
        // The points in evalPointList are also in subspace, thanks to the cacheInterface.
        if (feasible)
        {
            barrier->clearXFeas();
            for (auto evalPoint : evalPointList)
            {
                barrier->addXFeas(std::make_shared<NOMAD::EvalPoint>(evalPoint), evalType);
            }
        }
        else
        {
            barrier->clearXInf();
            for (auto evalPoint : evalPointList)
            {
                barrier->addXInf(std::make_shared<NOMAD::EvalPoint>(evalPoint));
            }
        }
    }

    else
    {
        // If the found points have the same value as the one in the barrier,
        // we want to keep the one that is already in the barrier.
        if (evalPointList.size() > 0)
        {
            NOMAD::EvalPointPtr newBest = std::make_shared<NOMAD::EvalPoint>(evalPointList[0]);
            NOMAD::EvalPointPtr refBest = (feasible) ? barrier->getFirstXFeas() : barrier->getFirstXInf();
            NOMAD::ComputeSuccessType computeSuccess;
            bool doUpdateBarrier = false;
            if (nullptr == refBest)
            {
                // No previous point. Update the barrier.
                doUpdateBarrier = true;
            }
            else if (computeSuccess(newBest, refBest) >= NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                // Dominance. Update the barrier.
                doUpdateBarrier = true;
            }
            else if (computeSuccess(refBest, newBest) == NOMAD::SuccessType::UNSUCCESSFUL)
            {
                // No point dominates the other. Do not update the Barrier - unless it has too many 
                // points (due to initialization).
                nbPoints = (feasible) ? barrier->getAllXFeas().size() : barrier->getAllXInf().size();
                if (nbPoints > 1)
                {
                    doUpdateBarrier = true;
                }
            }
            else
            {
                // If we end up here, it means that refBest dominates newBest.
                // There is a problem.
                s = "Warning: Previous ";
                s += ((feasible) ? "feasible " : "infeasible ");
                s += "frame center is better than current.";
                s += "\nPrevious frame center:   " + refBest->display();
                s += "\nCurrent  frame center:   " + newBest->display();
                AddOutputWarning(s);
            }

            if (doUpdateBarrier)
            {
                if (feasible)
                {
                    barrier->clearXFeas();
                    barrier->addXFeas(newBest, evalType);
                }
                else
                {
                    barrier->clearXInf();
                    barrier->addXInf(newBest);
                }
            }
        }
    }
}

