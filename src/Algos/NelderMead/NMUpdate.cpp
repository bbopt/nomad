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
#include "../../Algos/NelderMead/NMIteration.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"
#include "../../Algos/NelderMead/NMUpdate.hpp"
#include "../../Output/OutputInfo.hpp"

void NOMAD::NMUpdate::init()
{
    _name = getAlgoName() + "Update";
    verifyParentNotNull();

}


bool NOMAD::NMUpdate::runImp()
{

    const NOMAD::NMIteration * iter = dynamic_cast<const NOMAD::NMIteration*>(getParentOfType<NOMAD::NMIteration*>());
    if (nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Update must have a NMIteration among its ancestors.");
    }

    if ( iter->getY()->size() == 0 )
    {
        AddOutputDebug("No need to update because NM simplex is empty");
        return true;
    }

    const NOMAD::NMMegaIteration * megaIter = dynamic_cast<const NOMAD::NMMegaIteration*>(getParentOfType<NOMAD::NMMegaIteration*>());
    
    // megaIter barrier is already in subproblem.
    // So no need to convert refBestFeas and refBestInf
    // from full dimension to subproblem.
    auto barrier = megaIter->getBarrier();

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
            AddOutputDebug("Update: improving feasible point");
            if (refBestFeas)
            {
                AddOutputDebug(" from    " + refBestFeas->display() );
            }
            AddOutputDebug( " to " + newBestFeas->display() );
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
                AddOutputDebug("Update: improving infeasible point");
                if (refBestInf)
                {
                    AddOutputDebug(" from    " + refBestInf->display() );
                }
                AddOutputDebug(" to " + newBestInf->display());

            }
        }
        if (success == NOMAD::SuccessType::UNSUCCESSFUL)
        {
            std::string s = "Update: no success found";
            AddOutputDebug(s);
        }

        if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            // Compute new direction for main mesh.
            // The direction is related to the frame center which generated
            // newBest.
            auto pointFromPtr = newBest->getPointFrom();
            // TODO: Make getX() return a shared_ptr!!
            auto pointNewPtr  = newBest->getX();
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

            AddOutputInfo("Last NM step Successful.");
        }
        else
        {
            AddOutputInfo("Last NM step Unsuccessful.");
        }
    }

    return true;
}


void NOMAD::NMUpdate::updateBestInBarrier(bool feasible)
{
    std::vector<NOMAD::EvalPoint> evalPointList;
    auto barrier = getMegaIterationBarrier();
    auto evalType = getEvalType();
    std::string s;

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
                s = "Warning: Previous frame center is better than current.";
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

