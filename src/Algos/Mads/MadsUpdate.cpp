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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::MadsUpdate::init()
{
    _name = getAlgoName() + "Update";
    verifyParentNotNull();

    auto megaIter = getParentOfType<NOMAD::MadsMegaIteration*>();
    if (nullptr == megaIter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Error: An instance of class MadsUpdate must have a MegaIteration among its ancestors");
    }

}


bool NOMAD::MadsUpdate::runImp()
{
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    auto evalType = evc->getEvalType();
    // megaIter barrier is already in subproblem.
    // So no need to convert refBestFeas and refBestInf
    // from full dimension to subproblem.
    auto megaIter = getParentOfType<NOMAD::MadsMegaIteration*>();
    auto barrier = megaIter->getBarrier();
    auto mesh = megaIter->getMesh();
    std::string s;  // for output

    OUTPUT_DEBUG_START
    s = "Running " + getName() + ". Barrier: ";
    AddOutputDebug(s);
    s = barrier->display(4);
    AddOutputDebug(s);
    OUTPUT_DEBUG_END

    // Barrier is already updated from previous steps.
    // Get ref best feasible and infeasible, and then update
    // reference values.
    auto refBestFeas = barrier->getRefBestFeas();
    auto refBestInf  = barrier->getRefBestInf();

    barrier->updateRefBests();

    NOMAD::EvalPointPtr newBestFeas = barrier->getFirstXFeas();
    NOMAD::EvalPointPtr newBestInf  = barrier->getFirstXInf();

    if (nullptr != refBestFeas || nullptr != refBestInf)
    {
        // Compute success
        // Get which of newBestFeas and newBestInf is improving
        // the solution. Check newBestFeas first.
        NOMAD::ComputeSuccessType computeSuccess(evalType);
        std::shared_ptr<NOMAD::EvalPoint> newBest;
        NOMAD::SuccessType success = computeSuccess(newBestFeas, refBestFeas);
        if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            // newBestFeas is the improving point.
            newBest = newBestFeas;
            OUTPUT_DEBUG_START
            // Output Warning: When using '\n', the computed indentation for the
            // Step will be ignored. Leaving it like this for now. Using an
            // OutputInfo with AddMsg() would resolve the output layout.
            s = "Update: improving feasible point";
            if (refBestFeas)
            {
                s += " from\n    " + refBestFeas->display() + "\n";
            }
            s += " to " + newBestFeas->display();
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
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
                OUTPUT_DEBUG_START
                s = "Update: improving infeasible point";
                if (refBestInf)
                {
                    s+= " from\n    " + refBestInf->display() + "\n";
                }
                s += " to " + newBestInf->display();
                AddOutputDebug(s);
                OUTPUT_DEBUG_END
            }
        }
        if (success == NOMAD::SuccessType::UNSUCCESSFUL)
        {
            OUTPUT_DEBUG_START
            s = "Update: no success found";
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }


        // Debug verification
        // Compare computed success with value from MegaIteration.
        // This is the value from the previous MegaIteration. If it
        // was not evaluated, ignore the test.
        // If queue is not cleared between runs, also ignore the test.
        const bool clearEvalQueue = evc->getEvaluatorControlGlobalParams()->getAttributeValue<bool>("CLEAR_EVAL_QUEUE");
        const bool megaIterEvaluated = (NOMAD::SuccessType::NOT_EVALUATED != megaIter->getSuccessType());
        if (!clearEvalQueue && megaIterEvaluated && (success != megaIter->getSuccessType()))
        {
            s = "Warning: MegaIteration success type: ";
            s += NOMAD::enumStr(megaIter->getSuccessType());
            s += ". Is different than computed success type: " + NOMAD::enumStr(success);
            if (refBestFeas)
            {
                s += "\nRef best feasible:   " + refBestFeas->displayAll();
            }
            if (newBestFeas)
            {
                s += "\nNew best feasible:   " + newBestFeas->displayAll();
            }
            if (refBestInf)
            {
                s += "\nRef best infeasible: " + refBestInf->displayAll();
            }
            if (newBestInf)
            {
                s += "\nNew best infeasible: " + newBestInf->displayAll();
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
                s = "Error: Update cannot compute new direction for successful point: pointFromPtr is NULL ";
                s += newBest->display();
                throw NOMAD::StepException(__FILE__,__LINE__, s, this);
            }
            // PointFrom is in full dimension. Convert it to subproblem
            // to compute direction.
            auto fixedVariable = NOMAD::SubproblemManager::getSubFixedVariable(this);
            auto pointFromSub = std::make_shared<NOMAD::Point>(pointFromPtr->makeSubSpacePointFromFixed(fixedVariable));

            NOMAD::Direction dir = NOMAD::Point::vectorize(*pointFromSub, *pointNewPtr);
            OUTPUT_INFO_START
            std::string dirStr = "New direction " + dir.display();
            AddOutputInfo(dirStr);
            if (success == NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                AddOutputInfo("Last Iteration Improving. Delta remains the same.");
            }
            else
            {
                AddOutputInfo("Last Iteration Successful.");
            }
            OUTPUT_INFO_END

            // This computed direction may be used to sort points. Update values.
            // Use full space.
            auto pointNewFull = std::make_shared<NOMAD::Point>(pointNewPtr->makeFullSpacePointFromFixed(fixedVariable));
            auto dirFull = std::make_shared<NOMAD::Direction>(NOMAD::Point::vectorize(*pointFromPtr, *pointNewFull));
            NOMAD::EvcInterface::getEvaluatorControl()->setLastSuccessfulDir(dirFull);

            if (success >= NOMAD::SuccessType::FULL_SUCCESS)
            {
                // Update frame size for main mesh
                auto anisotropyFactor = _runParams->getAttributeValue<NOMAD::Double>("ANISOTROPY_FACTOR");
                bool anistropicMesh = _runParams->getAttributeValue<bool>("ANISOTROPIC_MESH");

                if (mesh->enlargeDeltaFrameSize(dir, anisotropyFactor, anistropicMesh))
                {
                    OUTPUT_INFO_START
                    AddOutputInfo("Delta is enlarged.");
                    OUTPUT_INFO_END
                }
                else
                {
                    OUTPUT_INFO_START
                    AddOutputInfo("Delta is not enlarged.");
                    OUTPUT_INFO_END
                }
            }
        }
        else
        {
            OUTPUT_INFO_START
            AddOutputInfo("Last Iteration Unsuccessful. Refine Delta.");
            OUTPUT_INFO_END
            mesh->refineDeltaFrameSize();
        }
    }

    mesh->updatedeltaMeshSize();
    OUTPUT_INFO_START
    AddOutputInfo("delta mesh  size = " + mesh->getdeltaMeshSize().display());
    AddOutputInfo("Delta frame size = " + mesh->getDeltaFrameSize().display());
    OUTPUT_INFO_END


    return true;
}


