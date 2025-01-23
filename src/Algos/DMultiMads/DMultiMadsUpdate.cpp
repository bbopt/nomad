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
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMegaIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsUpdate.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/SubproblemManager.hpp"

void NOMAD::DMultiMadsUpdate::init()
{
    setStepType(NOMAD::StepType::UPDATE);
    verifyParentNotNull();
}


std::string NOMAD::DMultiMadsUpdate::getName() const
{
    return getAlgoName() + NOMAD::stepTypeToString(_stepType);
}

bool NOMAD::DMultiMadsUpdate::runImp()
{
    std::string s;
    // Select the current incumbent (See Algo 4)
    // Two cases depending on previous iteration success:
    // - Success (full and partial success): select a possibly new incumbent point.
    // - Failure: decrease the mesh and select the current incumbent point

    // megaIter is already in subproblem.
    // So no need to convert from full dimension to subproblem.
    auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>();
    auto iter = getParentOfType<NOMAD::DMultiMadsIteration*>();
    auto barrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(megaIter->getBarrier());

    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "DMultiMadsUpdate: No barrier available");
    }
    
    if (NOMAD::EvalType::BB != barrier->getEvalType())
    {
        s = "DMultiMadsUpdate: Only BB eval type is handled";
    }

    if (NOMAD::ComputeType::STANDARD != barrier->getComputeType())
    {
        s = "DMultiMadsUpdate: Only STANDARD compute type is handled";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    OUTPUT_DEBUG_START
    s = "Running " + getName();
    AddOutputDebug(s);
    OUTPUT_DEBUG_END

    // If the previous iteration is not a success, the meshes of current incumbents are refined (no success). If not, the current incumbents may change.
    NOMAD::SuccessType previousSuccess = iter->getPreviousSuccessType();
    
    // If previous success type is not a success (partial or full), reduce the mesh associated to the current frame center.
    // The frame center can be set if it does not already exist (just after initialization)
    if (previousSuccess < NOMAD::SuccessType::PARTIAL_SUCCESS)
    {
        OUTPUT_DEBUG_START
        s = "Update: previous iter, NO success found";
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
        
        // If the iteration has a frameCenter, it will remain unchanged, except for the mesh.
        if ( nullptr != iter->getFrameCenter())
        {
            iter->getFrameCenter()->getMesh()->refineDeltaFrameSize();
        }
    }
    else
    {
        // No need to update incumbents if mesh of incumbents has not changed.
        // Incumbents have been updated when updated the barrier with eval points
        
        OUTPUT_DEBUG_START
        s = "Update: previous iter, success found";
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
    }

    OUTPUT_DEBUG_START
    s = "Barrier: ";
    AddOutputDebug(s);
    std::vector<std::string> vs = barrier->display(100, true);
    for (const auto & si : vs)
    {
        AddOutputDebug(si);
    }
    OUTPUT_DEBUG_END

    // Before selecting new frame centers, update current incumbents.
    barrier->updateCurrentIncumbents();


    // Select frame center. Can change if iteration is successful or not.
    // Current incumbents may have changed because barrier is updated by adding points or by because the frame size of a barrier point has changed (see above).
    auto currentBestFeas =  barrier->getCurrentIncumbentFeas();
    auto currentBestInf = barrier->getCurrentIncumbentInf();
    
    if (nullptr == currentBestFeas && nullptr == currentBestInf)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Update cannot set iteration frame center (no current best feas or infeas).");
    }
    
    if ( nullptr != currentBestFeas )
    {
        iter->setFrameCenter(currentBestFeas);
    }
    else if ( nullptr != currentBestInf )
    {
        iter->setFrameCenter(currentBestInf);
    }
    OUTPUT_INFO_START
    AddOutputInfo("Frame center: " + iter->getFrameCenter()->display());
    AddOutputInfo("Number of points in Lk (feas+inf): " + std::to_string(barrier->nbXFeas()+barrier->nbXInf()));
    AddOutputInfo("delta mesh size = " + iter->getFrameCenter()->getMesh()->getdeltaMeshSize().display());
    AddOutputInfo("Delta frame size = " + iter->getFrameCenter()->getMesh()->getDeltaFrameSize().display());
    OUTPUT_INFO_END
    
    
    return true;
}
