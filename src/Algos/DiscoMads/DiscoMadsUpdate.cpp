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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/DiscoMads/DiscoMadsMegaIteration.hpp"
#include "../../Algos/DiscoMads/DiscoMadsUpdate.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::DiscoMadsUpdate::init()
{
    setStepType(NOMAD::StepType::UPDATE);
    verifyParentNotNull();

    auto megaIter = getParentOfType<NOMAD::DiscoMadsMegaIteration*>();
    if (nullptr == megaIter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Error: An instance of class DiscoMadsUpdate must have a DiscoMegaIteration among its ancestors");
    }
    
    _clearEvalQueue = true;
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (nullptr != evc)
    {
        _clearEvalQueue = evc->getEvaluatorControlGlobalParams()->getAttributeValue<bool>("EVAL_QUEUE_CLEAR");
    }

}




bool NOMAD::DiscoMadsUpdate::runImp()
{
    bool lastIterationRevealing = getParentOfType<NOMAD::DiscoMadsMegaIteration*>()->isRevealing();

    // Number corresponding to current iteration
    size_t numIter = getParentOfType<NOMAD::DiscoMadsMegaIteration*>()->getK();

    // If last iteration was revealing, special update 
    if (lastIterationRevealing){
        // Barrier is already updated from previous steps (IterationUtils::postProcessing).
    
        // megaIter barrier is already in subproblem.
        // So no need to convert refBestFeas and refBestInf
        // from full dimension to subproblem.
        auto megaIter = getParentOfType<NOMAD::MadsMegaIteration*>();
        auto barrier = megaIter->getBarrier();
        auto mesh = megaIter->getMesh();
        std::string s;  // for output

        //Step 1. Get the best feasible and infeasible reference points, and then update reference values.
        auto refBestFeas = barrier->getRefBestFeas();
        auto refBestInf  = barrier->getRefBestInf();
        barrier->updateRefBests();

        // Display of current barrier at the beginning of the iteration considering updated ref points
        OUTPUT_DEBUG_START
        s = "Running " + getName() + ". Barrier: ";
        AddOutputDebug(s);
        std::vector<std::string> vs = barrier->display(4);
        for (const auto & si : vs)
        {
            AddOutputDebug(si);
        }
        s = "Update: revealing iteration";
        AddOutputDebug(s);
        OUTPUT_DEBUG_END

        //Step 2. Update mesh size parameter and frame parameter
        s = "Last Iteration revealing (iteration " + std::to_string(numIter-1) +"). Delta remains the same.";
        AddOutputInfo(s);
        mesh->updatedeltaMeshSize();   // update by security 

        OUTPUT_INFO_START
        AddOutputInfo("delta mesh size = " + mesh->getdeltaMeshSize().display());
        AddOutputInfo("Delta frame size = " + mesh->getDeltaFrameSize().display());
        OUTPUT_INFO_END

        //Step 3. Reset revealing status of  megaIteration in  DiscoMadsMegaIteration::RunImp

    }
    // If las iteration was not revealing, usual update of Mads
    else{
        MadsUpdate::runImp();
    }

    return true;
}


