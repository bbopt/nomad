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

#include <algorithm>

#include "../../Algos/SimpleMads/SimpleMads.hpp"
#include "../../Algos/SubproblemManager.hpp"


void NOMAD::SimpleMads::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_MADS);
    
    // We can accept Mads with more than one objective when doing a PhaseOneSearch of DMultiMads optimization.
    if (!_runParams->getAttributeValue<bool>("DMULTIMADS_OPTIMIZATION") && NOMAD::Algorithm::getNbObj() > 1)
    {
        throw NOMAD::InvalidParameter(__FILE__,__LINE__,"Mads solves single objective problems. To handle several objectives please use DMultiMads: DMULTIMADS_OPTIMIZATION yes");
    }
    
    
    if (nullptr == _poll.getBarrier())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"A simple mads must have a poll barrier.");
    }
    

}

bool NOMAD::SimpleMads::runImp()
{
    size_t k = 0;   // Iteration number (incremented at start)
    
    bool pollSuccess = true;
    
    // Termination: 1- Reach max model eval, 2- fail to create trial points on mesh
    while ( _poll.getNbEval() < _maxEval && pollSuccess)
    {
        _poll.start();
        pollSuccess = _poll.run();
        _poll.end();
        k++;
    }

    return true;
}

const NOMAD::SimpleEvalPoint & NOMAD::SimpleMads::getBestSimpleSolution(bool bestFeas) const
{
    
    const auto & barrier = _poll.getBarrier();
    bool runInPhaseOneSearch = _poll.getPhaseOneSearch();

    if (bestFeas && !runInPhaseOneSearch)
    {
        // If poll is in PhaseOneSearch the incumbentFeas is
        // PhaseOne feasible but not Standard feasible
        return barrier->getCurrentIncumbentFeas();
    }
    else
    {
        // If poll is in PhaseOneSearch the incumbentFeas is
        // Standard infeasible
        if (!runInPhaseOneSearch)
        {
            return barrier->getCurrentIncumbentInf();
        }
        else
        {
            return barrier->getCurrentIncumbentFeas();
        }
    }
    
}


NOMAD::EvalPoint NOMAD::SimpleMads::getBestSolution(bool bestFeas) const
{
    NOMAD::EvalPoint bestSol;

    const auto & barrier = _poll.getBarrier();
    bool runInPhaseOneSearch = _poll.getPhaseOneSearch();
    if (nullptr != barrier)
    {
        if (bestFeas)
        {
            // If poll is in PhaseOneSearch the incumbentFeas is
            // PhaseOne feasible but not Standard feasible
            if (!runInPhaseOneSearch)
            {
                bestSol = static_cast<NOMAD::EvalPoint>(barrier->getCurrentIncumbentFeas());
            }
        }
        else
        {
            // If poll is in PhaseOneSearch the incumbentFeas is
            // Standard infeasible
            if (!runInPhaseOneSearch)
            {
                bestSol = static_cast<NOMAD::EvalPoint>(barrier->getCurrentIncumbentInf());
            }
            else
            {
                bestSol = static_cast<NOMAD::EvalPoint>(barrier->getCurrentIncumbentFeas());
            }
        }
        
        if (bestSol.isComplete())
        {
            auto fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this);
            bestSol = bestSol.makeFullSpacePointFromFixed(fixedVariable);
        }
    }
    return bestSol;
}


void NOMAD::SimpleMads::endImp()
{
    
    if ( _endDisplay )
    {
        endDisplay();
    }
}


void NOMAD::SimpleMads::endDisplay() const
{
    // Display the best feasible solutions.
    std::string sFeas;
    
    // Output level is info if this algorithm is a sub part of another algorithm.
    NOMAD::OutputLevel outputLevel = NOMAD::OutputLevel::LEVEL_INFO;
    
    if ( ! NOMAD::OutputQueue::GoodLevel(outputLevel) )
    {
        return;
    }
    
    NOMAD::FHComputeTypeS computeType; /*default initializer is used*/
    
    auto solFormat = NOMAD::OutputQueue::getInstance()->getSolFormat();
    
    NOMAD::OutputInfo displaySolFeas(getName(), sFeas, outputLevel);
    
    sFeas = "Best feasible solution";
    
    auto bestFeas = getBestSolution(true);
    
    if (!bestFeas.isComplete())
    {
        sFeas += ":     Undefined.";
        displaySolFeas.addMsg(sFeas);
    }
    else
    {
        sFeas += ":     ";
        displaySolFeas.addMsg(sFeas + bestFeas.display(computeType,
                                                       solFormat,
                                                       NOMAD::DISPLAY_PRECISION_FULL,
                                                       false));
    }
    
    NOMAD::OutputQueue::Add(std::move(displaySolFeas));
    
    // Display the best infeasible solutions.
    std::string sInf;
    NOMAD::OutputInfo displaySolInf(getName(), sInf, outputLevel);
    sInf = "Best infeasible solution";
    
    auto bestInf = getBestSolution(false);
    
    if ( !bestInf.isComplete())
    {
        sInf += ":   Undefined.";
        displaySolInf.addMsg(sInf);
    }
    else
    {
        sInf += ":   ";
        displaySolInf.addMsg(sInf + bestInf.display(computeType,
                                                    solFormat,
                                                    NOMAD::DISPLAY_PRECISION_FULL,
                                                    false));
    }
    NOMAD::OutputQueue::Add(std::move(displaySolInf));
        
    std::string sNbEval = "Function evaluations: " + NOMAD::itos(_poll.getNbEval());
    NOMAD::OutputQueue::Add(sNbEval);
}
