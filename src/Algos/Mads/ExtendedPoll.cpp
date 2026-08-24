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


#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/CatMads/CatMads.hpp"
#include "../../Algos/CatMads/CatExtendedPollMethod.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/ExtendedPoll.hpp"
#include "../../Output/OutputQueue.hpp"

#ifdef TIME_STATS
#include "../../Util/Clock.hpp"

// Initialize static variables
double NOMAD::ExtendedPoll::_extendedPollTime=0.0;
double NOMAD::ExtendedPoll::_extendedPollEvalTime=0.0;
#endif // TIME_STATS

void NOMAD::ExtendedPoll::init()
{
    setStepType(NOMAD::StepType::EXTENDED_POLL);
    verifyParentNotNull();
    
    // Default extended poll is not enabled because it does nothing
    _isEnabled = false;
    
    auto *catAlgo = getParentOfType<NOMAD::CatAlgoUtils*>(NOMAD::Step::_parentStep);
    if (nullptr != catAlgo )
    {
        _extendedPollMethod = std::make_shared<NOMAD::CatExtendedPollMethod>(NOMAD::Step::_parentStep);
        _isEnabled = true;
    }
    
}


void NOMAD::ExtendedPoll::startImp()
{
   // Sanity check.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    // Reset the current trial stats before run is called and increment the number of calls
    // This is managed by iteration utils when using generateTrialPoint instead of the (start, run, end) sequence.
    _trialPointStats.resetCurrentStats();
    _trialPointStats.incrementNbCalls();
    
}


bool NOMAD::ExtendedPoll::runImp()
{
    bool extendedPollSuccessful = false;
    std::string s;
    
    // Sanity check. The runImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);
    
    if (!_isEnabled)
    {
        // Early out - Return false: No new success found.
        OUTPUT_DEBUG_START
        AddOutputDebug("Extended poll method not enabled.");
        OUTPUT_DEBUG_END
        return false;
    }
    
    // A local user stop is requested. Do not perform remaining search methods. Stop type reset is done at the end of iteration/megaiteration and algorithm.
    if (_stopReasons->testIf(NOMAD::IterStopType::USER_ITER_STOP) || _stopReasons->testIf(NOMAD::IterStopType::USER_ALGO_STOP) ||
        _stopReasons->testIf(NOMAD::EvalGlobalStopType::CUSTOM_GLOBAL_STOP)) 
    {
        return false;
    }
    
    
    
#ifdef TIME_STATS
    double extendedPollStartTime = NOMAD::Clock::getCPUTime();
    double extendedPollEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS
    
    _extendedPollMethod->start();
    _extendedPollMethod->run();
    _extendedPollMethod->end();
    
#ifdef TIME_STATS
    _extendedPollTime += NOMAD::Clock::getCPUTime() - extendedPollStartTime;
    _extendedPollEvalTime += NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime() - extendedPollEvalStartTime;
#endif // TIME_STATS
    
    // Manually set the SuccessType from _extendedPollMethod to ExtendedPoll
    setSuccessType(_extendedPollMethod->getSuccessType());

    // Search is successful only if full success type.
    extendedPollSuccessful = (_extendedPollMethod->getSuccessType() >= NOMAD::SuccessType::FULL_SUCCESS);
    if (extendedPollSuccessful)
    {
        // Do not go through the other search methods if a search is
        // successful.
        OUTPUT_INFO_START
        s = _extendedPollMethod->getName();
        s += " is successful. Stop reason: ";
        s += _stopReasons->getStopReasonAsString() ;
        
        AddOutputInfo(s);
        OUTPUT_INFO_END
    }
    
    
    return extendedPollSuccessful;
}


void NOMAD::ExtendedPoll::endImp()
{
    // Sanity check. The endImp function should be called only when trial points are generated and evaluated for extended poll method separately.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);
    
    if (!_isEnabled)
    {
        // Early out
        return;
    }
    
    // Update the trial stats of the parent (Mads)
    // This is directly managed by iteration utils when using generateTrialPoint instead of the (start, run, end) sequence done here.
    _trialPointStats.updateParentStats();


    // Need to reset the EvalStopReason if a sub optimization is used during Search and the max bb is reached for this sub optimization
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (evc->testIf(NOMAD::EvalMainThreadStopType::LAP_MAX_BB_EVAL_REACHED))
    {
        evc->setStopReason(NOMAD::getThreadNum(), NOMAD::EvalMainThreadStopType::STARTED);
    }

}

