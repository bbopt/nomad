/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
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

#include "../Algos/Algorithm.hpp"
#include "../Algos/EvcInterface.hpp"
#include "../Algos/Termination.hpp"
#include "../Util/Clock.hpp"

void NOMAD::Termination::init()
{
    setStepType(NOMAD::StepType::TERMINATION);
    verifyParentNotNull();
}


void NOMAD::Termination::startImp()
{
}


bool NOMAD::Termination::runImp()
{
    return _stopReasons->checkTerminate() ;
}


bool NOMAD::Termination::terminate(size_t iteration)
{
    bool stop = _stopReasons->checkTerminate();
    if (stop)
    {
        // A stop condition was already reached.
        return stop;
    }

    // Set stopReason due to criterions other than AlgoStopReasons<>
    auto maxIterations = _runParams->getAttributeValue<size_t>("MAX_ITERATIONS");
    auto maxTime = _runParams->getAttributeValue<size_t>("MAX_TIME");

    // Termination conditions go here.
    // This is also tested in EvaluatorControl
    if (NOMAD::Step::getUserTerminate())
    {
        // Force quit (by pressing CTRL-C):
        _stopReasons->set(NOMAD::BaseStopType::CTRL_C);
    }
    else if (maxIterations < NOMAD::INF_SIZE_T && iteration > maxIterations)
    {
        // Max iterations reached
        _stopReasons->set(NOMAD::IterStopType::MAX_ITER_REACHED);
    }
    else if (maxTime < NOMAD::INF_SIZE_T && NOMAD::Clock::getTimeSinceStart() >= maxTime)
    {
        // Max time reached
        _stopReasons->set(NOMAD::BaseStopType::MAX_TIME_REACHED);
    }
    else if (_runParams->getAttributeValue<bool>("STOP_IF_FEASIBLE") && solHasFeas())
    {
        _stopReasons->set(NOMAD::IterStopType::STOP_ON_FEAS);
    }
    else if (_runParams->getAttributeValue<bool>("STOP_IF_PHASE_ONE_SOLUTION") && hasPhaseOneSolution())
    {
        _stopReasons->set(NOMAD::IterStopType::PHASE_ONE_COMPLETED);
    }
    else
    {
        // Need to check on MaxEval and MaxBBEval a last time because in evaluatorControl
        // the stop reason may have been set due to all queue points evaluated.
        stop = NOMAD::EvcInterface::getEvaluatorControl()->reachedMaxEval();
    }

    stop = stop || _stopReasons->checkTerminate() ;
    return stop;
}


void NOMAD::Termination::endImp()
{
    const NOMAD::Algorithm* currentAlgo = getParentOfType<NOMAD::Algorithm*>();
    NOMAD::OutputLevel outputLevel = currentAlgo->isSubAlgo() ? NOMAD::OutputLevel::LEVEL_INFO
                                                              : NOMAD::OutputLevel::LEVEL_HIGH;

    if (_stopReasons->checkTerminate())
    {
        std::string terminationInfo = "A termination criterion is reached: ";
        terminationInfo += _stopReasons->getStopReasonAsString();
        auto evc = NOMAD::EvcInterface::getEvaluatorControl();
        if (_stopReasons->testIf(NOMAD::EvalGlobalStopType::MAX_BB_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getBbEval());
        }
        else if (_stopReasons->testIf(NOMAD::EvalGlobalStopType::MAX_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getNbEval());
        }
        else if (_stopReasons->testIf(NOMAD::EvalGlobalStopType::MAX_BLOCK_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getBlockEval());
        }
        else if (evc->testIf(NOMAD::EvalMainThreadStopType::MAX_MODEL_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getTotalModelEval());
        }
        else if (evc->testIf(NOMAD::EvalMainThreadStopType::LAP_MAX_BB_EVAL_REACHED))
        {
            terminationInfo += " " + NOMAD::itos(evc->getLapBbEval());
        }
        AddOutputInfo(terminationInfo, outputLevel);
    }
    else
    {
        std::string terminationInfo = "No termination criterion reached";
        AddOutputInfo(terminationInfo, outputLevel);
    }

}
