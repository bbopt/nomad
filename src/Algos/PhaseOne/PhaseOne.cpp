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
#include "../../Algos/PhaseOne/PhaseOne.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Output/OutputDirectToFile.hpp"

void NOMAD::PhaseOne::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_PHASE_ONE);
    verifyParentNotNull();

}


void NOMAD::PhaseOne::startImp()
{
    // Temporarily disable solution file (restored in endImp())
    NOMAD::OutputDirectToFile::getInstance()->disableSolutionFile();

    // Default algorithm start
    NOMAD::Algorithm::startImp();
    
    // Setup the run parameters to stop once a point that satisfies EB constraints is obtained
    _runParams = std::make_shared<NOMAD::RunParameters>(*_runParams);
    _runParams->setAttributeValue("STOP_IF_PHASE_ONE_SOLUTION", true);
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();
    _runParams->checkAndComply(evcParams, _pbParams);

    // Setup Mads
    _madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
    _mads = std::make_shared<NOMAD::Mads>(this, _madsStopReasons, _runParams, _pbParams, false /*false: barrier is not initialized from cache*/);
}


bool NOMAD::PhaseOne::runImp()
{
    bool madsSuccess = false;

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    auto previousComputeType = evc->getComputeType();
    evc->setComputeType(NOMAD::ComputeType::PHASE_ONE);

    // Run Mads on Phase One.
    _mads->start();
    madsSuccess = _mads->run();
    _mads->end();

    evc->setComputeType(previousComputeType);
    evc->resetBestIncumbent(-1); // Reset for display (-1 for all main threads)
    // Update PhaseOne stop reasons
    if ( !_mads->hasPhaseOneSolution())
    {
        auto phaseOneStopReasons = NOMAD::AlgoStopReasons<NOMAD::PhaseOneStopType>::get(_stopReasons);
        if (!madsSuccess)
        {
            phaseOneStopReasons->set(NOMAD::PhaseOneStopType::MADS_FAIL);
        }
        else
        {
            phaseOneStopReasons->set(NOMAD::PhaseOneStopType::NO_FEAS_PT);
        }
        
        // Phase one has failed to get a solution with EB constraints feasible.
        return false;
    }
    
    return true;
}


void NOMAD::PhaseOne::endImp()
{
    NOMAD::Algorithm::endImp();
    
    // Ensure evaluation of queue will continue
    NOMAD::EvcInterface::getEvaluatorControl()->restart();
    NOMAD::EvcInterface::getEvaluatorControl()->setLastSuccessfulFeasDir(nullptr);
    NOMAD::EvcInterface::getEvaluatorControl()->setLastSuccessfulInfDir(nullptr);

    // Re-enable writing in Solution file
    NOMAD::OutputDirectToFile::getInstance()->enableSolutionFile();

    if (solHasFeas())
    {
        std::vector<NOMAD::EvalPoint> evalPointList;
        NOMAD::CacheInterface cacheInterface(this);
        size_t numFeas = cacheInterface.findBestFeas(evalPointList, NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD);
        if (numFeas > 0)
        {
            // Evaluation info for output
            NOMAD::StatsInfo info;

            info.setBBO(evalPointList[0].getBBO(NOMAD::EvalType::BB));
            info.setSol(*(evalPointList[0].getX()));

            NOMAD::OutputDirectToFile::Write(info,true,false); // Write in solution (if solution_file exists) but not in history file
        }
    }
}
