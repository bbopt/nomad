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
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Util/fileutils.hpp"

// QPSolver algo specifics
#include "../../Algos/QPSolverAlgo/QPSolverAlgo.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverAlgoMegaIteration.hpp"

// QuadModel specifics
#include "../../Algos/QuadModel/QuadModelInitialization.hpp"


void NOMAD::QPSolverAlgo::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_QPSOLVER);

    bool qpsolverAlgoOpt = _runParams->getAttributeValue<bool>("QP_OPTIMIZATION"); // true if standalone
    
    if (!qpsolverAlgoOpt)
    {
        throw NOMAD::InvalidParameter(__FILE__,__LINE__,"QP algo is intended for standalone optimization. Set QP_OPTIMIZATION true.");
    }
    _initialization = std::make_unique<NOMAD::QuadModelInitialization>(this);

}

bool NOMAD::QPSolverAlgo::runImp()
{
    _algoSuccessful = false;
    

    if ( ! _stopReasons->checkTerminate() )
    {
        size_t k = 0;   // Iteration number
        
        // Barrier was computed by Initialization (if X0 provided).
        // Barrier is used for MegaIteration management.
        std::shared_ptr<NOMAD::BarrierBase> barrier = _initialization->getBarrier();
        if (nullptr == barrier)
        {
            // Barrier constructor automatically finds the best points in the cache.
            
            auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
            auto hNormType = _runParams->getAttributeValue<NOMAD::HNormType>("H_NORM");
            
            // Compute type for this optim
            FHComputeTypeS computeType; // Default from struct initializer
            computeType.hNormType = hNormType; // REM: No PhaseOne search for this algo!
        
            // Eval type for this optim
            auto evalType = NOMAD::EvcInterface::getEvaluatorControl()->getCurrentEvalType();
            
            // Create a single objective progressive barrier
            barrier = std::make_shared<NOMAD::ProgressiveBarrier>(hMax,
                                                       NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                                       evalType,
                                                       computeType);
        }

        // Create a single MegaIteration: manage multiple iterations.
        NOMAD::QPSolverAlgoMegaIteration megaIteration(this, k, barrier,NOMAD::SuccessType::UNDEFINED);
        while (!_termination->terminate(k))
        {
            megaIteration.start();
            bool currentMegaIterSuccess = megaIteration.run();
            megaIteration.end();

            _algoSuccessful = _algoSuccessful || currentMegaIterSuccess;

            k       = megaIteration.getK();
            // NOMAD::SuccessType megaIterSuccess = megaIteration.getSuccessType();
            
            if (getUserInterrupt())
            {
                hotRestartOnUserInterrupt();
            }
        }

        // _refMegaIteration is used for hot restart (read
        // and write), as well as to keep values used in Mads::end(). Update it here.
        _refMegaIteration = std::make_shared<NOMAD::QPSolverAlgoMegaIteration>(this, k, barrier, _success);

        _termination->start();
        _termination->run();
        _termination->end();
    }

    return _algoSuccessful;
}

void NOMAD::QPSolverAlgo::readInformationForHotRestart()
{
    // Restart from where we were before.
    // For this, we need to read some files.
    // Note: Cache file is treated independently from hot restart file.

    if (_runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES"))
    {
        // Verify the files exist and are readable.
        const std::string& hotRestartFile = _runParams->getAttributeValue<std::string>("HOT_RESTART_FILE");
        if (NOMAD::checkReadFile(hotRestartFile))
        {
            std::cout << "Read hot restart file " << hotRestartFile << std::endl;

            auto barrier = _initialization->getBarrier();
            int k = 0;
            NOMAD::SuccessType success = NOMAD::SuccessType::UNDEFINED;

            _refMegaIteration = std::make_shared<NOMAD::QPSolverAlgoMegaIteration>(this, k, barrier, success);

            // Here we use QPSolverAlgo::operator>>
            NOMAD::read<QPSolverAlgo>(*this, hotRestartFile);
        }
    }
}
