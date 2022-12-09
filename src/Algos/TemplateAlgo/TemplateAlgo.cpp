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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Util/fileutils.hpp"

// Template algo specifics
#include "../../Algos/TemplateAlgo/TemplateAlgo.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoInitialization.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoMegaIteration.hpp"

void NOMAD::TemplateAlgo::init()
{

    setStepType(NOMAD::StepType::ALGORITHM_RANDOM);

    // Instantiate random algorithm Initialization class (called automatically by start)
    _initialization = std::make_unique<NOMAD::TemplateAlgoInitialization>( this );
}

bool NOMAD::TemplateAlgo::runImp()
{
    _algoSuccessful = false;
    
    bool randomAlgoOpt = _runParams->getAttributeValue<bool>("RANDOM_ALGO_OPTIMIZATION");

    if ( ! _stopReasons->checkTerminate() )
    {
        size_t k = 0;   // Iteration number

        std::shared_ptr<NOMAD::BarrierBase> barrier = nullptr;

        // Handle two situations for calling this: standalone optimization or search method
        if (randomAlgoOpt)
        {
            // Barrier was computed by Initialization (if X0 provided).
            barrier = _initialization->getBarrier();            
        }
        else
        {
            // Get barrier from upper MadsMegaIteration, if available (it is the case when Mads calls this as part of a search method).
            auto madsMegaIter = getParentOfType<NOMAD::MadsMegaIteration*>(false);
            if (nullptr != madsMegaIter)
            {
                barrier = madsMegaIter->getBarrier();
            }
        }

        

        // Create a single MegaIteration: manage multiple iterations.
        NOMAD::TemplateAlgoMegaIteration megaIteration(this, k, barrier,NOMAD::SuccessType::UNDEFINED);
        while (!_termination->terminate(k))
        {
            megaIteration.start();
            bool currentMegaIterSuccess = megaIteration.run();
            megaIteration.end();

            _algoSuccessful = _algoSuccessful || currentMegaIterSuccess;

            k       = megaIteration.getK();
            NOMAD::SuccessType megaIterSuccess = megaIteration.getSuccessType();

            if (!randomAlgoOpt && megaIterSuccess !=NOMAD::SuccessType::FULL_SUCCESS) // Search method stops if not full success
            {
                auto algoStopReason = NOMAD::AlgoStopReasons<NOMAD::RandomAlgoStopType>::get ( _stopReasons );
                algoStopReason->set( NOMAD::RandomAlgoStopType::SINGLE_PASS_COMPLETED ); // This will stop iterations.
            }
            
            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }
        }

        // _refMegaIteration is used for hot restart (read
        // and write), as well as to keep values used in Mads::end(). Update it here.
        _refMegaIteration = std::make_shared<NOMAD::TemplateAlgoMegaIteration>(this, k, barrier, _success);

        _termination->start();
        _termination->run();
        _termination->end();
    }

    return _algoSuccessful;
}


void NOMAD::TemplateAlgo::readInformationForHotRestart()
{
    // Restart from where we were before.
    // For this, we need to read some files.
    // Note: Cache file is treated independently from hot restart file.

    if (_runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES"))
    {
        // Verify the files exist and are readable.
        std::string hotRestartFile = _runParams->getAttributeValue<std::string>("HOT_RESTART_FILE");
        if (NOMAD::checkReadFile(hotRestartFile))
        {
            std::cout << "Read hot restart file " << hotRestartFile << std::endl;
            
            auto barrier = std::make_shared<NOMAD::Barrier>();
            int k = 0;
            NOMAD::SuccessType success = NOMAD::SuccessType::UNDEFINED;

            _refMegaIteration = std::make_shared<NOMAD::TemplateAlgoMegaIteration>(this, k, barrier, success);

            // Here we use TemplateAlgo::operator>>
            NOMAD::read<TemplateAlgo>(*this, hotRestartFile);
        }
    }
}
