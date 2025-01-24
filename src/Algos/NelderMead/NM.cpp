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
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Util/fileutils.hpp"

// NM specific
#include "../../Algos/NelderMead/NM.hpp"
#include "../../Algos/NelderMead/NMInitialization.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"

// DMultiMads Specific
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"

void NOMAD::NM::init()
{

    setStepType(NOMAD::StepType::ALGORITHM_NM);

    // Instantiate NM initialization class
    _initialization = std::make_unique<NOMAD::NMInitialization>( this );
}

bool NOMAD::NM::runImp()
{
    _algoSuccessful = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        size_t k = 0;   // Iteration number

        std::shared_ptr<NOMAD::BarrierBase> barrier = nullptr;
        
        auto nmOpt = _runParams->getAttributeValue<bool>("NM_OPTIMIZATION");

        if (nmOpt)
        {
            // Barrier was computed by Initialization.
            barrier = _initialization->getBarrier();
        }
        else
        {
            // Get barrier from upper MegaIteration, if available.
            auto megaIter = getParentOfType<NOMAD::MegaIteration*>(false);
            if (nullptr != megaIter)
            {
                barrier = megaIter->getBarrier();
            }
            
        }
        
        NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::UNDEFINED;
        
        // Create a MegaIteration: manage multiple iterations.
        NOMAD::NMMegaIteration megaIteration(this, k, barrier, megaIterSuccess);
        while (!_termination->terminate(k))
        {
            megaIteration.start();
            bool currentMegaIterSuccess = megaIteration.run();
            megaIteration.end();

            _algoSuccessful = _algoSuccessful || currentMegaIterSuccess;

            // For passing to _refMegaIteration
            k       = megaIteration.getK();
            megaIterSuccess = megaIteration.getSuccessType();

            if (getUserInterrupt())
            {
                hotRestartOnUserInterrupt();
            }
        }
        
        // Issue #372: For hot restart make sure to save the simplex (maybe as X0s)
        // _refMegaIteration is used for hot restart (read
        // and write), as well as to keep values used in Mads::end(). Update it here.
        _refMegaIteration = std::make_shared<NOMAD::NMMegaIteration>(this, k, barrier, megaIterSuccess);

        _termination->start();
        _termination->run();
        _termination->end();
    }

    return _algoSuccessful;
}


void NOMAD::NM::readInformationForHotRestart()
{
    // Restart from where we were before.
    // For this, we need to read some files.
    // Note: Cache file is treated independently of hot restart file.

    if (_runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES"))
    {
        // Verify the files exist and are readable.
        const std::string& hotRestartFile = _runParams->getAttributeValue<std::string>("HOT_RESTART_FILE");
        if (NOMAD::checkReadFile(hotRestartFile))
        {
            std::cout << "Read hot restart file " << hotRestartFile << std::endl;

            // Create a GMesh and a MegaIteration with default values, to be filled
            // by istream is.
            // Issue #372: Fix potential bug with Hot Restart
            // Note: Assuming the progressive barrier read is in the same subspace as the current subspace.
            // This could be fixed if we write and read the progressive barrier in full subspace.
            
            // Create a single objective barrier
            auto barrier = std::make_shared<NOMAD::ProgressiveBarrier>();
            int k = 0;
            NOMAD::SuccessType success = NOMAD::SuccessType::UNDEFINED;


            _refMegaIteration = std::make_shared<NOMAD::NMMegaIteration>(this, k, barrier, success);

            // Here we use NM::operator>>
            NOMAD::read<NM>(*this, hotRestartFile);
        }
    }
}
