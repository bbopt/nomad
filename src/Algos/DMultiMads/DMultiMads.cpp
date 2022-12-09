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

// Algo specifics
#include "../../Algos/DMultiMads/DMultiMads.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMegaIteration.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"

void NOMAD::DMultiMads::init()
{

    setStepType(NOMAD::StepType::ALGORITHM_DMULTIMADS);

    // Instantiate algorithm Initialization class (Start function automatically called)
   // The Mads initialization manages Mesh and X0
    _initialization = std::make_unique<NOMAD::MadsInitialization>( this, true /*use Cache for barrier init*/, true /*initialization for DMultiMads*/ );
    
    if (NOMAD::Algorithm::getNbObj() < 2)
    {
        throw NOMAD::InvalidParameter(__FILE__,__LINE__,"DMultiMads is intended to solve problems with more than one objective.");
    }
    
}

bool NOMAD::DMultiMads::runImp()
{
    _algoSuccessful = false;
    
    if ( !_runParams->getAttributeValue<bool>("DMULTIMADS_OPTIMIZATION") )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMads is a standalone optimization algo. Cannot be used as a Mads search method.");
    }

    if ( ! _stopReasons->checkTerminate() )
    {
        size_t k = 0;   // Iteration number

        // DMultiMadsBarrier created during Initialization (with X0).
        std::shared_ptr<NOMAD::BarrierBase> barrier = _initialization->getBarrier();
        
        // Mesh created during Initialization
        NOMAD::MeshBasePtr initialMesh = dynamic_cast<NOMAD::MadsInitialization*>(_initialization.get())->getMesh();

        // Create a single MegaIteration: manage multiple iterations.
        NOMAD::DMultiMadsMegaIteration megaIteration(this, k, barrier, initialMesh, NOMAD::SuccessType::UNDEFINED);
        while (!_termination->terminate(k))
        {
            megaIteration.start();
            megaIteration.run();
            megaIteration.end();

            k       = megaIteration.getK();
            
            if (!_algoSuccessful && megaIteration.getSuccessType() >= NOMAD::SuccessType::FULL_SUCCESS)
            {
                _algoSuccessful = true;
            }
            
//            if (_userInterrupt)
//            {
//                hotRestartOnUserInterrupt();
//            }
        }

//        // _refMegaIteration is used for hot restart (read
//        // and write), as well as to keep values used in Mads::end(). Update it here.
        _refMegaIteration = std::make_shared<NOMAD::DMultiMadsMegaIteration>(this, k, barrier, nullptr, _success);

        _termination->start();
        _termination->run();
        _termination->end();
    }

    return _algoSuccessful;
}


void NOMAD::DMultiMads::readInformationForHotRestart()
{
    // Not implemented
}
