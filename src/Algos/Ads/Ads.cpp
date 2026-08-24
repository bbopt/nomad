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


#include "../../Algos/Ads/Ads.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/Ads/AdsMegaIteration.hpp"


void NOMAD::Ads::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_ADS);
    
}

bool NOMAD::Ads::runImp()
{
    size_t k = 0;   // Iteration number (incremented at start)
    
    NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::UNDEFINED;
    
    if (!_termination->terminate(k))
    {
        std::shared_ptr<NOMAD::MeshBase> mesh;
        std::shared_ptr<NOMAD::BarrierBase> barrier;
        
        if (nullptr != _refMegaIteration)
        {
            // Case hot restart
            k       = _refMegaIteration->getK();
            barrier = _refMegaIteration->getBarrier();
            
            // Downcast from MegaIteration to MadsMegaIteration
            mesh    = (std::dynamic_pointer_cast<NOMAD::MadsMegaIteration> (_refMegaIteration ))->getMesh();
            megaIterSuccess = _refMegaIteration->getSuccessType();
            _success = megaIterSuccess;
        }
        else
        {
            mesh = dynamic_cast<NOMAD::MadsInitialization*>(_initialization.get())->getMesh();
            barrier = _initialization->getBarrier();
        }
        
        // Mads member _refMegaIteration is used for hot restart (read and write),
        // as well as to keep values used in Mads::end(), and may be used for _termination.
        // Update it here.
        _refMegaIteration = std::make_shared<NOMAD::AdsMegaIteration>(this, k, barrier, mesh, megaIterSuccess);
        
        // Create a MegaIteration for looping: manage multiple iterations on different
        // meshes and with different frame centers at the same time.
        NOMAD::AdsMegaIteration megaIteration(this, k, barrier, mesh, megaIterSuccess);
        while (!_termination->terminate(k))
        {
            
            megaIteration.start();
            megaIteration.run();
            megaIteration.end();
            
            // Counter is incremented when calling mega iteration end()
            k       = megaIteration.getK();
            
            if (!_algoSuccessful && megaIteration.getSuccessType() >= NOMAD::SuccessType::FULL_SUCCESS)
            {
                _algoSuccessful = true;
            }
            
            if (getUserInterrupt())
            {
                hotRestartOnUserInterrupt();
            }
        }
    }
    
    _termination->start();
    _termination->run();
    _termination->end();
    
    return _algoSuccessful;
}
