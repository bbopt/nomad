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

#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/SSDMads/SSDMads.hpp"
#include "../../Algos/SSDMads/SSDMadsMegaIteration.hpp"


void NOMAD::SSDMads::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_SSD_MADS);
    verifyParentNotNull();

    // Instantiate Mads initialization class
    _initialization = std::make_unique<NOMAD::MadsInitialization>( this );
}


void NOMAD::SSDMads::readInformationForHotRestart()
{
}


bool NOMAD::SSDMads::runImp()
{
    size_t k = 0;   // Iteration number
    NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::NOT_EVALUATED;

    bool runOk = true;

    // Note: _initialization is run in Algorithm::startImp().

    if (!_termination->terminate(k))
    {
        std::shared_ptr<NOMAD::MeshBase> mesh = dynamic_cast<NOMAD::MadsInitialization*>(_initialization.get())->getMesh();


        auto barrier = _initialization->getBarrier();

        // Mads member _megaIteration is used for hot restart (read and write),
        // as well as to keep values used in Mads::end(), and may be used for _termination.
        // Update it here.
        _megaIteration = std::make_shared<NOMAD::SSDMadsMegaIteration>(this, k, barrier, mesh, megaIterSuccess);

        while (!_termination->terminate(k))
        {

            // Create a MegaIteration to manage the pollster worker and the regular workers.
            NOMAD::SSDMadsMegaIteration ssdMegaIteration(this, k, barrier, mesh, megaIterSuccess);
            ssdMegaIteration.start();
            ssdMegaIteration.run();
            ssdMegaIteration.end();

            // Remember these values to construct the next MegaIteration.
            k       = ssdMegaIteration.getNextK();
            barrier = ssdMegaIteration.getBarrier();
            mesh    = ssdMegaIteration.getMesh();
            megaIterSuccess = ssdMegaIteration.getSuccessType();

            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }
        }
    }

    else
    {
        runOk = false;
    }

    _termination->start();
    _termination->run();
    _termination->end();

    return runOk;
}
