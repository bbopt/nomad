/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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

#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/fileutils.hpp"

void NOMAD::Mads::init()
{
    _name = "MADS";

    // Instantiate Mads initialization class
    _initialization = std::make_unique<NOMAD::MadsInitialization>( this );

}


bool NOMAD::Mads::runImp()
{
    size_t k = 0;   // Iteration number
    NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::NOT_EVALUATED;

    bool successFound = false;

    if (!_termination->terminate(k))
    {
        std::shared_ptr<NOMAD::MeshBase> mesh;
        std::shared_ptr<NOMAD::Barrier> barrier;

        if (nullptr != _megaIteration)
        {
            // Case hot restart
            k       = _megaIteration->getK();
            barrier = _megaIteration->getBarrier();

            // Downcast from MegaIteration to MadsMegaIteration
            mesh    = (std::dynamic_pointer_cast<NOMAD::MadsMegaIteration> (_megaIteration ))->getMesh();
            megaIterSuccess = _megaIteration->getSuccessType();
        }
        else
        {
            mesh = dynamic_cast<NOMAD::MadsInitialization*>(_initialization.get())->getMesh();
            barrier = _initialization->getBarrier();
        }

        // Mads member _megaIteration is used for hot restart (read and write),
        // as well as to keep values used in Mads::end(), and may be used for _termination.
        // Update it here.
        _megaIteration = std::make_shared<NOMAD::MadsMegaIteration>(this, k, barrier, mesh, megaIterSuccess);


        while (!_termination->terminate(k))
        {
            // Create an MegaIteration: manage multiple iterations on different
            // meshes and with different frame centers at the same time.
            NOMAD::MadsMegaIteration megaIteration(this, k, barrier, mesh, megaIterSuccess);
            megaIteration.start();
            megaIteration.run();
            megaIteration.end();

            // Remember these values to construct the next MegaIteration.
            k       = megaIteration.getNextK();
            barrier = megaIteration.getBarrier();
            mesh    = megaIteration.getMesh();
            megaIterSuccess = megaIteration.getSuccessType();

            if (!successFound && megaIterSuccess >= NOMAD::SuccessType::FULL_SUCCESS)
            {
                successFound = true;
            }

            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }
        }
    }

    _termination->start();
    _termination->run();
    _termination->end();

    return successFound;
}


void NOMAD::Mads::hotRestartOnUserInterrupt()
{
#ifdef TIME_STATS
    if (isRootAlgo())
    {
        _totalCPUAlgoTime += NOMAD::Clock::getCPUTime() - _startTime;
    }
#endif // TIME_STATS
    hotRestartBeginHelper();

    // Reset mesh because parameters have changed.
    std::stringstream ss;
    auto mesh = getIterationMesh();
    if (nullptr != mesh)
    {
        ss << *mesh;
        // Reset pointer
        mesh.reset();
        mesh = std::make_shared<NOMAD::GMesh>(_pbParams);
        // Get old mesh values
        ss >> *mesh;
    }

    hotRestartEndHelper();
#ifdef TIME_STATS
    if (isRootAlgo())
    {
        _startTime = NOMAD::Clock::getCPUTime();
    }
#endif // TIME_STATS
}


void NOMAD::Mads::readInformationForHotRestart()
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
            std::string s = "Read hot restart file " + hotRestartFile;
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_NORMAL);

            // Create a GMesh and an MadsMegaIteration with default values, to be filled
            // by istream is.
            // NOTE: Working in full dimension
            auto barrier = std::make_shared<NOMAD::Barrier>(NOMAD::INF, NOMAD::Point(_pbParams->getAttributeValue<size_t>("DIMENSION")), NOMAD::EvalType::BB);
            std::shared_ptr<NOMAD::MeshBase> mesh = std::make_shared<NOMAD::GMesh>(_pbParams);

            _megaIteration = std::make_shared<NOMAD::MadsMegaIteration>(this, 0, barrier, mesh, NOMAD::SuccessType::NOT_EVALUATED);

            // Here we use Algorithm::operator>>
            NOMAD::read<NOMAD::Mads>(*this, hotRestartFile);
        }
    }
}
