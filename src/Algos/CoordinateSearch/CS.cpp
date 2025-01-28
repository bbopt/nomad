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
#include "../../Algos/CoordinateSearch/CSMesh.hpp"
#include "../../Algos/CoordinateSearch/CS.hpp"
#include "../../Algos/CoordinateSearch/CSInitialization.hpp"
#include "../../Algos/CoordinateSearch/CSIteration.hpp"
#include "../../Algos/CoordinateSearch/CSUpdate.hpp"
#include "../../Algos/CoordinateSearch/CSMegaIteration.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"


#include "../../Algos/EvcInterface.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Util/fileutils.hpp"
#ifdef TIME_STATS
#include "../../Util/Clock.hpp"
#endif

void NOMAD::CS::init(bool barrierInitializedFromCache)
{
    setStepType(NOMAD::StepType::ALGORITHM_CS);
    
    // Instantiate Mads initialization class
    _initialization = std::make_unique<NOMAD::CSInitialization>( this , barrierInitializedFromCache);

}

bool NOMAD::CS::runImp()
{
    size_t k = 1;   // Iteration number
    
    bool successFound = false;
    
    NOMAD::SuccessType megaIterSuccess;
    
    if (!_termination->terminate(k))
    {
        std::shared_ptr<NOMAD::MeshBase> mesh;
        std::shared_ptr<NOMAD::BarrierBase> barrier;
        
        if (nullptr != _refMegaIteration)
        {
            // Case hot restart
            k       = _refMegaIteration->getK();
            barrier = _refMegaIteration->getBarrier();
            
            // Downcast from MegaIteration to CSMegaIteration
            mesh    = (std::dynamic_pointer_cast<NOMAD::CSMegaIteration> (_refMegaIteration ))->getMesh();
            megaIterSuccess = _refMegaIteration->getSuccessType();
        }
        else
        {
            mesh = dynamic_cast<NOMAD::CSInitialization*>(_initialization.get())->getMesh();
            barrier = _initialization->getBarrier();
        }
        
        // Mads member _megaIteration is used for hot restart (read and write),
        // as well as to keep values used in CS::end(), and may be used for _termination.
        // Update it here.
        _refMegaIteration = std::make_shared<NOMAD::CSMegaIteration>(this, k, barrier, mesh, megaIterSuccess);
        
        NOMAD::CSMegaIteration megaIteration(this, k, barrier, mesh, megaIterSuccess);
        while (!_termination->terminate(k))
        {
            megaIteration.start();
            megaIteration.run();
            megaIteration.end();
            
            // Counter is incremented when calling mega iteration end()
            k       = megaIteration.getK();
            megaIterSuccess = megaIteration.getSuccessType();
            
            if (!successFound && megaIterSuccess >= NOMAD::SuccessType::FULL_SUCCESS)
            {
                successFound = true;
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
    
    return successFound;
}


void NOMAD::CS::hotRestartOnUserInterrupt()
{
    if (_stopReasons->checkTerminate())
    {
        return;
    }
#ifdef TIME_STATS
    if (isRootAlgo())
    {
        _totalCPUAlgoTime += NOMAD::Clock::getCPUTime() - _startTime;
    }
#endif // TIME_STATS
    hotRestartBeginHelper();

    // Reset mesh because parameters have changed.
    std::stringstream ss;
    const NOMAD::Iteration* iteration = getParentOfType<NOMAD::Iteration*>();
    if (nullptr != iteration)
    {
        auto mesh = getIterationMesh();
        ss << *mesh;
        // Reset pointer
        mesh.reset();
        mesh = std::make_shared<NOMAD::CSMesh>(iteration->getPbParams());
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



void NOMAD::CS::readInformationForHotRestart()
{
    // Restart from where we were before.
    // For this, we need to read some files.
    // Note: Cache file is treated independently of hot restart file.

    if (_runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES"))
    {
        // Verify the files exist and are readable.
        std::string hotRestartFile = _runParams->getAttributeValue<std::string>("HOT_RESTART_FILE");
        if (NOMAD::checkReadFile(hotRestartFile))
        {
            std::string s = "Read hot restart file " + hotRestartFile;
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_NORMAL);

            // Create a CSMesh and an MadsMegaIteration with default values, to be filled
            // by istream is.
            // NOTE: Working in full dimension
            auto barrier = std::make_shared<NOMAD::ProgressiveBarrier>(NOMAD::INF, NOMAD::Point(_pbParams->getAttributeValue<size_t>("DIMENSION")), NOMAD::EvalType::BB);
            std::shared_ptr<NOMAD::MeshBase> mesh = std::make_shared<NOMAD::CSMesh>(_pbParams);

            _refMegaIteration = std::make_shared<NOMAD::CSMegaIteration>(this, 0, barrier, mesh, NOMAD::SuccessType::UNDEFINED);

            // Here we use Algorithm::operator>>
            NOMAD::read<NOMAD::CS>(*this, hotRestartFile);
        }
    }
}
