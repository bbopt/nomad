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
/**
 \file   DiscoMads.cpp
 \brief  The DiscoMads algorithm (main): implementation
 \author Solene Kojtych
 \see    DiscoMads.hpp
 */

#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/DiscoMads/DiscoMads.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/DiscoMads/DiscoMadsMegaIteration.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/fileutils.hpp"
#ifdef TIME_STATS
#include "../../Util/Clock.hpp"
#endif

void NOMAD::DiscoMads::init(bool barrierInitializedFromCache)
{
    setStepType(NOMAD::StepType::ALGORITHM_DISCO_MADS);
    
    verifyParentNotNull();

    // Instantiate Mads initialization class
    _initialization = std::make_unique<NOMAD::MadsInitialization>( this, barrierInitializedFromCache, false /*initialization for DMultiMads*/, true /*initialization for DiscodMAds*/ );



    // -- Display discoMads parameters
    bool detectHiddConst = _runParams->getAttributeValue<bool>("DISCO_MADS_HID_CONST");                   // only for hidden constraints revaluation

    if(detectHiddConst)
    {   // use to reveal hidden constraints
        NOMAD::Double highValue = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_HID_CONST_OUTPUT_VALUE");   
        OUTPUT_INFO_START
            AddOutputInfo("DiscoMads used to reveal hidden constraints.",true,false);
            AddOutputInfo("Value attributed to OBJ/PB output of failed evaluations: "+highValue.tostring());
        OUTPUT_INFO_END
    }
    else{
        // use to reveal discontinuities
        NOMAD::Double detectionRadius = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_DETECTION_RADIUS"); 
        NOMAD::Double limitRate = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_LIMIT_RATE");
        OUTPUT_INFO_START
            AddOutputInfo("DiscoMads used to reveal discontinuities.",true,false);
            AddOutputInfo("Discontinuities characterized by detection radius "+detectionRadius.tostring()+" and limit rate "+limitRate.tostring());
        OUTPUT_INFO_END
    }
        // Common parameters
    NOMAD::Double exclusionRadius = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_EXCLUSION_RADIUS");
    size_t revealingPollnbPoints = _runParams->getAttributeValue<size_t>("DISCO_MADS_REVEALING_POLL_NB_POINTS"); 
    NOMAD::Double revealingPollRadius = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_REVEALING_POLL_RADIUS");    
    OUTPUT_INFO_START
        AddOutputInfo("Exclusion radius: "+exclusionRadius.tostring());
        AddOutputInfo("Revealing poll: "+to_string(revealingPollnbPoints)+" points, radius = "+revealingPollRadius.tostring());
        AddOutputInfo("",false,true);
    OUTPUT_INFO_END    


}




bool NOMAD::DiscoMads::runImp()
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
            mesh    = (std::dynamic_pointer_cast<NOMAD::DiscoMadsMegaIteration> (_refMegaIteration ))->getMesh();
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
        _refMegaIteration = std::make_shared<NOMAD::DiscoMadsMegaIteration>(this, k, barrier, mesh, megaIterSuccess);

        
            
        // Create a DiscoMegaIteration for looping: manage multiple iterations on different
        // meshes and with different frame centers at the same time.
        NOMAD::DiscoMadsMegaIteration megaIteration(this, k, barrier, mesh, megaIterSuccess);
       
        
        
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

        }
    }
    
    _termination->start();
    _termination->run();
    _termination->end();

    return _algoSuccessful;
}


void NOMAD::DiscoMads::hotRestartOnUserInterrupt()
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
        
        mesh = std::make_shared<NOMAD::GMesh>(iteration->getPbParams(),iteration->getRunParams());
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

void NOMAD::DiscoMads::readInformationForHotRestart()
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

            // Create a GMesh and an DiscoMadsMegaIteration with default values, to be filled
            // by istream is.
            // NOTE: Working in full dimension
            auto barrier = std::make_shared<NOMAD::ProgressiveBarrier>(NOMAD::INF, NOMAD::Point(_pbParams->getAttributeValue<size_t>("DIMENSION")), NOMAD::EvalType::BB);
            
            std::shared_ptr<NOMAD::MeshBase> mesh = std::make_shared<NOMAD::GMesh>(_pbParams,_runParams);

            _refMegaIteration = std::make_shared<NOMAD::DiscoMadsMegaIteration>(this, 0, barrier, mesh, NOMAD::SuccessType::UNDEFINED);

            // Here we use Algorithm::operator>>
            NOMAD::read<NOMAD::DiscoMads>(*this, hotRestartFile);
        }
    }
}
