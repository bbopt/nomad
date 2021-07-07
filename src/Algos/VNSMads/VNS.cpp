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
#include "../../Algos/Mads/SinglePollMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/LHSearchType.hpp"
#include "../../Util/fileutils.hpp"

// VNS specific
#include "../../Algos/VNSMads/VNS.hpp"


void NOMAD::VNS::init()
{
    /*
    if ( _runParams->getAttributeValue<bool>("MEGA_SEARCH_POLL") )
    {
        _name += " One Iteration";
    }
     */
    setStepType(NOMAD::StepType::ALGORITHM_VNS_MADS);
    
    if (nullptr == _frameCenter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"VNS Mads needs a frame center");
    }

}


void NOMAD::VNS::startImp()
{
    // All stop reasons are reset.
    _stopReasons->setStarted();

    // Reset the lapBbEval counter for this sub-optimization
    NOMAD::EvcInterface::getEvaluatorControl()->resetLapBbEval();

    // Setup Mads
    _madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
    
}

bool NOMAD::VNS::runImp()
{
    _algoSuccessful = false;
    
    _algoBestSuccess = NOMAD::SuccessType::NOT_EVALUATED;
    
    auto VNSStopReasons = NOMAD::AlgoStopReasons<NOMAD::VNSStopType>::get(_stopReasons);
    
    if ( _stopReasons->checkTerminate() )
    {
        return _algoSuccessful;
    }
    
    if (_runParams->getAttributeValue<bool>("VNS_MADS_OPTIMIZATION"))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"VNS_MADS_OPTIMIZATION not yet implemented");
    }
    
    // Get the parent Mads Mega iteration and its associated mesh
    auto parentMadsMegaIter = getParentOfType<NOMAD::MadsMegaIteration*>(false);
    if (nullptr == parentMadsMegaIter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"VNS Mads needs a MadsMegaIteration");
    }
    
    // Shaking direction: use single poll method
    NOMAD::SinglePollMethod pollMethod(this, *_frameCenter);
    std::list<NOMAD::Direction> scaledDirection = pollMethod.generateFullSpaceScaledDirections(false,parentMadsMegaIter->getMesh());
    
    // Get the single direction
    NOMAD::Direction dir = scaledDirection.front();
    if (!dir.isDefined())
    {
      throw NOMAD::Exception(__FILE__,__LINE__,"VNS_MADS_OPTIMIZATION: single scaled direction not defined");
    }
    
    // Multiply shake direction by VNS neighborhood parameter
    dir *= (double)NOMAD::EvcInterface::getEvaluatorControl()->getVNSMadsNeighParameter();
    
    OUTPUT_INFO_START
    AddOutputInfo("Shaking direction: " + dir.display());
    OUTPUT_INFO_END
    
    // shaking: the perturbation is tried twice with dir and -dir
    //          (in case x == x + dir after snapping)
    NOMAD::Point shakePoint = *(_frameCenter->getX()) + dir; // pun intended;
    shakePoint.snapToBounds(_pbParams->getAttributeValue<ArrayOfDouble>("LOWER_BOUND"),_pbParams->getAttributeValue<ArrayOfDouble>("UPPER_BOUND"));
    for ( int nbt = 0 ; nbt < 2 ; ++nbt )
    {
        if ( shakePoint == *(_frameCenter->getX()))
        {
            // no third try: the search fails
            if ( nbt == 1 )
            {
                OUTPUT_INFO_START
                AddOutputInfo("VNS: Shaking failed");
                OUTPUT_INFO_END
                
                VNSStopReasons->set(NOMAD::VNSStopType::SHAKING_FAILED);
                
                return _algoSuccessful;
            }
            
            // 2nd try (-dir instead of dir):
            shakePoint =  *(_frameCenter->getX()) - dir;
            shakePoint.snapToBounds(_pbParams->getAttributeValue<ArrayOfDouble>("LOWER_BOUND"),_pbParams->getAttributeValue<ArrayOfDouble>("UPPER_BOUND"));
            
        }
    }
    
    // Generate base point (X0 for Mads)
    auto evalPoint = NOMAD::EvalPoint(shakePoint);
    
    // Meta information on eval point
    evalPoint.setPointFrom(_frameCenter, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
    evalPoint.addGenStep(getStepType());

    // Get the current Mads frame size to be used as min frame size for suboptimization
    NOMAD::ArrayOfDouble currentMadsFrameSize = parentMadsMegaIter->getMesh()->getDeltaFrameSize();
    
    setupPbParameters(shakePoint,currentMadsFrameSize);
    setupRunParameters();
    
    NOMAD::Mads mads(this, _madsStopReasons, _optRunParams, _optPbParams, false /*false: Barrier not initialized from cache */ );
    
    // Run Mads.
    mads.start();
    mads.run();
    mads.end();
    
    if ( _madsStopReasons->testIf(NOMAD::MadsStopType::X0_FAIL) )
    {
        VNSStopReasons->set(NOMAD::VNSStopType::X0_FAILED);
    }
    else
    {
        _barrier = mads.getMegaIterationBarrier();
        _algoSuccessful = true;
        _algoBestSuccess = NOMAD::SuccessType::FULL_SUCCESS;
    }
    
    
    _termination->start();
    _termination->run();
    _termination->end();
    
    return _algoSuccessful;
}




void NOMAD::VNS::endImp()
{
    NOMAD::Algorithm::endImp();


}


void NOMAD::VNS::setupRunParameters()
{
    _optRunParams = std::make_shared<NOMAD::RunParameters>(*_runParams);

    _optRunParams->setAttributeValue("MAX_ITERATIONS", INF_SIZE_T);

    // VNS do not perform VNS search
    _optRunParams->setAttributeValue("VNS_MADS_SEARCH", false);
    
    // No LH search
    _optRunParams->setAttributeValue("LH_SEARCH", NOMAD::LHSearchType("0 0"));
    
    
    auto vnsFactor = _runParams->getAttributeValue<size_t>("VNS_MADS_SEARCH_MAX_TRIAL_PTS_NFACTOR");
    auto dim = _pbParams->getAttributeValue<size_t>("DIMENSION");
    if (vnsFactor < NOMAD::INF_SIZE_T)
    {
        NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( dim*vnsFactor );
    }
    
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();

    _optRunParams->checkAndComply(evcParams, _optPbParams);
}


void NOMAD::VNS::setupPbParameters(const NOMAD::Point & center, const NOMAD::ArrayOfDouble & currentMadsFrameSize)
{
    _optPbParams = std::make_shared<NOMAD::PbParameters>(*_pbParams);

    // Reset initial mesh and frame sizes
    // The initial mesh and frame sizes will be calculated from bounds and X0
    _optPbParams->resetToDefaultValue("INITIAL_MESH_SIZE");
    _optPbParams->resetToDefaultValue("INITIAL_FRAME_SIZE");
    
    // set min frame size (min mesh size will be updated)
    _optPbParams->resetToDefaultValue("MIN_MESH_SIZE");
    _optPbParams->resetToDefaultValue("MIN_FRAME_SIZE");
    _optPbParams->setAttributeValue("MIN_FRAME_SIZE", currentMadsFrameSize);

    NOMAD::ArrayOfPoint x0s{center};
    _optPbParams->setAttributeValue("X0", x0s);

    // We do not want certain warnings appearing in sub-optimization.
    _optPbParams->doNotShowWarnings();

    _optPbParams->checkAndComply();

}
