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

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/VNSSearchMethod.hpp"
#include "../../Algos/VNSMads/VNS.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"
//
// Reference: File VNS_Search.cpp in NOMAD 3.9.1
// Author: Christophe Tribes


// Initialize static members (this is not very good, see Issue #604)
size_t NOMAD::VNSSearchMethod::_bbEvalByVNS=0;
size_t NOMAD::VNSSearchMethod::_nbVNSSearchRuns=0;
NOMAD::Point NOMAD::VNSSearchMethod::_refFrameCenter= NOMAD::Point();

void NOMAD::VNSSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_VNS_MADS);
    verifyParentNotNull();

    const auto parentSearch = getParentStep()->getParentOfType<NOMAD::VNSSearchMethod*>(false);
    
    // Do not perform if EVAL_SURROGATE_OPTIMIZATION is true
    bool bBEval = ( NOMAD::EvcInterface::getEvaluatorControl()->getEvalType() == EvalType::BB ) ;
    
    setEnabled((nullptr == parentSearch) && _runParams->getAttributeValue<bool>("VNS_MADS_SEARCH") && bBEval);

}

void NOMAD::VNSSearchMethod::reset()
{
    _refFrameCenter = NOMAD::Point();
    _bbEvalByVNS = 0;
    _nbVNSSearchRuns = 0;
}

bool NOMAD::VNSSearchMethod::runImp()
{
    bool foundBetter = false;
    
    _nbVNSSearchRuns++;
    
    if (isEnabled())
    {
        // check the VNS_trigger criterion:
        size_t bbEval = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
        auto trigger = _runParams->getAttributeValue<NOMAD::Double>("VNS_MADS_SEARCH_TRIGGER");
        if (bbEval == 0 || double(_bbEvalByVNS)/bbEval < trigger.todouble())
        {
            
            std::shared_ptr<EvalPoint> frameCenter = nullptr;
            
            NOMAD::EvcInterface::getEvaluatorControl()->incrementVNSMadsNeighParameter();
            
            // Frame center for VNS Mads
            std::shared_ptr<NOMAD::Barrier> barrier = nullptr;
            
            // Get barrier from upper MadsMegaIteration, if available.
            auto madsMegaIter = getParentOfType<NOMAD::MadsMegaIteration*>(false);
            if (nullptr != madsMegaIter)
            {
                barrier = madsMegaIter->getBarrier();
            }
            else
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"VNS Mads needs a barrier");
            }
            
            // MegaIteration's barrier member is already in sub dimension.
            auto bestXFeas = barrier->getFirstXFeas();
            auto bestXInf  = barrier->getFirstXInf();
            
            auto evalType = NOMAD::EvcInterface::getEvaluatorControl()->getEvalType();
            auto computeType = NOMAD::EvcInterface::getEvaluatorControl()->getComputeType();
            if (nullptr != bestXFeas
                && bestXFeas->getF(evalType, computeType).isDefined()
                && bestXFeas->getF(evalType, computeType) < MODEL_MAX_OUTPUT)
            {
                frameCenter = bestXFeas;
            }
            else if (nullptr != bestXInf
                     && bestXInf->getF(evalType, computeType).isDefined()
                     && bestXInf->getF(evalType, computeType) < MODEL_MAX_OUTPUT
                     && bestXInf->getH(evalType, computeType).isDefined()
                     && bestXInf->getH(evalType, computeType) < MODEL_MAX_OUTPUT)
            {
                frameCenter = bestXInf;
            }
            
            
            if ( nullptr != frameCenter )
            {
                if ( !_refFrameCenter.isDefined() || *(frameCenter->getX()) != _refFrameCenter)
                {
                    NOMAD::EvcInterface::getEvaluatorControl()->resetVNSMadsNeighParameter();
                    _refFrameCenter = *(frameCenter->getX());
                }
                NOMAD::EvcInterface::getEvaluatorControl()->incrementVNSMadsNeighParameter();
                                
                // VNS is an algorithm with its own stop reasons.
                auto vnsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::VNSStopType>>();
                
                // Create the NM algorithm with its own stop reason
                NOMAD::VNS vnsAlgo (this,
                                    vnsStopReasons ,
                                    _runParams,
                                    _pbParams,
                                    frameCenter);
                vnsAlgo.setEndDisplay(false);
                
                vnsAlgo.start();
                vnsAlgo.run();
                vnsAlgo.end();
                
                // Get the success type and update Mads barrier with VNS Mads barrier
                auto vnsBarrier = vnsAlgo.getBarrier();
                
                if (nullptr != vnsBarrier)
                {
                    auto vnsBestFeas = vnsBarrier->getFirstXFeas();
                    auto vnsBestInf = vnsBarrier->getFirstXInf();
                    NOMAD::SuccessType success = barrier->getSuccessTypeOfPoints(vnsBestFeas,
                                                                                 vnsBestInf,
                                                                                 NOMAD::EvalType::BB,
                                                                                 NOMAD::ComputeType::STANDARD);
                    setSuccessType(success);
                    if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
                    {
                        foundBetter = true;
                    }
                    
                    // Update the barrier
                    barrier->updateWithPoints(vnsBarrier->getAllPoints(),
                                                                    NOMAD::EvalType::BB,
                                                                    NOMAD::ComputeType::STANDARD,
                                                                    _runParams->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"));
                    
                }
                
                size_t bbEvalAfterVNS = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
                
                _bbEvalByVNS += bbEvalAfterVNS - bbEval;
                
            }
            
        }
    }
    return foundBetter;
}


void NOMAD::VNSSearchMethod::generateTrialPointsImp()
{
    std::string s;
    NOMAD::EvalPointSet trialPoints;

    // The trial points of one iteration of VNS are generated (not evaluated).
    // The trial points are obtained by shuffle + mads poll

    // auto madsIteration = getParentOfType<MadsIteration*>();

    /*
    // Note: Use first point of barrier as simplex center.
    NOMAD::VNSSingle singleVNS(this,
                            std::make_shared<NOMAD::EvalPoint>(getMegaIterationBarrier()->getFirstPoint()),
                            madsIteration->getMesh());
    singleVNS.start();
    singleVNS.end();

    // Pass the generated trial pts to this
    auto trialPts = singleVNS.getTrialPoints();
    for (auto point : trialPts)
    {
        insertTrialPoint(point);
    }
     */

}   // end generateTrialPoints


