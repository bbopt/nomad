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

#include "../Algos/EvcInterface.hpp"
#include "../Algos/QuadModel/QuadModelEvaluator.hpp"
#include "../Algos/AlgoStopReasons.hpp"
#include "../Algos/SubproblemManager.hpp"
#include "../Algos/SurrogateEvaluation.hpp"
#include "../Algos/SurrogateEvaluator.hpp"
#include "../Output/OutputQueue.hpp"

void NOMAD::SurrogateEvaluation::init()
{

    if (EvalType::SURROGATE == _evalType)
    {
        setStepType(NOMAD::StepType::SURROGATE_EVALUATION);
    }
    else if (EvalType::MODEL == _evalType)
    {
        setStepType(NOMAD::StepType::MODEL_EVALUATION);
    }
    verifyParentNotNull();
    
    // For now, no need to reset evaluator. SurrogateEvaluation is build from new every time.
}


void NOMAD::SurrogateEvaluation::startImp()
{
     auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    if (EvalType::SURROGATE == _evalType)
    {
        _evaluator = std::make_shared<NOMAD::SurrogateEvaluator>(evc->getEvalParams());
    }
    else if(EvalType::MODEL == _evalType)
    {
        auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");

        auto fullFixedVar = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this);
        OUTPUT_INFO_START
        std::string s = "Create QuadModelEvaluator with fixed variable = ";
        s += fullFixedVar.display();
        AddOutputInfo(s);
        OUTPUT_INFO_END

        _quadModelIteration = std::make_unique<NOMAD::QuadModelIteration>(_parentStep, _frameCenter, 0, nullptr, _trialPoints /* trial points are used for sorting. */);
        _quadModelIteration->start();
        
        // No Run of QuadModelIteration in this case.
        
        // QuadModelIteration is used for holding the model used to create the evaluator
        auto model = _quadModelIteration->getModel();
        if (nullptr!=model && model->is_ready())
        {
            _evaluator = std::make_shared<NOMAD::QuadModelEvaluator>(evc->getEvalParams(),
                                                                     model,
                                                                     modelDisplay,
                                                                     fullFixedVar);
        }
        _quadModelIteration->end();
        
    }

    
}


bool NOMAD::SurrogateEvaluation::runImp()
{
    // Evaluation using static surrogate or model. The evaluation will be used for sorting afterwards.
    // Setup evaluation for SURROGATE or MODEL:
    //  - Set opportunistic evaluation to false
    //  - Set the Evaluator to SURROGATE or MODEL
    //  - The point's evalType will be set to SURROGATE or MODEL in evalTrialPoints().
    // Evaluate the points using the surrogate or model
    // Reset for BB:
    //  - Reset opportunism
    //  - Reset Evaluator
    // And proceed - the sort using surrogate or model will be done afterwards.
    
    if(nullptr == _evaluator)
    {
        OUTPUT_INFO_START
        std::string s = "Evaluator is not ready for evaluation";
        AddOutputInfo(s);
        OUTPUT_INFO_END
        return false;
    }
    
    NOMAD::EvcInterface evcInterface(_parentStep);
    
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    auto previousOpportunism = evc->getOpportunisticEval();
    evc->setOpportunisticEval(false);
    
    // Replace the EvaluatorControl's evaluator with this one
    // we just created
    auto previousEvaluator = evc->setEvaluator(_evaluator);
    if (nullptr == previousEvaluator)
    {
        std::cerr << "Warning: Could not set " << evalTypeToString(_evalType) << " Evaluator" << std::endl;
        return false;
    }
    
    // No barrier need for surrogate (MODEL or SURROGATE) evaluations for sort because not comparison is neeeded to detect success
    evc->setBarrier(nullptr);
    
    evc->lockQueue();
    
    evcInterface.keepPointsThatNeedEval(_trialPoints, false /* do not use the mesh */);
    
    // The number of points in the queue for evaluation
    // Watch out, this number must be fetched before unlocking the queue,
    // otherwise there is a risk that the evaluations already restarted
    // and the queue may be empty.
    size_t nbEvalPointsThatNeedEval = evc->getQueueSize(NOMAD::getThreadNum());
    
    // Arguments to unlockQueue:
    // true: do sort
    // keepN: keep a maximum of N points
    // removeStepType: only remove points of this StepType.
    evc->unlockQueue(true, true, getStepType());
    
    if (nbEvalPointsThatNeedEval > 0)
    {
        evcInterface.startEvaluation();
        
        
        // Update trial points with evaluated trial points.
        // Note: If cache is not used, Points that are not evaluated yet
        // will be forgotten.
        NOMAD::EvalPointSet evalPointSet;
        for (auto evalPoint : evcInterface.retrieveAllEvaluatedPoints())
        {
            evalPointSet.insert(evalPoint);
        }
        OUTPUT_DEBUG_START
        std::string s;
        s = "Number of trial points: " + std::to_string(_trialPoints.size());
        _parentStep->AddOutputDebug(s);
        s = "Number of trial points that needed eval: " + std::to_string(nbEvalPointsThatNeedEval);
        _parentStep->AddOutputDebug(s);
        s = "Number of evaluated points: " + std::to_string(evalPointSet.size());
        _parentStep->AddOutputDebug(s);
        OUTPUT_DEBUG_END
        
        _trialPoints.clear();
        _trialPoints = evalPointSet;
    }
    
    
    evc->setEvaluator(previousEvaluator);
    evc->setOpportunisticEval(previousOpportunism);
    
    // Points are now all evaluated using SURROGATE or MODEL.
    // Points are still in the parent step's trial points. There is no update to do here.
    // They are ready to be sorted using their surrogate values, and then
    // evaluated using the BB.
    
    return true;
}


void NOMAD::SurrogateEvaluation::endImp()
{
}
