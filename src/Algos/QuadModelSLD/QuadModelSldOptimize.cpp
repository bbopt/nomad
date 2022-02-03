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
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldAlgo.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldEvaluator.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldIteration.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldOptimize.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"
#include "../../Type/EvalSortType.hpp"

void NOMAD::QuadModelSldOptimize::init()
{
    setStepType(NOMAD::StepType::OPTIMIZE);
    verifyParentNotNull();

    if (nullptr == _iterAncestor )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,getName() + " must have an Iteration ancestor.");
    }

}


void NOMAD::QuadModelSldOptimize::startImp()
{

    auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");
    _displayLevel = (std::string::npos != modelDisplay.find("O"))
        ? NOMAD::OutputLevel::LEVEL_INFO
        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

    OUTPUT_INFO_START
    std::string s;
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();
    s = "QUAD_MODEL_MAX_EVAL: " + std::to_string(evcParams->getAttributeValue<size_t>("QUAD_MODEL_MAX_EVAL"));
    AddOutputInfo(s, _displayLevel);
    s = "BBOT: " + NOMAD::BBOutputTypeListToString(NOMAD::QuadModelSldAlgo::getBBOutputType());
    AddOutputInfo(s, _displayLevel);
    OUTPUT_INFO_END

    generateTrialPoints();
    
}


bool NOMAD::QuadModelSldOptimize::runImp()
{
    std::string s;
    bool foundBetter = false;
    
    throw NOMAD::Exception(__FILE__,__LINE__,"Quad Model SLD stand alone optimization not implemented");
    
//    if ( ! _stopReasons->checkTerminate() )
//    {
//        foundBetter = evalTrialPoints(this);
//
//        if (_modelFixedVar.nbDefined() > 0)
//        {
//            NOMAD::EvalPointSet evalPointSet;
//            for (auto trialPoint : _trialPoints)
//            {
//                evalPointSet.insert(trialPoint.makeFullSpacePointFromFixed(_modelFixedVar));
//            }
//            _trialPoints.clear();
//            _trialPoints = evalPointSet;
//        }
//        // Update barrier
//        postProcessing();
//
//        // If the oracle point cannot be evaluated the optimization has failed.
//        if (_success==NOMAD::SuccessType::NOT_EVALUATED)
//        {
//            auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );
//            qmsStopReason->set( NOMAD::ModelStopType::NO_NEW_POINTS_FOUND);
//        }
//
//    }
//
//
    return foundBetter;
}


void NOMAD::QuadModelSldOptimize::endImp()
{
    // Clean up the cache of points having only EvalType::SGTE
    NOMAD::CacheBase::getInstance()->deleteModelEvalOnly(NOMAD::getThreadNum());
    
    
}


void NOMAD::QuadModelSldOptimize::setupRunParameters()
{
    _optRunParams = std::make_shared<NOMAD::RunParameters>(*_runParams);

    _optRunParams->setAttributeValue("MEGA_SEARCH_POLL", false);
    
    _optRunParams->setAttributeValue("MAX_ITERATIONS", INF_SIZE_T);

    // Ensure there is no model, no NM and no VNS used in model optimization.
    _optRunParams->setAttributeValue("QUAD_MODEL_SLD_SEARCH", false);
    _optRunParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    _optRunParams->setAttributeValue("SGTELIB_MODEL_SEARCH", false);
    _optRunParams->setAttributeValue("NM_SEARCH", false);
    _optRunParams->setAttributeValue("SPECULATIVE_SEARCH", true);
    
    _optRunParams->setAttributeValue("ANISOTROPIC_MESH", false);
    
    _optRunParams->setAttributeValue("VNS_MADS_SEARCH", false);

    // Set direction type to Ortho 2n
    _optRunParams->setAttributeValue("DIRECTION_TYPE",NOMAD::DirectionType::ORTHO_2N);

    // No hMax in the context of QuadModel
    _optRunParams->setAttributeValue("H_MAX_0", NOMAD::Double(NOMAD::INF));

    // Disable user calls
    _optRunParams->setAttributeValue("USER_CALLS_ENABLED", false);

    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();

    _optRunParams->checkAndComply(evcParams, _optPbParams);
}


void NOMAD::QuadModelSldOptimize::setupPbParameters()
{
    const NOMAD::QuadModelSldIteration * iter = getParentOfType<NOMAD::QuadModelSldIteration*>();
    
    if (nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Optimize must have a Iteration among its ancestors.");
    }
    
    auto model = iter->getModel();
    
    size_t n = _refPbParams->getAttributeValue<size_t>("DIMENSION");
    if (n != model->get_n())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Dimensions do not match");
    }
    
    // Copy refPbParams for the parameters of this optimization
    _optPbParams = std::make_shared<NOMAD::PbParameters>(*_refPbParams);
    
    // 1/2
    NOMAD::Point X0( n , 0.0 ) ;
    NOMAD::ArrayOfPoint x0s{X0};
    
    // 2/2: model center if different than (0,..,0) and if in [-1;1]:
    NOMAD::Point x1 = model->get_center();
    if ( x1.size() == n && x1.isComplete() )
    {
        
        model->scale ( x1 );
        
        bool diff   = false;
        bool bnd_ok = true;
        
        for (size_t i = 0 ; i < n ; ++i )
        {
            
            if ( x1[i] != 0 )
                diff = true;
            
            if ( x1[i].abs() > 1.0 || model->variable_is_fixed(static_cast<int>(i)))
            {
                bnd_ok = false;
                break;
            }
            x1[i] *= 1000.0;
        }
        
        if ( diff && bnd_ok )
            x0s.push_back(x1);
            
    }
    _optPbParams->setAttributeValue("X0", x0s);
        
    
    // fixed variables:
    NOMAD::Point fixedVar(n);
    for (size_t i = 0 ; i < n ; ++i )
    {
        if ( model->variable_is_fixed(static_cast<int>(i)) )
        {
            fixedVar[i]=0.0;
        }
    }
    _optPbParams->setAttributeValue("FIXED_VARIABLE",fixedVar);

    
    // bounds:
    NOMAD::Point lb ( n , -1.0 );
    NOMAD::Point ub ( n ,  1.0 );
    
    auto LB = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto UB = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    
    if ( LB.isComplete() || UB.isComplete() )
    {
        
        model->unscale ( lb );
        model->unscale ( ub );
        
        for (size_t i = 0 ; i < n ; ++i )
        {
            
            if ( LB[i].isDefined() && LB[i] > lb[i] )
                lb[i] = LB[i];
            
            if ( UB[i].isDefined() && UB[i] < ub[i] )
                ub[i] = UB[i];
        }
        
        model->scale ( lb );
        model->scale ( ub );
        
        for (size_t i = 0 ; i < n ; ++i )
        {
            if ( ub[i] < lb[i] || lb[i] > 0.0 || ub[i] < 0.0 )
            {
                throw NOMAD::Exception(__FILE__, __LINE__,
                                       "Bounds do not match");
            }
            lb[i] *= 1000.0;
            ub[i] *= 1000.0;
        }
    }
    else
    {
        lb *= 1000.0;
        ub *= 1000.0;
    }
    _optPbParams->setAttributeValue("LOWER_BOUND", static_cast<NOMAD::ArrayOfDouble>(lb));
    _optPbParams->setAttributeValue("UPPER_BOUND", static_cast<NOMAD::ArrayOfDouble>(ub));

    
    
    // Reset initial mesh and frame sizes
    // The initial mesh and frame sizes will be calculated from bounds and X0
    _optPbParams->resetToDefaultValue("INITIAL_MESH_SIZE");
    _optPbParams->resetToDefaultValue("INITIAL_FRAME_SIZE");

    // Hardcoded stoping criterion. Will be changed in the future
    _optPbParams->setAttributeValue("MIN_MESH_SIZE",NOMAD::ArrayOfDouble(n,1E-30));
    _optPbParams->setAttributeValue("MIN_FRAME_SIZE",NOMAD::ArrayOfDouble(n,1E-30));
    
    // Granularity is set to 0 and bb_input_type is set to all continuous variables. Candidate points are projected on the mesh before evaluation.
    _optPbParams->resetToDefaultValue("GRANULARITY");
    _optPbParams->resetToDefaultValue("BB_INPUT_TYPE");

    // No variable groups are considered for suboptimization
    _optPbParams->resetToDefaultValue("VARIABLE_GROUP");
    


    // We do not want certain warnings appearing in sub-optimization.
    _optPbParams->doNotShowWarnings();

    _optPbParams->checkAndComply();

}


void NOMAD::QuadModelSldOptimize::generateTrialPointsImp()
{
    // Setup Pb parameters just before optimization.
    setupPbParameters();

    // Set and verify run parameter values
    setupRunParameters();

    OUTPUT_INFO_START
    std::ostringstream oss;
    oss << "Run Parameters for QuadModelOptimize:" << std::endl;
    _optRunParams->display(oss, false);
    AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
    OUTPUT_INFO_END
    
    auto fixedVar = _optPbParams->getAttributeValue<NOMAD::Point>("FIXED_VARIABLE");

    if ( fixedVar.nbDefined() == fixedVar.size() )
    {
        OUTPUT_INFO_START
        std::ostringstream oss;
        oss << "Effective dimension is null. No QuadModelOptimize" << std::endl;
        AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
        OUTPUT_INFO_END

        return;
    }
    // Set specific evaluator control
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    // Enforce no opportunism and use no cache.
    auto previousOpportunism = evc->getOpportunisticEval();
    auto previousUseCache = evc->getUseCache();
    evc->setOpportunisticEval(true);
    evc->setUseCache(false);
    
    // Enforce specific sort type
    auto evcParams = evc->getEvaluatorControlGlobalParams();
    auto previousEvalSortType = evcParams->getAttributeValue<NOMAD::EvalSortType>("EVAL_QUEUE_SORT");
    evcParams->setAttributeValue("EVAL_QUEUE_SORT", NOMAD::EvalSortType::DIR_LAST_SUCCESS);
    evcParams->checkAndComply();

    auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");

    OUTPUT_INFO_START
    std::string s = "Create QuadModelEvaluator with fixed variable = ";
    s += fixedVar.display();
    AddOutputInfo(s);
    OUTPUT_INFO_END

    // Evaluations are in the quad model (local) full space: fixed variables detected when creating the training set (lb==ub) are not modified by the optimizer but evaluator receives points in local full space. The "local full space" is used because fixed variables from the parent problem are not seen/considered.
    auto ev = std::make_shared<NOMAD::QuadModelSldEvaluator>(evc->getEvalParams(),
                                                          _model,
                                                          modelDisplay);

    // Replace the EvaluatorControl's evaluator with this one
    // we just created
    auto previousEvaluator = evc->setEvaluator(ev);
    if (nullptr == previousEvaluator)
    {
        std::cerr << "Warning: QuadModelOptimize: Could not set MODEL Evaluator" << std::endl;
        return;
    }

    auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();


    // Create a Mads step
    // Parameters for mads (_optRunParams and _optPbParams) are already updated.
    // We force mads to use local fixed variable to make sure we only work on fixed variables identified during construction of the training set. The evaluator works on the quad model, we don't need to map to the global full space for evaluation (not BB eval).
    auto mads = std::make_shared<NOMAD::Mads>(this, madsStopReasons, _optRunParams, _optPbParams, true /* use local fixed variables */);
    
    mads->setEndDisplay(true);
    
    evc->resetModelEval();
    
    mads->start();
    bool runOk = mads->run();    
    mads->end();
    
    //NOMAD::OutputQueue::getInstance()->setMaxOutputLevel(OutputLevel::LEVEL_NORMAL);
    
    evc->resetModelEval();
    evc->setEvaluator(previousEvaluator);

    // Note: No need to check the Mads stop reason: It is not a stop reason
    // for QuadModel.

    // Reset opportunism to previous values.
    evc->setOpportunisticEval(previousOpportunism);
    evc->setUseCache(previousUseCache);
    
    // Reset eval sort type to previous values
    evcParams->setAttributeValue("EVAL_QUEUE_SORT", previousEvalSortType);
    evcParams->checkAndComply();

    if (!runOk)
    {
        auto modelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        modelStopReasons->set(NOMAD::ModelStopType::MODEL_OPTIMIZATION_FAIL);
    }
    else
    {
        // Get the best points in their reference dimension
        _bestXFeas = std::make_shared<NOMAD::EvalPoint>(mads->getBestSolution(true));
        _bestXInf  = std::make_shared<NOMAD::EvalPoint>(mads->getBestSolution(false));
        if (_bestXFeas->isComplete())
        {
            // New EvalPoint to be evaluated.
            // Add it to the list (local or in Search method).
            *_bestXFeas *= 0.001;  // Scale from [-1000;1000] to [-1;1]
            bool inserted = insertTrialPoint(*_bestXFeas);

            OUTPUT_INFO_START
            std::string s = "Partly unscaled xt(feas):";
            s += (inserted) ? " " : " not inserted: ";
            s += _bestXFeas->display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
        else
        {
            // Reset shared_ptr to default
            _bestXFeas.reset();
        }
        if (_bestXInf->isComplete())
        {
            // New EvalPoint to be evaluated.
            // Add it to the lists (local or in Search method).
            *_bestXInf *= 0.001;  // Scale from [-1000;1000] to [-1;1]
            bool inserted = insertTrialPoint(*_bestXInf);
            
            OUTPUT_INFO_START
            std::string s = "Partly unscaled xt(infeas):";
            s += (inserted) ? " " : " not inserted: ";
            s += _bestXInf->display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        
        }
        else
        {
            // Reset shared_ptr to default
            _bestXInf.reset();
        }
    }

}



