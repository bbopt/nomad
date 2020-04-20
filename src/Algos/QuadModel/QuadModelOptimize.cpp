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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/QuadSearchMethod.hpp"
#include "../../Algos/QuadModel/QuadModelEvaluator.hpp"
#include "../../Algos/QuadModel/QuadModelOptimize.hpp"
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"


void NOMAD::QuadModelOptimize::init()
{
    _name = getAlgoName() + "Optimize";
    verifyParentNotNull();

    if (nullptr == _iterAncestor )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,_name + " must have an Iteration ancestor.");
    }
    
}


void NOMAD::QuadModelOptimize::startImp()
{
    
    
    auto modelDisplay = _runParams->getAttributeValue<std::string>("MODEL_DISPLAY");
    _displayLevel = (std::string::npos != modelDisplay.find("O"))
                        ? NOMAD::OutputLevel::LEVEL_INFO
                        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

    OUTPUT_INFO_START
    std::string s;
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlParams();
    s = "MAX_SGTE_EVAL: " + std::to_string(evcParams->getAttributeValue<size_t>("MAX_SGTE_EVAL"));
    AddOutputInfo(s, _displayLevel);
    s = "BBOT: " + NOMAD::BBOutputTypeListToString(NOMAD::QuadModelAlgo::getBBOutputType());
    AddOutputInfo(s, _displayLevel);
    OUTPUT_INFO_END
    
    generateTrialPoints();

}

bool NOMAD::QuadModelOptimize::runImp()
{
    std::string s;
    bool foundBetter = false;
    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
        
        // Update barrier
        postProcessing(getEvalType());
        
        // If the oracle point cannot be evaluated the optimization has failed.
        if (_success==NOMAD::SuccessType::NOT_EVALUATED)
        {
            auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );
            qmsStopReason->set( NOMAD::ModelStopType::NO_NEW_POINTS_FOUND);
        }
        
    }

    
    return foundBetter;
}

void NOMAD::QuadModelOptimize::endImp()
{
    // Clean up the cache of points having only EvalType::SGTE
    NOMAD::CacheBase::getInstance()->deleteSgteOnly();
}

void NOMAD::QuadModelOptimize::setupRunParameters()
{
    _optRunParams = std::make_shared<NOMAD::RunParameters>(*_runParams);
    
    _optRunParams->setAttributeValue("MAX_ITERATIONS", INF_SIZE_T);
    
    // Ensure there is no model used in model optimization.
    _optRunParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    _optRunParams->setAttributeValue("SGTELIB_SEARCH", false);
    _optRunParams->setAttributeValue("NM_SEARCH", false);
    
    // No hMax in the context of QuadModel
    _optRunParams->setAttributeValue("H_MAX_0", NOMAD::Double(NOMAD::INF));

    // Disable user calls
    _optRunParams->setAttributeValue("USER_CALLS_ENABLED", false);

    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlParams();

    _optRunParams->checkAndComply(evcParams, _pbParams);
    
}


void NOMAD::QuadModelOptimize::setupPbParameters()
{
    _optPbParams = std::make_shared<NOMAD::PbParameters>(*_refPbParams);
    _optPbParams->setAttributeValue("LOWER_BOUND", _modelLowerBound);
    _optPbParams->setAttributeValue("UPPER_BOUND", _modelUpperBound);
    _optPbParams->setAttributeValue("FIXED_VARIABLE",_modelFixedVar);
    
    NOMAD::ArrayOfPoint x0s;
    auto frameCenter = _iterAncestor->getFrameCenter();
    if (frameCenter->inBounds(_modelLowerBound, _modelUpperBound))
    {
        x0s.push_back(*(frameCenter->getX()));
    }
    else
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"A frameCenter must be available and within bounds to set X0 for quad optimization.");
    }
    _optPbParams->setAttributeValue("X0", x0s);

    // Ensure that optimization will be done using SGTE evals. All sub
    // steps will use SGTE.
    _optPbParams->setAttributeValue("EVAL_TYPE", NOMAD::EvalType::SGTE);

    // We do not want certain warnings appearing in sub-optimization.
    _optPbParams->doNotShowWarnings();

    _optPbParams->checkAndComply();

}

void NOMAD::QuadModelOptimize::generateTrialPoints()
{
    // Clear the previous trial points
    clearTrialPoints();
    
    // Set specific evaluator control
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlParams();
    _previousEvcParamOpportunism = evcParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL");
    _previousEvcParamBbMaxBlockSize = evcParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
    int previousNbThreadsOpenMP = _runParams->getAttributeValue<int>("NB_THREADS_OPENMP");

    evcParams->setAttributeValue("BB_MAX_BLOCK_SIZE", INF_SIZE_T );
    evcParams->setAttributeValue("OPPORTUNISTIC_EVAL", false);
    evcParams->checkAndComply();
    _runParams->setAttributeValue("NB_THREADS_OPENMP", 1);
    _runParams->checkAndComply(evcParams, _pbParams);
    
    auto evalParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvalParams();
    // Enforce no opportunism and use no cache.

    auto previousOpportunism = evcParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL");
    auto previousUseCache = evcParams->getAttributeValue<bool>("USE_CACHE");
    evcParams->setAttributeValue("OPPORTUNISTIC_EVAL", false);
    evcParams->setAttributeValue("USE_CACHE", false);
    evcParams->checkAndComply();

    auto modelDisplay = _runParams->getAttributeValue<std::string>("MODEL_DISPLAY");
    auto ev = std::make_unique<NOMAD::QuadModelEvaluator>(evalParams,
                                                          _model,
                                                          modelDisplay,
                                                          getSubFixedVariable());
    
    // Replace the EvaluatorControl's evaluator with this one
    // we just created
    auto mainEvaluator = EvcInterface::getEvaluatorControl()->getEvaluatorUPtr();
    EvcInterface::getEvaluatorControl()->setEvaluator(std::move(ev));
    
    // Setup EvalPoint success computation to be based on sgte rather than bb.
    NOMAD::ComputeSuccessType::setComputeSuccessTypeFunction(
                                NOMAD::ComputeSuccessType::computeSuccessTypeSgte);

    
    auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
    
    // Set the bounds and fixed variables from the model
    setModelBoundsAndFixedVar();
    
    
    // Set and verify run parameter values
    setupRunParameters();
    
    // Setup Pb parameters just before optimization.
    setupPbParameters();
    
    OUTPUT_INFO_START
    std::ostringstream oss;
    oss << "Run Parameters for QuadModelOptimize:" << std::endl;
    _optRunParams->display(oss, false);
    AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
    OUTPUT_INFO_END
    
    
    // Create a Mads step
    // Parameters for mads (_optRunParams and _optPbParams) are already updated.
    auto mads = std::make_shared<NOMAD::Mads>(this, madsStopReasons, _optRunParams, _optPbParams);
    mads->setName(mads->getName() + " (QuadModelOptimize)");
    EvcInterface::getEvaluatorControl()->resetSgteEval();
    NOMAD::MainStep::setAlgoComment("(QuadModelOptimize)");
    mads->start();
    bool runOk = mads->run();
    mads->end();
    NOMAD::MainStep::resetPreviousAlgoComment();
    
    // Note: No need to check the Mads stop reason: It is not a stop reason
    // for QuadModel.
    
    // Reset opportunism to previous values.
    evcParams->setAttributeValue("OPPORTUNISTIC_EVAL", previousOpportunism);
    evcParams->setAttributeValue("USE_CACHE", previousUseCache);
    evcParams->checkAndComply();
    
    // When we are done, restore mainEvaluator
    EvcInterface::getEvaluatorControl()->setEvaluator(std::move(mainEvaluator));
    
    // Reset success computation function --> use the default
    NOMAD::ComputeSuccessType::setComputeSuccessTypeFunction(
                                            NOMAD::ComputeSuccessType::defaultComputeSuccessType);
    
    evcParams->setAttributeValue("BB_MAX_BLOCK_SIZE", _previousEvcParamBbMaxBlockSize );
    evcParams->setAttributeValue("OPPORTUNISTIC_EVAL", _previousEvcParamOpportunism);
    evcParams->checkAndComply();
    _runParams->setAttributeValue("NB_THREADS_OPENMP", previousNbThreadsOpenMP);
    _runParams->checkAndComply(evcParams, _pbParams);
    
    
    if (!runOk)
    {
        auto modelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        modelStopReasons->set(NOMAD::ModelStopType::MODEL_OPTIMIZATION_FAIL);
    }
    else
    {
        auto bestXFeas = mads->getMegaIterationBarrier()->getFirstXFeas();
        auto bestXInf  = mads->getMegaIterationBarrier()->getFirstXInf();
        if (nullptr != bestXFeas)
        {
            // New EvalPoint to be evaluated.
            // Add it to the list (local or in Search method).
            bool inserted = insertTrialPoint(NOMAD::EvalPoint(*bestXFeas));
            
            OUTPUT_INFO_START
            std::string s = "xt:";
            s += (inserted) ? " " : " not inserted: ";
            s += bestXFeas->display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
        if (nullptr != bestXInf)
        {
            // New EvalPoint to be evaluated.
            // Add it to the lists (local or in Search method).
            insertTrialPoint(NOMAD::EvalPoint(*bestXInf));
            
        }
    }
}


// Set the bounds and the extra fixed variables (when lb=ub) of the model.
void NOMAD::QuadModelOptimize::setModelBoundsAndFixedVar()
{
    const SGTELIB::Matrix & X = _trainingSet->get_matrix_X();
    
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    if (n != (size_t)X.get_nb_cols())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "QuadModel::setModelBounds() dimensions do not match");
    }

    int nbDim = X.get_nb_cols();
    int nbPoints = X.get_nb_rows();

    // Build model bounds & detect fixed variables
    NOMAD::Double lb;
    NOMAD::Double ub;

    for (int j = 0; j < nbDim; j++)
    {
        lb = _modelLowerBound[j];
        ub = _modelUpperBound[j];
        for (int p = 0; p < nbPoints; p++)
        {
            lb = NOMAD::min(lb, NOMAD::Double(X.get(p,j)));
            ub = NOMAD::max(ub, NOMAD::Double(X.get(p,j)));
        }
        if (lb==ub)
        {
            _modelFixedVar[j] = ub;
            lb = -INF;
            ub = +INF;
        }
            
        _modelLowerBound[j] = lb;
        _modelUpperBound[j] = ub;

    }
    OUTPUT_INFO_START
    std::string s = "model lower bound:";
    s += _modelLowerBound.display();
    s += "model upper bound:";
    s += _modelUpperBound.display();
    AddOutputInfo(s);
    OUTPUT_INFO_END

} // end setModelBounds
