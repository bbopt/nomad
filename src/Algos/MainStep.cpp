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
 \file   MainStep.cpp
 \brief  Main Step to hold MADS
 \author Viviane Rochon Montplaisir
 \date   June 2018
 */

// Generic
#include "../Algos/AlgoStopReasons.hpp"
#include "../Algos/EvcInterface.hpp"
#include "../Algos/MainStep.hpp"
#include "../Algos/SubproblemManager.hpp"
#include "../Cache/CacheSet.hpp"
#include "../Eval/ProgressiveBarrier.hpp"
#include "../Math/LHS.hpp"
#include "../Math/RNG.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Output/OutputDirectToFile.hpp"
#include "../Util/Clock.hpp"
#include "../Type/EvalSortType.hpp"
#include "../Util/Clock.hpp"
#include "../Util/fileutils.hpp"

// Specific algos
#include "../Algos/LatinHypercubeSampling/LH.hpp"
#include "../Algos/CoordinateSearch/CS.hpp"
#include "../Algos/DMultiMads/DMultiMads.hpp"
#include "../Algos/Mads/Mads.hpp"
#include "../Algos/Mads/MadsIteration.hpp"
#include "../Algos/Mads/Search.hpp"
#include "../Algos/Mads/VNSSearchMethod.hpp"
#include "../Algos/NelderMead/NM.hpp"
#ifdef _OPENMP
#include "../Algos/PSDMads/PSDMads.hpp"
#endif
#ifdef USE_SGTELIB
#include "../Algos/QPSolverAlgo/QPSolverAlgo.hpp"
#include "../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../Algos/SgtelibModel/SgtelibModel.hpp"
#endif
#include "../Algos/SSDMads/SSDMads.hpp"
#include "../Algos/DiscoMads/DiscoMads.hpp"
#include "../Algos/TemplateAlgo/TemplateAlgo.hpp"

// For hardware thread number
#ifdef _OPENMP
#include <thread>
#endif

void NOMAD::MainStep::setAllParameters(const std::shared_ptr<NOMAD::AllParameters> &allParams)
{
    // AllParameters keep the whole bunch of parameters.
    _allParams = allParams;

    // run and pb parameters are members of the Step class and must also be set here.
    _runParams = allParams->getRunParams();
    _pbParams  = allParams->getPbParams();
}


void NOMAD::MainStep::init()
{
    _allParams = std::make_shared<NOMAD::AllParameters>();

    // run and pb parameters are members of the Step class and must also be set here.
    _runParams = _allParams->getRunParams();
    _pbParams  = _allParams->getPbParams();

    setStepType(NOMAD::StepType::MAIN);

    // Start the clock
    NOMAD::Clock::reset();
}


NOMAD::MainStep::~MainStep()
{
    _algos.clear();
}

NOMAD::ArrayOfPoint NOMAD::MainStep::suggest()
{
    NOMAD::ArrayOfPoint suggestedPoints;
    AddOutputInfo("Start step " + getName(), true, false);

    // No X0 should be provided
    auto x0s = _allParams->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
    if (!x0s.empty() && !x0s[0].toBeDefined())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Using suggest with x0 provided. Use x0 instead of calling suggest.");
    }

    // Display attributes and check attribute consistency
    if (_allParams->getAttributeValue<int>("DISPLAY_DEGREE") >= (int)NOMAD::OutputLevel::LEVEL_DEBUG)
    {
        _allParams->display( std::cout );
    }

    NOMAD::OutputQueue::getInstance()->initParameters( _allParams->getDispParams() );
    NOMAD::OutputDirectToFile::getInstance()->init( _allParams->getDispParams() );

    // Create internal cache
    // The attribute passed to the function must be true for a cold suggest.
    // A cold suggest has a new MainStep created everytime suggest is called.
    // We rerun the algo using the points provided up to a certain state from a previous suggest. Then new points can be suggested.
    // For warm suggest the MainStep is reused and algo iterates from the state it was before.
    createCache(_allParams->getAttributeValue<bool>("USE_CACHE_FILE_FOR_RERUN"));

    size_t nbLHEval = _allParams->getAttributeValue<size_t>("LH_EVAL");

    if (0 != nbLHEval)
    {
        suggestedPoints = suggestFromLH(nbLHEval);
    }
    else if (_allParams->getAttributeValue<bool>("MEGA_SEARCH_POLL"))
    {
        auto cacheFile = _allParams->getCacheParams()->getAttributeValue<std::string>("CACHE_FILE");
        if (cacheFile.empty() && 0 == NOMAD::CacheBase::getInstance()->size())
        {
            std::string err = "Cache file is not provided or is empty. A cache is required to obtain Suggest points from a Mads MegaSearchPoll. To create a cache file, use suggest with LH_EVAL.";
            throw NOMAD::StepException(__FILE__,__LINE__, err, this);
        }
        auto lhSearchType = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
        if (0 != lhSearchType.getNbInitial())
        {
            std::string err = "LH_SEARCH's first value should be set to zero when calling Suggest with Mads MegaSearchPoll ";
            throw NOMAD::StepException(__FILE__,__LINE__, err, this);
        }

        // Update X0 from CacheFile - to satisfy MadsInitialization().
        updateX0sFromCacheAndFromLHSInit();
        auto x0s = _allParams->getPbParams()->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
        if (x0s.empty() || x0s[0].toBeDefined())
        {
            AddOutputInfo("No X0 is available. Cannot suggest any new point with MegaSearchPoll");
            return suggestedPoints;
        }

        // Need to start evaluator control.
        // Even though there are no BB evaluations done by suggest(),
        // the EvaluatorControl still holds information,
        // and is especially useful for Model searches.
        auto evaluatorControl = NOMAD::EvcInterface::getEvaluatorControl();
        if (nullptr == evaluatorControl)
        {
            if (_evaluators.empty())
            {
                // Batch mode. Create Evaluator on the go.
                // We need an evaluator to store eval params. This is used to get bb output type and other info.
                // A fake (EcalXDefined::UNDEFINED) bb evaluator is provided. This evaluator cannot evaluate a point. An exception is triggered when calling
                _evaluators.push_back(std::make_shared<NOMAD::Evaluator>(_allParams->getEvalParams(),
                                                                         NOMAD::EvalType::BB,
                                                                         NOMAD::EvalXDefined::UNDEFINED ));

            }
            else
            {
                    std::string err = "An evaluator has been set in the main step. Suggest does not require one. This can be problematic. ";
                    throw NOMAD::StepException(__FILE__,__LINE__, err, this);
            }

            std::unique_ptr<NOMAD::EvaluatorControlParameters> evalContParams(new NOMAD::EvaluatorControlParameters(*_allParams->getEvaluatorControlParams()));
            evalContParams->checkAndComply();
            evaluatorControl = std::make_shared<NOMAD::EvaluatorControl>(_allParams->getEvaluatorControlGlobalParams(),
                                                          std::move(evalContParams));

            // Add a fake evaluator
            evaluatorControl->addEvaluator(_evaluators[0]);


            NOMAD::EvcInterface::setEvaluatorControl(std::move(evaluatorControl));
        }

        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
        auto mads = std::make_shared<NOMAD::Mads>(this,
                                              stopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams(),
                                              false /* false: barrier not initialized from cache, use X0s*/);

        suggestedPoints = mads->suggest();
    }
    else
    {
            std::string err = "Suggest currently supports only Mads MEGA_SEARH_POLL or LH_EVAL. LH_EVAL should be used only when no cache is available. ";
            throw NOMAD::StepException(__FILE__,__LINE__, err, this);
    }

    AddOutputInfo("End step " + getName(), false, true);

    return suggestedPoints;
}


void NOMAD::MainStep::observe(const std::vector<NOMAD::EvalPoint>& evalPointList)
{
    AddOutputInfo("Start step " + getName() , true, false);

    // Display attributes and check attribute consistency
    if (_allParams->getAttributeValue<int>("DISPLAY_DEGREE") >= (int)NOMAD::OutputLevel::LEVEL_DEBUG)
    {
        _allParams->display(std::cout);
    }

    NOMAD::OutputQueue::getInstance()->initParameters(_allParams->getDispParams());
    NOMAD::OutputDirectToFile::getInstance()->init(_allParams->getDispParams());

    // The attribute passed to the function is always false (no cache file read for rerun). The points observed are passed as list to be put into the cache file
    createCache(false);
    if (evalPointList.size() > 0)
    {
        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
        auto mads = std::make_shared<NOMAD::Mads>(this,
                                              stopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams(),
                                                  false /* false: cache file not used*/);
        mads->observe(evalPointList);
        // Update interesting parameters
        _allParams->setAttributeValue("INITIAL_FRAME_SIZE", mads->getPbParams()->getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE"));
        _allParams->setAttributeValue("H_MAX_0", mads->getRunParams()->getAttributeValue<NOMAD::Double>("H_MAX_0"));
        _allParams->getPbParams()->doNotShowWarnings();
        _allParams->checkAndComply();
    }

    AddOutputInfo("End step " + getName(), false, true);
}


// This version to be used for the PyNomad interface
std::vector<std::string> NOMAD::MainStep::observe(const NOMAD::ArrayOfPoint& xs,
                                                  const vector<NOMAD::ArrayOfDouble>& fxs, const std::string & destinationCacheFileName )
{
    // Convert input xs and fxs to a vector of EvalPoints.
    std::vector<NOMAD::EvalPoint> evalPointList;
    if (xs.size() != fxs.size())
    {
        throw NOMAD::StepException(__FILE__,__LINE__,
                    "Observe: Input points and input values should have the same size.",
                    this);
    }

    auto bbOutputType = _allParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    for (size_t i = 0; i < xs.size(); i++)
    {
        NOMAD::EvalPoint evalPoint(xs[i]);
        evalPoint.setBBO(fxs[i].display(), bbOutputType, NOMAD::EvalType::BB);
        evalPointList.push_back(evalPoint);
    }
    observe(evalPointList);

    std::vector<std::string> updatedParams; // Return parameters as vector of strings, more convenient for PyNomad
    updatedParams.push_back("INITIAL_FRAME_SIZE ( " + _allParams->getPbParams()->getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE").display() + " )");
    updatedParams.push_back("H_MAX_0 " + _allParams->getRunParams()->getAttributeValue<NOMAD::Double>("H_MAX_0").display());

    // Save cache.
    // Note: This is usually done in Algorithm::end() via setInformationForHotRestart().
    // Since mads->end() is not called, we have to call it here.
    if(destinationCacheFileName.size() != 0)
        NOMAD::CacheBase::getInstance()->setFileName(destinationCacheFileName);

    if(NOMAD::CacheBase::getInstance()->getFileName().size() != 0)
        NOMAD::CacheBase::getInstance()->write();

    return updatedParams;
}


void NOMAD::MainStep::startImp()
{
#ifdef TIME_STATS
    _startTime = NOMAD::Clock::getCPUTime();
#endif // TIME_STATS


    if (nullptr == _allParams)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Using Library mode. Parameters must be set prior to running MainStep step.");
    }

    // Clear Subproblem map before a new run.
    // This is especially useful in the case of the Runner.
    NOMAD::SubproblemManager::getInstance()->reset();

    // reset PROBLEM_DIR attribute
    // (not obtained by reading paramFile):
    // ------------------------------------
    std::string problemDir;
    size_t k = _paramFileName.find_last_of ( NOMAD::DIR_SEP );
    // Reset only if param file is not in working dir
    if ( !_paramFileName.empty() && k < _paramFileName.size() )
    {
        problemDir = _paramFileName.substr (0,k) + NOMAD::DIR_SEP;
        _allParams->setAttributeValue("PROBLEM_DIR", problemDir);
    }

    // Batch mode / Read parameters file
    if (!_paramFileName.empty())
    {
        AddOutputInfo("Parameters file: " + _paramFileName);

        // Read all entries and set attribute values for each type of parameters
        _allParams->read( _paramFileName );
    }

    // Display attributes and check attribute consistency
    _allParams->checkAndComply();
    if (_allParams->getAttributeValue<int>("DISPLAY_DEGREE") >= (int)NOMAD::OutputLevel::LEVEL_DEBUG)
    {
        _allParams->display( std::cout );
    }

    NOMAD::OutputQueue::getInstance()->initParameters( _allParams->getDispParams() );
    NOMAD::OutputDirectToFile::getInstance()->init( _allParams->getDispParams() );

    // Solution file is written only at the end
    if (_allParams->getAttributeValue<bool>("SOLUTION_FILE_FINAL"))
    {
        // Disable solution file until the end
        NOMAD::OutputDirectToFile::getInstance()->disableSolutionFile();
    }

    createCache(_allParams->getAttributeValue<bool>("USE_CACHE_FILE_FOR_RERUN"));
    updateX0sFromCacheAndFromLHSInit();
    auto x0s = _allParams->getPbParams()->getAttributeValue<NOMAD::ArrayOfPoint>("X0");

    auto nbLHEval = _allParams->getRunParams()->getAttributeValue<size_t>("LH_EVAL");
    auto doRandomAlgo = _allParams->getRunParams()->getAttributeValue<bool>("RANDOM_ALGO_OPTIMIZATION");
    auto doDMultiMadsAlgo = _allParams->getRunParams()->getAttributeValue<bool>("DMULTIMADS_OPTIMIZATION");


    if ( (x0s.empty() || x0s[0].toBeDefined()) && nbLHEval == 0 && !doRandomAlgo && !doDMultiMadsAlgo)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Need X0 to continue.");
    }


    // Setup EvaluatorControl
    setNumThreads();

    if (_evaluators.empty())
    {
        // Down the road: manage multi fidelity surrogates. Add evaluators for each fidelity or pass the fidelity level to the constructor.

        // Batch mode. Create Evaluator on the go.
        bool surrogateAsBB = _allParams->getAttributeValue<bool>("EVAL_SURROGATE_OPTIMIZATION");
        auto evalType = (surrogateAsBB) ? NOMAD::EvalType::SURROGATE : NOMAD::EvalType::BB;
        _evaluators.push_back(std::make_shared<NOMAD::Evaluator>(_allParams->getEvalParams(),
                                                          evalType,
                                                          NOMAD::EvalXDefined::USE_BB_EVAL));
        if (!surrogateAsBB)
        {
            auto evalSortType = _allParams->getAttributeValue<NOMAD::EvalSortType>("EVAL_QUEUE_SORT");
            bool surrogateForVNS = _allParams->getAttributeValue<bool>("VNS_MADS_SEARCH") && _allParams->getAttributeValue<bool>("VNS_MADS_SEARCH_WITH_SURROGATE");
            if (NOMAD::EvalSortType::SURROGATE == evalSortType || surrogateForVNS)
            {
                _evaluators.push_back(std::make_shared<NOMAD::Evaluator>(_allParams->getEvalParams(),
                                                                         NOMAD::EvalType::SURROGATE,
                                                                         NOMAD::EvalXDefined::USE_BB_EVAL));

            }
        }

    }

    auto evaluatorControl = NOMAD::EvcInterface::getEvaluatorControl();
    if (nullptr == evaluatorControl)
    {
        std::unique_ptr<NOMAD::EvaluatorControlParameters> evalContParams(new NOMAD::EvaluatorControlParameters(*_allParams->getEvaluatorControlParams()));
        evalContParams->checkAndComply();
        evaluatorControl = std::make_shared<NOMAD::EvaluatorControl>(_allParams->getEvaluatorControlGlobalParams(),
                                                      std::move(evalContParams));
        for( const auto & evaluator : _evaluators)
        {
            evaluatorControl->addEvaluator(evaluator);
        }
        // Need to force selection of BB evaluator. Otherwise the last added evaluator is used.
        evaluatorControl->setCurrentEvaluatorType(NOMAD::EvalType::BB);

        NOMAD::EvcInterface::setEvaluatorControl(std::move(evaluatorControl));
    }

    // Currently this does nothing.
    NOMAD::EvcInterface::getEvaluatorControl()->start();

    // Create the Algorithms that we want to solve.
    // Currently available: PhaseOne, Mads, LH, NM, SgtelibModelEval, CS
    // QuadModelOptimization
    // Note: These algorithm parameters are mutually exclusive:
    // LH_EVAL NM_OPTIMIZATION PSD_MADS_OPTIMIZATION QUAD_MODEL_OPTIMIZATION
    // SGTELIB_MODEL_EVAL SSD_MADS_OPTIMIZATION.
    // This is caught by checkAndComply().
    _algos.clear();

    auto doCSoptimization = _allParams->getRunParams()->getAttributeValue<bool>("CS_OPTIMIZATION");
    auto doNMOptimization = _allParams->getRunParams()->getAttributeValue<bool>("NM_OPTIMIZATION");
#ifdef _OPENMP
    bool doPSDMads = _allParams->getRunParams()->getAttributeValue<bool>("PSD_MADS_OPTIMIZATION");
#endif
#ifdef USE_SGTELIB
    bool doQuadModelOpt = _allParams->getRunParams()->getAttributeValue<bool>("QUAD_MODEL_OPTIMIZATION");
    bool doSgtelibModelEval = _allParams->getRunParams()->getAttributeValue<bool>("SGTELIB_MODEL_EVAL");
    bool doQPSolverQuadModelOpt = _allParams->getRunParams()->getAttributeValue<bool>("QP_OPTIMIZATION");
#endif
    bool doSSDMads = _allParams->getRunParams()->getAttributeValue<bool>("SSD_MADS_OPTIMIZATION");
    bool doRandomAlgoOpt = _allParams->getRunParams()->getAttributeValue<bool>("RANDOM_ALGO_OPTIMIZATION");
    bool doDMultimadsOpt = _allParams->getRunParams()->getAttributeValue<bool>("DMULTIMADS_OPTIMIZATION");
    bool doDiscoMads = _allParams->getRunParams()->getAttributeValue<bool>("DISCO_MADS_OPTIMIZATION");

    bool evalOpportunistic = _allParams->getAttributeValue<bool>("EVAL_OPPORTUNISTIC");
    // LH_EVAL can be done before another algo (not after!)
    if ( nbLHEval > 0 )
    {
        auto lhStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::LHStopType>>();
        auto lh = std::make_shared<NOMAD::LH>(this,
                                              lhStopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams());

        // All the LH sample points must be evaluated. No opportunism.
        if ( evalOpportunistic )
            AddOutputInfo("Opportunistic evaluation is disabled for LH when ran as a standalone algorithm.");
        lh->setEvalOpportunistic(false);

        _algos.push_back(lh);
    }


    if ( doNMOptimization )
    {
        auto nmStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::NMStopType>>();

        NOMAD::EvcInterface::getEvaluatorControl()->setOpportunisticEval(false);

        auto nm = std::make_shared<NOMAD::NM>(this,
                                              nmStopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams());
        // All the NM points must be evaluated. No opportunism.
        if (evalOpportunistic)
        {
            AddOutputInfo("Opportunistic evaluation is disabled for NM when ran as a single algorithm.");
        }
        nm->setEvalOpportunistic(false);

        _algos.push_back(nm);
    }
    else if( doCSoptimization )
    {
        auto csStopReason = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::CSStopType>>();

        auto cs = std::make_shared<NOMAD::CS>(this,
                                              csStopReason ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams());

        cs->setEvalOpportunistic(evalOpportunistic);

        _algos.push_back(cs);
    }

#ifdef _OPENMP
    else if (doPSDMads)
    {
        auto PSDMadsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>(); // A PSD-MADS has a MADS type stop reason.
        auto psd = std::make_shared<NOMAD::PSDMads>(this,
                                                    _evaluators,
                                                    _allParams->getEvaluatorControlParams(),
                                                    PSDMadsStopReasons,
                                                    _allParams->getRunParams(),
                                                    _allParams->getPbParams());
        psd->setEvalOpportunistic(evalOpportunistic);

        _algos.push_back(psd);
    }
#endif
#ifdef USE_SGTELIB
    else if (doQuadModelOpt)
    {
        auto quadModelStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::ModelStopType>>();

        // Deactivate mega search poll
        _allParams->setAttributeValue("MEGA_SEARCH_POLL", false);
        _allParams->checkAndComply();

        auto quadModelAlgo = std::make_shared<NOMAD::QuadModelAlgo>(this,
                                                        quadModelStopReasons,
                                                        _allParams->getRunParams(),
                                                        _allParams->getPbParams());

        // All the Sgtelib Model sample points are evaluated sequentially. No opportunism.
        quadModelAlgo->setEvalOpportunistic(false);

        _algos.push_back(quadModelAlgo);
    }
    else if (doSgtelibModelEval)
    {
        auto sgtelibModelStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::ModelStopType>>();

        std::shared_ptr<NOMAD::BarrierBase> barrier = nullptr;
        if (NOMAD::CacheBase::getInstance()->size() > 0)
        {
            // Create a single objective progressive barrier
            barrier = std::make_shared<NOMAD::ProgressiveBarrier>();
        }
        auto sgtelibModel = std::make_shared<NOMAD::SgtelibModel>(this,
                                                        sgtelibModelStopReasons,
                                                        barrier,
                                                        _allParams->getRunParams(),
                                                        _allParams->getPbParams(),
                                                        nullptr);   // no mesh
        // All the Sgtelib Model sample points are evaluated sequentially. No opportunism.
        sgtelibModel->setEvalOpportunistic(false);

        _algos.push_back(sgtelibModel);
    }
    else if (doQPSolverQuadModelOpt)
    {
        auto quadModelStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::ModelStopType>>();

        // All the Sgtelib Model sample points are evaluated sequentially. No opportunism.
        NOMAD::EvcInterface::getEvaluatorControl()->setOpportunisticEval(false);
        _allParams->setAttributeValue("MEGA_SEARCH_POLL", false);
        _allParams->checkAndComply();

        auto qpSolverQuadModelOptAlgo = std::make_shared<NOMAD::QPSolverAlgo>(this,
                                                        quadModelStopReasons,
                                                        _allParams->getRunParams(),
                                                        _allParams->getPbParams());
        _algos.push_back(qpSolverQuadModelOptAlgo);
    }
#endif // USE_SGTELIB
    else if (doSSDMads)
    {
        auto SSDMadsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>(); // A SSD-MADS has MADS Stop type.
        auto ssd = std::make_shared<NOMAD::SSDMads>(this,
                                                    SSDMadsStopReasons,
                                                    _allParams->getRunParams(),
                                                    _allParams->getPbParams());
        ssd->setEvalOpportunistic(evalOpportunistic);
        _algos.push_back(ssd);
    }
    else if (doRandomAlgoOpt)
    {
        auto randomAlgoStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::RandomAlgoStopType>>(); // A template algo has a special Stop type.
        auto randomAlgo = std::make_shared<NOMAD::TemplateAlgo>(this,
                                                                  randomAlgoStopReasons,
                                                                  _allParams->getRunParams(),
                                                                  _allParams->getPbParams());
        randomAlgo->setEvalOpportunistic(evalOpportunistic);

        _algos.push_back(randomAlgo);
    }
    else if (doDMultimadsOpt)
    {
        auto dMultiMadsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>(); // DMultiMads  Stop type.
        auto dMultiMadsAlgo = std::make_shared<NOMAD::DMultiMads>(this,
                                                                  dMultiMadsStopReasons,
                                                                  _allParams->getRunParams(),
                                                                  _allParams->getPbParams());
        dMultiMadsAlgo->setEvalOpportunistic(evalOpportunistic);

        _algos.push_back(dMultiMadsAlgo);
    }
    else if (doDiscoMads)
    {
        auto discoMadsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
        auto discoMadsAlgo = std::make_shared<NOMAD::DiscoMads>(this,
                                                                  discoMadsStopReasons,
                                                                  _allParams->getRunParams(),
                                                                  _allParams->getPbParams());
        _algos.push_back(discoMadsAlgo);
    }
    else
    {
        // Mads is the default algorithm (PhaseOne search is managed in MadsInitialization)

        if (nbLHEval > 0)
        {
            std::cout << "Warning: LH_EVAL is performed but Mads is disabled. To perform LH initialization for Mads use LH_SEARCH."<<std::endl ;
            return;

        }

        // The stop reasons for mads
        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

        // Default behavior: create Mads and add it to _algos.
        auto mads = std::make_shared<NOMAD::Mads>(this,
                                                  stopReasons ,
                                                  _allParams->getRunParams(),
                                                  _allParams->getPbParams(),
                                                  false /* false: do no initialize barrier from cache, X0s is used instead */);
        mads->setEvalOpportunistic(evalOpportunistic);
        _algos.push_back(mads);
    }

    bool useIbex = _allParams->getRunParams()->getAttributeValue<bool>("USE_IBEX");
	if (useIbex)
	{

#ifdef USE_IBEX
		// A Set file determine the feasible domain

		// If the Set file is already created, we can load it
		bool setFile = _allParams->getRunParams()->getAttributeValue<bool>("SET_FILE");
		if (setFile)
		{
			string setFileName = _allParams->getAttributeValue<string>("SET_FILE_NAME");
			const char * c = setFileName.c_str();
			_set = std::make_shared<ibex::Set>(c);
		}

		// Else, we can load the system file of the problem (containing the variables, constraints, ...) and then create the Set
		else
		{
			string constraintsFileName = _allParams->getAttributeValue<string>("SYSTEM_FILE_NAME");
		        const char * c = constraintsFileName.c_str();

			ibex::System sys(c);
			ibex::IntervalVector box = sys.box;

			const int n = sys.nb_ctr;
			ibex::Array<ibex::NumConstraint> constraints = sys.ctrs;
			ibex::Array<ibex::Sep> separators(n);

			for (int i = 0; i<n; i++)
			{
				separators.set_ref(i, *new ibex::SepFwdBwd(constraints[i]));
			}

			ibex::SepInter sep_poly_in(separators);
			_set = std::make_shared<ibex::Set>(box);
			sep_poly_in.contract(*_set, 0.1);
		}
#else
    throw NOMAD::Exception(__FILE__, __LINE__, "IBEX projection requires to configure and build NOMAD with the option USE_IBEX.");
#endif
    }

}


bool NOMAD::MainStep::runImp()
{
    bool ret = false;
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    for (auto algo : _algos)
    {
        // Get ready to evaluate new points for each new algo.
        evc->restart();

        // Begin parallel region.
        // Note: always use default(none) and never default(shared), which is much
        // too risky.
#ifdef _OPENMP
#pragma omp parallel default(none) shared(algo,evc,ret)
#endif // _OPENMP
        {

            // Set opportunism for main thread only
            if (evc->isMainThread(NOMAD::getThreadNum()))
            {
                evc->setOpportunisticEval(algo->getEvalOpportunistic());
            }

            // Start evaluatorControl on all threads
            evc->run();

            if (evc->isMainThread(NOMAD::getThreadNum()))
            {

                // Note: Algo start has been moved from outside to inside the parallel region
                // so that evaluation of multiple X0s is done in parallel.
                // This may cause trouble. Still under investigation.
                algo->start();

                // Algo run is done in main thread(s) only.
                ret = algo->run();

                // When algo is done, evaluatorControl will wait for all main threads to be done.
                evc->stop();
            }
        }   // End of parallel region.
        algo->end();

     }

    return ret;
}


void NOMAD::MainStep::endImp()
{
    // Pass last algo stop reason to MainStep (before clearing the algos)
    _stopReasons = _algos.back()->getAllStopReasons();


#ifdef TIME_STATS
   _totalCPUTime += NOMAD::Clock::getCPUTime() - _startTime;
#endif // TIME_STATS

    displayDetailedStats();

    writeFinalSolutionFile();

    _algos.clear();

}


int NOMAD::MainStep::getNumThreads() const
{
    int nbThreadsParam = 1;
#ifdef _OPENMP
    // Number of threads on which to do the evaluation.
    nbThreadsParam = _allParams->getAttributeValue<int>("NB_THREADS_OPENMP");
    int nbThreadsHard = static_cast<int>(std::thread::hardware_concurrency());
    if (nbThreadsParam < 1)
    {
        nbThreadsParam = nbThreadsHard;
    }
    else if (nbThreadsParam > nbThreadsHard)
    {
        std::string s = "NB_THREADS_OPENMP exceeds the number of threads registered for this hardware: ";
        s += NOMAD::itos(nbThreadsHard);
        s += ". If this is true, it is not efficient. Let's continue anyway.";
        NOMAD::OutputQueue::Add(s,NOMAD::OutputLevel::LEVEL_NORMAL);
    }
#endif // _OPENMP

    return nbThreadsParam;
}


void NOMAD::MainStep::setNumThreads() const
{
#ifdef _OPENMP
    // Set number of threads on which to do the evaluation.
    int nbThreadsParam = getNumThreads();
    omp_set_num_threads(nbThreadsParam);

    std::string s = "OpenMP will use " + NOMAD::itos(nbThreadsParam) + " thread";
    if (nbThreadsParam > 1) { s += "s"; }
    s += " in the parallel regions.\n";
    NOMAD::OutputQueue::Add(s,NOMAD::OutputLevel::LEVEL_NORMAL);

#endif // _OPENMP
}


void NOMAD::MainStep::printNumThreads() const
{
#ifdef _OPENMP
    // Print the actual number of threads used (once if in parallel region).
#pragma omp single nowait
    {
        int nbThreads = omp_get_num_threads();
        std::string s = "Number of threads currently in use in the region: " + NOMAD::itos(nbThreads) ;
        NOMAD::OutputQueue::Add(s,NOMAD::OutputLevel::LEVEL_NORMAL);
    }
#endif // _OPENMP
}


// Detect that we must do a Phase One.
bool NOMAD::MainStep::detectPhaseOne()
{
    bool hasEBConstraints = false;
    bool hasNoFeas = !NOMAD::CacheBase::getInstance()->hasFeas();

    auto bbOutputTypeList = _allParams->getEvalParams()->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    if (std::find(bbOutputTypeList.begin(), bbOutputTypeList.end(), NOMAD::BBOutputType("EB"))
        != bbOutputTypeList.end())
    {
        hasEBConstraints = true;
    }

    return hasEBConstraints && hasNoFeas;
}


// Helper for start
void NOMAD::MainStep::createCache(bool useCacheForRerun) const
{
    // Creation of an instance of CacheSet with CacheParameters
    // This must be done ONCE before accessing the singleton using NOMAD::CacheBase::getInstance()
    try
    {
        NOMAD::CacheBase::getInstance();
    }
    catch (...)
    {
        NOMAD::CacheSet::setInstance(_allParams->getCacheParams(),
                                     _allParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"),
                                     _allParams->getAttributeValue<NOMAD::ArrayOfDouble>("BB_EVAL_FORMAT"));
        if (useCacheForRerun)
        {
            // Swap cache and cacheForRerun
            NOMAD::CacheBase::getInstance()->moveEvalPointToCacheForRerun();
        }

    }
}


void NOMAD::MainStep::updateX0sFromCacheAndFromLHSInit() const
{
    // Update X0s, if needed.
    auto x0s = _allParams->getPbParams()->getAttributeValue<NOMAD::ArrayOfPoint>("X0");

    bool updatedX0s = false;

    if (x0s.empty() || x0s[0].toBeDefined())
    {
        // X0 not provided directly by user. Find them in cache.
        x0s.clear();

        bool canUseCache = (NOMAD::CacheBase::getInstance()->size() > 0);

        if (canUseCache)
        {
            // Use best points in cache for x0.
            std::vector<NOMAD::EvalPoint> evalPointList;
            // Note: We are working in full dimension here, not in subproblem.
            // For this reason, use cache instance directly, not CacheInterface.
            auto fixedVariable = _allParams->getPbParams()->getAttributeValue<NOMAD::Point>("FIXED_VARIABLE");
            auto evalType = NOMAD::EvalType::BB;
            auto computeType = NOMAD::ComputeType::STANDARD;
            NOMAD::CacheBase::getInstance()->findBestFeas(evalPointList,
                                                          fixedVariable, evalType,
                                                          computeType);
            if (0 == evalPointList.size())
            {
                auto hMax = _allParams->getRunParams()->getAttributeValue<NOMAD::Double>("H_MAX_0");
                NOMAD::CacheBase::getInstance()->findBestInf(evalPointList,
                                                             hMax, fixedVariable,
                                                             evalType, computeType);
            }
            updatedX0s = (evalPointList.size() > 0)? true:false;
            for (const auto & evalPoint: evalPointList)
            {
                x0s.push_back(evalPoint);
            }
        }
    }

    // Complete X0 with LHS sampled points
    auto lhSearchType = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
    auto fixedVariables = _pbParams->getAttributeValue<NOMAD::Point>("FIXED_VARIABLE");
    bool canUseLH = (lhSearchType.isEnabled() && lhSearchType.getNbInitial() > 0);
    if (canUseLH)
    {
        NOMAD::ArrayOfPoint sampleEvalPoints = suggestFromLH(lhSearchType.getNbInitial());
        for (size_t i = 0; i < sampleEvalPoints.size(); i++)
        {
            if (fixedVariables.nbDefined() >= 1)
            {
                sampleEvalPoints[i] = sampleEvalPoints[i].makeFullSpacePointFromFixed(fixedVariables);
            }
            x0s.push_back(sampleEvalPoints[i]);
        }
        updatedX0s = updatedX0s || (sampleEvalPoints.size() > 0)? true:false;
    }
    if (updatedX0s)
    {
        _allParams->getPbParams()->setAttributeValue("X0", x0s);
        _allParams->getPbParams()->checkAndComply();
    }

}


NOMAD::ArrayOfPoint NOMAD::MainStep::suggestFromLH(const size_t nbPoints) const
{
    auto lhStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::LHStopType>>();
    auto lhRunParams = std::make_shared<NOMAD::RunParameters>(*_allParams->getRunParams());
    lhRunParams->setAttributeValue("LH_EVAL", nbPoints);

    lhRunParams->checkAndComply(_allParams->getEvaluatorControlGlobalParams(), _allParams->getPbParams());

    NOMAD::LH lh(this, lhStopReasons, lhRunParams, _allParams->getPbParams());

    return lh.suggest();
}


/*------------------------------------------*/
/*             display NOMAD usage          */
/*------------------------------------------*/
void NOMAD::MainStep::displayUsage( const char* exeName )
{
    // Get executable name without path
    std::string strExeName(exeName);
    std::size_t i = strExeName.rfind("\\");
    if (i == std::string::npos)
    {
        i = strExeName.rfind("/");
    }
    if (i != std::string::npos)
    {
        strExeName.replace(0, i+1, "");
    }

    std::string usage;
    usage += \
        "Run NOMAD      : " + strExeName + " parameters_file\n" \
      + "Info           : " + strExeName + " -i\n" \
      + "Help           : " + strExeName + " -h [keyword]\n" \
      + "Version        : " + strExeName + " -v\n" \
      + "Usage          : " + strExeName + " -u\n\n";

    NOMAD::OutputQueue::Add(usage, NOMAD::OutputLevel::LEVEL_ERROR);

}


/*------------------------------------------*/
/*          display NOMAD version           */
/*------------------------------------------*/
void NOMAD::MainStep::displayVersion()
{
    std::string version = "Version ";
    version += NOMAD_VERSION_NUMBER;
#ifdef DEBUG
    version += " Debug.";
#else
    version += " Release.";
#endif
#ifdef _OPENMP
    version += " Using OpenMP.";
#else
    version += " Not using OpenMP.";
#endif // _OPENMP

#ifdef USE_SGTELIB
    version += " Using SGTELIB.";
#else
    version += " Not using SGTELIB.";
#endif

    NOMAD::OutputQueue::Add(version, NOMAD::OutputLevel::LEVEL_VERY_HIGH);

}


/*------------------------------------------*/
/*             display NOMAD info           */
/*------------------------------------------*/
void NOMAD::MainStep::displayInfo()
{

    std::string info;

    info += "Version " ;
    info += NOMAD_VERSION_NUMBER ;
    
#ifdef _OPENMP
    info += ". Using OpenMP.";
#else
    info += ". Not using OpenMP.";
#endif // _OPENMP
    info += " \n \n";

    info += "NOMAD 4 has been created by \n";
    info += "    Viviane Rochon Montplaisir \n";
    info += "    Christophe Tribes \n\n";

    info += "The copyright of NOMAD 4 is owned by \n";
    info += "    Charles Audet \n";
    info += "    Sebastien Le Digabel \n";
    info += "    Viviane Rochon Montplaisir \n";
    info += "    Christophe Tribes \n\n";

    info += "NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada, \n";
    info += " NSERC (Natural Sciences and Engineering Research Council of Canada), \n";
    info += " InnovÉÉ (Innovation en Énergie Électrique) and \n";
    info += " IVADO (The Institute for Data Valorization) \n\n";

    info += "Download  : https://www.gerad.ca/nomad or \n";
    info += " https://github.com/bbopt/nomad \n" ;
    info += "License   : see LICENSE file \n";
    info += "User guide: https://nomad-4-user-guide.readthedocs.io \n";
    info += "Help      : run nomad -h KEYWORD on the command line \n";
    info += "Examples  : see \'examples\' directory \n\n";

    info += "Please report bugs to nomad@gerad.ca or \n";
    info += "create an issue at https://github.com/bbopt/nomad \n\n";


    NOMAD::OutputQueue::Add(info, NOMAD::OutputLevel::LEVEL_NORMAL);
}


/*------------------------------------------------------*/
/*             display NOMAD help on a subject          */
/*------------------------------------------------------*/
void NOMAD::MainStep::displayHelp( const std::string & helpSubject , bool devHelp )
{
    _allParams->displayHelp( helpSubject, devHelp, std::cout );
}

/*------------------------------------------------------*/
/*     display all parameters in CSV format for doc     */
/*------------------------------------------------------*/
void NOMAD::MainStep::displayCSVDoc()
{
    _allParams->displayCSVDoc( std::cout );
}

// What to do when user interrupts NOMAD
void NOMAD::MainStep::hotRestartOnUserInterrupt()
{
    // Update AllParams ...

    // User asked for interruption.
    // Read new parameters and then continue.

    // Read parameter file, or read new parameters inline.

    // Note that at this point it is not verified if the new parameters
    // are allowed to change. For instance, changing Problem parameters
    // might cause issues.
    // Changing GRANULARITY for a lower value might not affect
    // MIN_MESH_SIZE that was previously set (although we try to update it). Etc.

    hotRestartBeginHelper();

    if (!getUserTerminate())
    {
        std::cout << "Hot restart" ;

        // Do not use a shared_ptr _evaluator because it is NULL in this function
        std::vector<std::string> paramLines;
        _cbHotRestart(paramLines);

        if ( paramLines.size() == 0 )
        {
            std::cout << std::endl << "Enter a parameter file name," << std::endl;
            std::cout << "or enter parameter values, ending with CTRL-D." << std::endl;

            std::string line;
            std::getline(std::cin, line);
            // Is this line a parameter file?
            if (NOMAD::checkReadFile(line))
            {
                std::cout << "Reading parameter file: " << line << std::endl;
                _allParams->read(line, true /* overwrite */);
            }
            else
            {
                _allParams->readParamLine(line);
                // Continue reading parameters
                while(!getUserTerminate() && std::getline(std::cin, line))
                {
                    // Read parameter
                    _allParams->readParamLine(line);
                }
            }
        }
        else
        {
            std::cout << ": read parameters update" << std::endl;
            for ( auto line : paramLines )
            {
                _allParams->readParamLine( line );
            }
        }

        _allParams->checkAndComply();
    }

    hotRestartEndHelper();
}


void NOMAD::MainStep::resetComponentsBetweenOptimization()
{
    // Make sure to clear the cache before the next run
    resetCache();

	// Reset user terminate. Important for Matlab and Python interface
	// with persistent memory of static variables
	NOMAD::Step::resetUserTerminate();
    NOMAD::Step::resetUserInterrupt();

    // Reset static tag counter
    NOMAD::EvalPoint::resetCurrentTag();

    // Reset SubproblemManager map
    NOMAD::SubproblemManager::getInstance()->reset();

    // Reset evaluator control
    NOMAD::EvcInterface::resetEvaluatorControl();

    // Reset the random number generator to its initial state
    NOMAD::RNG::reset();

    // Reset parameter entries
    NOMAD::Parameters::eraseAllEntries();

}

void NOMAD::MainStep::resetCache()
{
    // Get a new cache
    NOMAD::CacheBase::resetInstance(); // Need to reset the singleton. When calling createCache there is no instance and we are sure to call NOMAD::CacheSet::setInstance from scratch. The cache file is read and the cache is set with its content.
}

void NOMAD::MainStep::resetEvaluatorControl()
{
    NOMAD::EvcInterface::resetEvaluatorControl();
}


void NOMAD::MainStep::displayDetailedStats() const
{
    // Display detailed stats
    std::string evalStatsFile = _allParams->getAttributeValue<std::string>("EVAL_STATS_FILE");

    if (evalStatsFile.empty() || evalStatsFile.compare("-") == 0)
    {
        return;
    }

    NOMAD::ArrayOfString s1,s2;

#ifndef TIME_STATS
    s1.add("Total real time (s):");
    s2.add(std::to_string(NOMAD::Clock::getTimeSinceStart()));
#endif

    s1.add("Blackbox evaluations:");
    size_t bbEval = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
    s2.add(NOMAD::itos(bbEval));

    s1.add("Blackbox evaluations (from cache file for rerun):");
    size_t bbEvalFromCacheForRerun = NOMAD::EvcInterface::getEvaluatorControl()->getBbEvalFromCacheForRerun();
    s2.add(NOMAD::itos(bbEvalFromCacheForRerun));

    s1.add("Cache hits:");
    size_t nbCacheHits = NOMAD::CacheBase::getNbCacheHits();
    s2.add(NOMAD::itos(nbCacheHits));

    s1.add("Total number of evaluations:");
    size_t nbEval = NOMAD::EvcInterface::getEvaluatorControl()->getNbEval();
    s2.add(NOMAD::itos(nbEval));

    s1.add("Blackbox evaluations that are not counted:");
    int nbEvalNoCount = static_cast<int>(nbEval - bbEval - nbCacheHits);
    s2.add(NOMAD::itos(nbEvalNoCount));

    s1.add("Blackbox evaluations that are not ok:");
    size_t nbEvalNotOk = NOMAD::EvcInterface::getEvaluatorControl()->getBbEvalNotOk();
    s2.add(NOMAD::itos(nbEvalNotOk));

    s1.add("Blackbox feasible evaluations:");
    size_t nbFeasBBEval = NOMAD::EvcInterface::getEvaluatorControl()->getFeasBbEval();
    s2.add(NOMAD::itos(nbFeasBBEval));

    s1.add("Blackbox feasible success evaluations:");
    size_t nbRelSucc = NOMAD::EvcInterface::getEvaluatorControl()->getNbRelativeSuccess();
    s2.add(NOMAD::itos(nbRelSucc));

    s1.add("Block evaluations:");
    size_t blkEval = NOMAD::EvcInterface::getEvaluatorControl()->getBlockEval();
    s2.add(NOMAD::itos(blkEval));

    s1.add("Total surrogate evaluations:");
    size_t totalSurrogateEval = NOMAD::EvcInterface::getEvaluatorControl()->getSurrogateEval();
    s2.add(NOMAD::itos(totalSurrogateEval));

    s1.add("Total surrogate evaluations from cache file (for rerun):");
    size_t totalSurrogateEvalFromCacheForRerun = NOMAD::EvcInterface::getEvaluatorControl()->getSurrogateEvalFromCacheForRerun();
    s2.add(NOMAD::itos(totalSurrogateEvalFromCacheForRerun));

    s1.add("Total model evaluations:");
    size_t totalModelEval = NOMAD::EvcInterface::getEvaluatorControl()->getTotalModelEval();
    s2.add(NOMAD::itos(totalModelEval));

    s1.add("PhaseOne success evaluations:");
    size_t nbPhaseOneSucc = NOMAD::EvcInterface::getEvaluatorControl()->getNbPhaseOneSuccess();
    s2.add(NOMAD::itos(nbPhaseOneSucc));

    s1.add("Index of the best feasible evaluation:");
    size_t indexFeasEval = NOMAD::EvcInterface::getEvaluatorControl()->getIndexFeasEval();
    s2.add(NOMAD::itos(indexFeasEval));

    s1.add("Index of the best infeasible evaluation:");
    size_t indexInfeasEval = NOMAD::EvcInterface::getEvaluatorControl()->getIndexInfeasEval();
    s2.add(NOMAD::itos(indexInfeasEval));

    s1.add("Index of evaluation block containing the best feasible solution:");
    size_t indexSuccBlock = NOMAD::EvcInterface::getEvaluatorControl()->getIndexSuccBlockEval();
    s2.add(NOMAD::itos(indexSuccBlock));

    size_t nbRevealingIter = NOMAD::EvcInterface::getEvaluatorControl()-> getNbRevealingIter(); // Curr. only used in DiscoMads
    if (nbRevealingIter>0)
    {
        s1.add("Revealing Iterations:");
        s2.add(NOMAD::itos(nbRevealingIter));
    }

#ifdef TIME_STATS
    s1.add("");
    s2.add("");

    s1.add("Time statistics (s)");
    s2.add("");

    s1.add("Total real time:");
    s2.add(std::to_string(NOMAD::Clock::getTimeSinceStart()));

    s1.add("Total CPU time:");
    s2.add(std::to_string(_totalCPUTime));

    s1.add("Total time in Evaluator:");
    s2.add(std::to_string(NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime()));

    s1.add("Total time in Mads Iterations:");
    s2.add(std::to_string(NOMAD::MadsIteration::getIterTime()));

    s1.add("Total time in Search:");
    s2.add(std::to_string(NOMAD::MadsIteration::getSearchTime())
        +" (Eval: "
        + std::to_string(NOMAD::MadsIteration::getSearchEvalTime())
        +")");

    s1.add("Total time in Poll:");
    s2.add(std::to_string(NOMAD::MadsIteration::getPollTime())
          +" (Eval: "
          + std::to_string(NOMAD::MadsIteration::getPollEvalTime())
          +")");

    std::vector<double> searchTime = NOMAD::Search::getSearchTime();
    std::vector<double> searchEvalTime = NOMAD::Search::getSearchEvalTime();
    s1.add("Total time in Speculative Search:");
    s2.add(std::to_string(searchTime[0])
          + " (Eval: "
          + std::to_string(searchEvalTime[0])
          + ")");

    s1.add("Total time in User Search:");
    s2.add(std::to_string(searchTime[1])
          + " (Eval: "
          + std::to_string(searchEvalTime[1])
          + ")");

    s1.add("Total time in Sgtelib Search:");
    s2.add(std::to_string(searchTime[2])
          + " (Eval: "
          + std::to_string(searchEvalTime[2])
          + ")");

    s1.add("Total time in LH Search:");
    s2.add(std::to_string(searchTime[3])
          + " (Eval: "
          + std::to_string(searchEvalTime[3])
          + ")");

    s1.add("Total time in NM Search:");
    s2.add(std::to_string(searchTime[4])
          + " (Eval: "
          + std::to_string(searchEvalTime[4])
          + ")");

#endif // TIME_STATS

    NOMAD::ArrayOfString paddedStats = NOMAD::ArrayOfString::combineAndAddPadding(s1,s2);
    std::ofstream evalStatsStream;
    // Open eval stats file and clear it (trunc)
    evalStatsStream.open(evalStatsFile.c_str(), std::ofstream::out | std::ios::trunc);
    if (evalStatsStream.fail())
    {
        std::cout << "Warning: could not open evaluation stats file " << evalStatsFile << std::endl;
    }
    evalStatsStream << paddedStats << std::endl;
    evalStatsStream.close();
}

void NOMAD::MainStep::writeFinalSolutionFile() const
{
    if (_allParams->getAttributeValue<bool>("SOLUTION_FILE_FINAL"))
    {
        // Enable solution file that was disabled at start
        NOMAD::OutputDirectToFile::getInstance()->enableSolutionFile();

        auto barrier = _algos.back()->getMegaIterationBarrier();

        if (nullptr != barrier)
        {
            const std::vector<EvalPointPtr>& xFeas = barrier->getAllXFeas();
            if (xFeas.size() > 1)
            {
                // If we have a success and we have muliple best feasible solution, we rewrite the solution file.

                bool append = false;
                for (const EvalPointPtr & ev: xFeas)
                {
                    NOMAD::StatsInfo info;

                    info.setBBO(ev->getBBO(NOMAD::EvalType::BB));
                    info.setSol(*(ev->getX()));

                    NOMAD::OutputDirectToFile::Write(info, true, false /* do no write in history file */, append /* append in solution file */);
                    append = true;
                }
            }
        }
    }
}


void NOMAD::MainStep::addEvaluator(const EvaluatorPtr ev)
{
    NOMAD::EvalType evalTypeAdded = ev->getEvalType();

    if (NOMAD::EvalType::MODEL == evalTypeAdded)
    {
        std::string s = "Error in evaluator management: cannot add evaluator Model eval type in main step";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    auto it = std::find_if(_evaluators.begin(),_evaluators.end(), [evalTypeAdded](NOMAD::EvaluatorPtr e){ return e->getEvalType() == evalTypeAdded; });

    if ( _evaluators.end() != it )
    {
        std::string s = "Error in evaluator management: evaluator with EvalType = " + NOMAD::evalTypeToString(evalTypeAdded);
        s += " has already been added.";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    if (NOMAD::EvalType::SURROGATE == evalTypeAdded &&
        NOMAD::EvalSortType::SURROGATE != _allParams->getAttributeValue<NOMAD::EvalSortType>("EVAL_QUEUE_SORT") )
    {
        if ( ! _allParams->getAttributeValue<bool>("VNS_MADS_SEARCH") ||
            ( _allParams->getAttributeValue<bool>("VNS_MADS_SEARCH") &&
             ! _allParams->getAttributeValue<bool>("VNS_MADS_SEARCH_WITH_SURROGATE") ) )
        {
            std::cout << "Warning: A SURROGATE evaluator is available but it will not be used. To use it, set EVAL_QUEUE_SORT to SURROGATE or set VNS_MADS_SEARCH_WITH_SURROGATE." << std::endl;
        }

    }
    _evaluators.push_back(ev);

}

void NOMAD::MainStep::setEvaluator(const EvaluatorPtr ev)
{
    _evaluators.clear();
    addEvaluator(ev);
}


int NOMAD::MainStep::getRunFlag() const
{

    bool hasFeas = NOMAD::CacheBase::getInstance()->hasFeas(NOMAD::EvalType::BB,NOMAD::ComputeType::STANDARD);
    bool hasInfeas = NOMAD::CacheBase::getInstance()->hasInfeas(NOMAD::EvalType::BB,NOMAD::ComputeType::STANDARD); // hasFeas and hasInfeas can be both false. Maybe cache contains no valid BB evaluations -> initialization error

    bool initializationError = AllStopReasons::testIf(BaseStopType::INITIALIZATION_FAILED);
    if (initializationError && !(hasFeas || hasInfeas))
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Failed initialization detected but cache contains a valid evaluation point.");
    }

    if (nullptr == _stopReasons)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Stop reasons is null. The function getRunFlag must be called after MainStep::end.");
    }

    bool stopByUser = (
                     AllStopReasons::testIf(BaseStopType::CTRL_C) ||
                     AllStopReasons::testIf(BaseStopType::USER_GLOBAL_STOP));
    bool nomadError = (
                       AllStopReasons::testIf(BaseStopType::ERROR) ||
                       AllStopReasons::testIf(BaseStopType::UNKNOWN_STOP_REASON));

    bool stopOnEvalIterCrit =  (
                            AllStopReasons::testIf(EvalGlobalStopType::MAX_EVAL_REACHED) ||
                            AllStopReasons::testIf(EvalGlobalStopType::MAX_BB_EVAL_REACHED) ||
                            AllStopReasons::testIf(EvalGlobalStopType::MAX_BLOCK_EVAL_REACHED) ||
                            _stopReasons->testIf(IterStopType::MAX_ITER_REACHED) );

    bool stopOnTimeLimit =  ( AllStopReasons::testIf(NOMAD::BaseStopType::MAX_TIME_REACHED));

    bool stopOnFeasible = _stopReasons->testIf(IterStopType::STOP_ON_FEAS);

    auto madsStopReason = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get (_stopReasons);
    bool madsStopOnMeshCrit = false;
    if (nullptr != madsStopReason)
    {
      madsStopOnMeshCrit = ( madsStopReason->testIf(NOMAD::MadsStopType::MESH_PREC_REACHED) ||
                             madsStopReason->testIf(NOMAD::MadsStopType::MIN_MESH_INDEX_REACHED) ||
                             madsStopReason->testIf(NOMAD::MadsStopType::MIN_MESH_SIZE_REACHED));
    }

    bool targetReachedCrit = false; // Stop when matching target not yet implemented

    int runFlag = -3 ;
    if (stopByUser) {  // CTRL-C or User stopped (callback function)
        runFlag = -5;
    }
    else if (nomadError) { // Error
        runFlag = -3;
    }
    else if (stopOnTimeLimit) {  // Time limit reached (user option)
        runFlag = -4;
    }
    else if (stopOnFeasible) {  // Stop on feasible point (user option)
        runFlag = -6;
    }
    else if (initializationError) { // Initial point failed to evaluate
        runFlag = -3;
    }
    else if (targetReachedCrit || (hasFeas && madsStopOnMeshCrit)) { // Objective target reached OR Mads converged (mesh criterion) to a feasible point. The true problem is considered (outputs from blackbox evaluations, not surrogate).
        runFlag = 1;
    }
    else if (hasFeas && stopOnEvalIterCrit) { // At least one feasible point obtained and evaluation budget (single bb or block of bb) spent or max iteration (user option) reached.
        runFlag = 0;
    }
    else if (hasInfeas && madsStopOnMeshCrit) { // Mads mesh converged but no feasible point obtained (only infeasible). The true problem is considered (outputs from blackbox evaluations, not surrogate).
        runFlag = -1;
    }
    else if (hasInfeas && stopOnEvalIterCrit) { // No feasible point obtained (only infeasible) and evaluation budget (single bb or block of bb) spent or max iteration (user option) reached
        runFlag = -2;
    }
    else
    {
      // Something else must have happened.
        runFlag = -3;
    }

    return runFlag;
}
