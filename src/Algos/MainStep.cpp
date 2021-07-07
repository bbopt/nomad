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
#include "../Math/LHS.hpp"
#include "../Math/RNG.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Output/OutputDirectToFile.hpp"
#include "../Util/Clock.hpp"
#include "../Util/fileutils.hpp"

// Specific algos
#include "../Algos/LatinHypercubeSampling/LH.hpp"
#include "../Algos/Mads/Mads.hpp"
#include "../Algos/Mads/MadsIteration.hpp"
#include "../Algos/Mads/Search.hpp"
#include "../Algos/Mads/VNSSearchMethod.hpp"
#include "../Algos/NelderMead/NM.hpp"
#include "../Algos/PhaseOne/PhaseOne.hpp"
#ifdef _OPENMP
#include "../Algos/PSDMads/PSDMads.hpp"
#endif
#ifdef USE_SGTELIB
#include "../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../Algos/SgtelibModel/SgtelibModel.hpp"
#endif
#include "../Algos/SSDMads/SSDMads.hpp"


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
    createCache();

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
        updateX0sFromCache();

        // Need to start evaluator control.
        // Even though there are no BB evaluations done by suggest(),
        // the EvaluatorControl still holds information,
        // and is especially useful for Model searches.
        auto evaluatorControl = NOMAD::EvcInterface::getEvaluatorControl();
        if (nullptr == evaluatorControl)
        {
            if (nullptr == _evaluator)
            {
                // Batch mode. Create Evaluator on the go.
                _evaluator = std::shared_ptr<NOMAD::Evaluator>(
                                   new NOMAD::Evaluator(_allParams->getEvalParams(),
                                                        NOMAD::EvalType::BB,
                                                        NOMAD::EvalXDefined::USE_BB_EVAL));
            }

            std::unique_ptr<NOMAD::EvaluatorControlParameters> evalContParams(new NOMAD::EvaluatorControlParameters(*_allParams->getEvaluatorControlParams()));
            evalContParams->checkAndComply();
            evaluatorControl = std::make_shared<NOMAD::EvaluatorControl>(_evaluator,
                                                     _allParams->getEvaluatorControlGlobalParams(),
                                                     std::move(evalContParams));
            NOMAD::EvcInterface::setEvaluatorControl(std::move(evaluatorControl));
        }

        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
        auto mads = std::make_shared<NOMAD::Mads>(this,
                                              stopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams());

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

    createCache();

    if (evalPointList.size() > 0)
    {
        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
        auto mads = std::make_shared<NOMAD::Mads>(this,
                                              stopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams());
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
        evalPoint.setBBO(fxs[i].display(), bbOutputType);
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

    createCache();
    updateX0sFromCache();

    
    // Setup EvaluatorControl
    setNumThreads();

    if (nullptr == _evaluator)
    {
        // Batch mode. Create Evaluator on the go.
        bool surrogateAsBB = _allParams->getAttributeValue<bool>("EVAL_SURROGATE_OPTIMIZATION");
        auto evalType = (surrogateAsBB) ? NOMAD::EvalType::SURROGATE : NOMAD::EvalType::BB;
        _evaluator = std::shared_ptr<NOMAD::Evaluator>(
                        new NOMAD::Evaluator(_allParams->getEvalParams(),
                                             evalType,
                                             NOMAD::EvalXDefined::USE_BB_EVAL));
    }

    auto evaluatorControl = NOMAD::EvcInterface::getEvaluatorControl();
    if (nullptr == evaluatorControl)
    {
        std::unique_ptr<NOMAD::EvaluatorControlParameters> evalContParams(new NOMAD::EvaluatorControlParameters(*_allParams->getEvaluatorControlParams()));
        evalContParams->checkAndComply();
        evaluatorControl = std::make_shared<NOMAD::EvaluatorControl>(_evaluator,
                                                                      _allParams->getEvaluatorControlGlobalParams(),
                                                                      std::move(evalContParams));
        NOMAD::EvcInterface::setEvaluatorControl(std::move(evaluatorControl));
    }

    // Currently this does nothing.
    NOMAD::EvcInterface::getEvaluatorControl()->start();

    // Create the Algorithms that we want to solve.
    // Currently available: PhaseOne, Mads, LH, NM, SgtelibModelEval,
    // QuadModelOptimization
    // Note: These algorithm parameters are mutually exclusive:
    // LH_EVAL NM_OPTIMIZATION PSD_MADS_OPTIMIZATION QUAD_MODEL_OPTIMIZATION
    // SGTELIB_MODEL_EVAL SSD_MADS_OPTIMIZATION.
    // This is caught by checkAndComply().
    _algos.clear();

    auto nbLHEval = _allParams->getRunParams()->getAttributeValue<size_t>("LH_EVAL");
    auto doNMOptimization = _allParams->getRunParams()->getAttributeValue<bool>("NM_OPTIMIZATION");
#ifdef _OPENMP
    bool doPSDMads = _allParams->getRunParams()->getAttributeValue<bool>("PSD_MADS_OPTIMIZATION");
#endif
#ifdef USE_SGTELIB
    bool doQuadModelOpt = _allParams->getRunParams()->getAttributeValue<bool>("QUAD_MODEL_OPTIMIZATION");
    bool doSgtelibModelEval = _allParams->getRunParams()->getAttributeValue<bool>("SGTELIB_MODEL_EVAL");
#endif
    bool doSSDMads = _allParams->getRunParams()->getAttributeValue<bool>("SSD_MADS_OPTIMIZATION");

    if ( nbLHEval > 0 )
    {
        auto lhStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::LHStopType>>();

        // All the LH sample points must be evaluated. No opportunism.
        if ( _allParams->getAttributeValue<bool>("EVAL_OPPORTUNISTIC") )
            AddOutputInfo("Opportunistic evaluation is disabled for LH when ran as a single algorithm.");

        NOMAD::EvcInterface::getEvaluatorControl()->setOpportunisticEval(false);

        auto lh = std::make_shared<NOMAD::LH>(this,
                                              lhStopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams());
        _algos.push_back(lh);
    }
    else if ( doNMOptimization )
    {
        auto nmStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::NMStopType>>();

        // All the NM points must be evaluated. No opportunism.
        if ( _allParams->getAttributeValue<bool>("EVAL_OPPORTUNISTIC") )
        {
            AddOutputInfo("Opportunistic evaluation is disabled for NM when ran as a single algorithm.");
        }

        NOMAD::EvcInterface::getEvaluatorControl()->setOpportunisticEval(false);

        auto nm = std::make_shared<NOMAD::NM>(this,
                                              nmStopReasons ,
                                              _allParams->getRunParams(),
                                              _allParams->getPbParams());
        _algos.push_back(nm);
    }
#ifdef _OPENMP
    else if (doPSDMads)
    {
        auto PSDMadsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>(); // A PSD-MADS is a MADS with respect to stop reasons.
        auto psd = std::make_shared<NOMAD::PSDMads>(this,
                                                    _evaluator,
                                                    _allParams->getEvaluatorControlParams(),
                                                    PSDMadsStopReasons,
                                                    _allParams->getRunParams(),
                                                    _allParams->getPbParams());
        _algos.push_back(psd);
    }
#endif
#ifdef USE_SGTELIB
    else if (doQuadModelOpt)
    {
        auto quadModelStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::ModelStopType>>();

        // All the Sgtelib Model sample points are evaluated sequentially. No opportunism.
        NOMAD::EvcInterface::getEvaluatorControl()->setOpportunisticEval(false);
        _allParams->setAttributeValue("MEGA_SEARCH_POLL", false);
        _allParams->checkAndComply();

        auto quadModelAlgo = std::make_shared<NOMAD::QuadModelAlgo>(this,
                                                        quadModelStopReasons,
                                                        _allParams->getRunParams(),
                                                        _allParams->getPbParams());
        _algos.push_back(quadModelAlgo);
    }
    else if (doSgtelibModelEval)
    {
        auto sgtelibModelStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::ModelStopType>>();

        // All the Sgtelib Model sample points are evaluated. No opportunism.
        NOMAD::EvcInterface::getEvaluatorControl()->setOpportunisticEval(false);


        std::shared_ptr<NOMAD::Barrier> barrier = nullptr;
        if (NOMAD::CacheBase::getInstance()->size() > 0)
        {
            barrier = std::make_shared<NOMAD::Barrier>();
        }
        auto sgtelibModel = std::make_shared<NOMAD::SgtelibModel>(this,
                                                        sgtelibModelStopReasons,
                                                        barrier,
                                                        _allParams->getRunParams(),
                                                        _allParams->getPbParams(),
                                                        nullptr);   // no mesh
        _algos.push_back(sgtelibModel);
    }
#endif // USE_SGTELIB
    else if (doSSDMads)
    {
        auto SSDMadsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>(); // A SSD-MADS is a MADS with respect to stop reasons.
        auto ssd = std::make_shared<NOMAD::SSDMads>(this,
                                                    SSDMadsStopReasons,
                                                    _allParams->getRunParams(),
                                                    _allParams->getPbParams());
        _algos.push_back(ssd);
    }
    else
    {
        // The stop reasons for mads
        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

        if (detectPhaseOne())
        {
            // First, run PhaseOne, which has its own Mads.
            // Then, run regular Mads.

            auto phaseOneStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::PhaseOneStopType>>();
            auto phaseOne = std::make_shared<NOMAD::PhaseOne>(this,
                                                              phaseOneStopReasons,
                                                              _allParams->getRunParams(),
                                                              _allParams->getPbParams());
            // Ensure PhaseOne does not show found solutions
            phaseOne->setEndDisplay(false);

            _algos.push_back(phaseOne);

            auto mads = std::make_shared<NOMAD::Mads>(this,
                                                      stopReasons ,
                                                      _allParams->getRunParams(),
                                                      _allParams->getPbParams());
            // Get MegaIteration from PhaseOne to continue optimization
            mads->setMegaIteration(phaseOne->getMegaIteration());
            _algos.push_back(mads);
        }
        else
        {
            // Default behavior: create Mads and add it to _algos.
            auto mads = std::make_shared<NOMAD::Mads>(this,
                                                      stopReasons ,
                                                      _allParams->getRunParams(),
                                                      _allParams->getPbParams());
            _algos.push_back(mads);
        }
    }
    


}


bool NOMAD::MainStep::runImp()
{
    bool ret = false;
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    evc->restart();

    for (auto algo : _algos)
    {
        // Note: Algo start could be moved from outside to inside the parallel region
        // so that evaluation of multiple X0s is done in parallel.
        // Currently this may cause an issue with evc parameters (may be get before
        // they are checked).
        algo->start();

        // Begin parallel region.
        // Note: always use default(none) and never default(shared), which is much
        // too risky.
#ifdef _OPENMP
#pragma omp parallel default(none) shared(algo,evc,ret)
#endif // _OPENMP
        {
            printNumThreads();

            // Start evaluatorControl on all threads
            evc->run();

            if (evc->isMainThread(NOMAD::getThreadNum()))
            {
                // Algo run is done in main thread(s) only.
                ret = algo->run();

                // When algo is done, evaluatorControl will wait for all main threads to be done.
                evc->stop();
            }
        }   // End of parallel region.
        algo->end();

        if (algo->getAllStopReasons()->checkTerminate())
        {
            break;
        }
     }

    return ret;
}


void NOMAD::MainStep::endImp()
{
    _algos.clear();

#ifdef TIME_STATS
   _totalCPUTime += NOMAD::Clock::getCPUTime() - _startTime;
#endif // TIME_STATS

    displayDetailedStats();
}


int NOMAD::MainStep::getNumThreads() const
{
    int nbThreadsParam = 1;
#ifdef _OPENMP
    // Number of threads on which to do the evaluation.
    nbThreadsParam = _allParams->getAttributeValue<int>("NB_THREADS_OPENMP");
    if (nbThreadsParam < 1)
    {
        nbThreadsParam = omp_get_max_threads();
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
#endif // _OPENMP
}


void NOMAD::MainStep::printNumThreads() const
{
#ifdef _OPENMP
    // Once we are in the parallel region, print the actual
    // number of threads used.
#pragma omp single nowait
    {
        int nbThreads = omp_get_num_threads();
        std::string s = "Using " + NOMAD::itos(nbThreads) + " thread";
        if (nbThreads > 1) { s += "s"; }
        s += ".";
        NOMAD::OutputQueue::Add(s);
    }
#endif // _OPENMP
}


// Detect that we must do a Phase One.
bool NOMAD::MainStep::detectPhaseOne()
{
    bool hasEBConstraints = false;
    bool hasNoFeas = !NOMAD::CacheBase::getInstance()->hasFeas();

    auto bbOutputTypeList = _allParams->getEvalParams()->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    if (std::find(bbOutputTypeList.begin(), bbOutputTypeList.end(), NOMAD::BBOutputType::EB)
        != bbOutputTypeList.end())
    {
        hasEBConstraints = true;
    }

    return hasEBConstraints && hasNoFeas;
}


// Helper for start
void NOMAD::MainStep::createCache() const
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
    }
}


void NOMAD::MainStep::updateX0sFromCache() const
{
    // Update X0s, if needed.
    auto x0s = _allParams->getPbParams()->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
    if (x0s.empty() || x0s[0].toBeDefined())
    {
        // X0 not provided directly by user. Find them in cache, or generate them using LHS.
        x0s.clear();

        bool canUseCache = (NOMAD::CacheBase::getInstance()->size() > 0);
        auto lhSearchType = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
        bool canUseLH = (lhSearchType.isEnabled() && lhSearchType.getNbInitial() > 0);

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
                                                computeType, nullptr);
            if (0 == evalPointList.size())
            {
                auto hMax = _allParams->getRunParams()->getAttributeValue<NOMAD::Double>("H_MAX_0");
                NOMAD::CacheBase::getInstance()->findBestInf(evalPointList,
                                                    hMax, fixedVariable,
                                                    evalType, computeType, nullptr);
            }
            if (0 == evalPointList.size())
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "Cache did not find a best point to initialize X0");
            }
            for (size_t i = 0; i < evalPointList.size(); i++)
            {
                x0s.push_back(evalPointList[i]);
            }
        }
        else if (canUseLH)
        {
            x0s = suggestFromLH(lhSearchType.getNbInitial());
        }

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
    // Note: The "Beta" information is not part of the NOMAD_VERSION_NUMBER.
    //version += " Beta 2";
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
    /*
    std::string info;
    // This file has been removed. See Issue #618.
    std::string filename = "Util/Copyright.hpp";
    if (NOMAD::readAllFile(info, filename))
    {
        NOMAD::OutputQueue::Add(info, NOMAD::OutputLevel::LEVEL_ERROR);
    }
    */
}


/*------------------------------------------------------*/
/*             display NOMAD help on a subject          */
/*------------------------------------------------------*/
void NOMAD::MainStep::displayHelp( const std::string & helpSubject , bool devHelp )
{
    _allParams->displayHelp( helpSubject, devHelp, std::cout );
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

    if (!_userTerminate)
    {
        std::cout << "Hot restart" ;

        // Do not use the shared_ptr _evaluator because it is NULL in this function
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
                while(!_userTerminate && std::getline(std::cin, line))
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

    // Reset static tag counter
    NOMAD::EvalPoint::resetCurrentTag();
    // Reset SubproblemManager map
    NOMAD::SubproblemManager::getInstance()->reset();
    // Reset evaluator control
    NOMAD::EvcInterface::resetEvaluatorControl();
    // Reset seed
    NOMAD::RNG::resetPrivateSeedToDefault();
    // Reset parameter entries
    NOMAD::Parameters::eraseAllEntries();
    // Reset VNS Mads search
    NOMAD::VNSSearchMethod::reset();
}

void NOMAD::MainStep::resetCache()
{
    // Get a new cache
    NOMAD::CacheBase::resetInstance(); // Need to reset the singleton. When calling createCache there is no instance and we are sure to call NOMAD::CacheSet::setInstance from scratch. The cache file is read and the cache is set with its content.
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
        std::cerr << "Warning: could not open evaluation stats file " << evalStatsFile << std::endl;
    }
    evalStatsStream << paddedStats << std::endl;
    evalStatsStream.close();
}

