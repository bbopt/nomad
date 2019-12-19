/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/**
 \file   MainStep.cpp
 \brief  Main Step to hold MADS
 \author Viviane Rochon Montplaisir
 \date   June 2018
 */
#include "../Cache/CacheSet.hpp"

// Generic
#include "../Algos/CacheInterface.hpp"
#include "../Algos/EvcInterface.hpp"
#include "../Algos/MainStep.hpp"

// Specific algos
#include "../Algos/LatinHypercubeSampling/LH.hpp"
#include "../Algos/Mads/Mads.hpp"
#ifdef USE_SGTELIB
#include "../Algos/SgtelibModel/SgtelibModel.hpp"
#endif
#include "../Algos/PhaseOne/PhaseOne.hpp"
#include "../Algos/NelderMead/NM.hpp"

#include "../Util/Clock.hpp"
#include "../Util/fileutils.hpp"


// Initialization of static members
std::string NOMAD::MainStep::_algoComment = "";
std::vector<std::string> NOMAD::MainStep::_prevAlgoComment = std::vector<std::string>();
bool NOMAD::MainStep::_forceAlgoComment = false;


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

    _name = "Main";

    // Start the clock
    NOMAD::Clock::reset();
}


NOMAD::MainStep::~MainStep()
{
    _algos.clear();
    _subproblems.clear();
}


// Set _algoComment, and push the current algo comment to _prevAlgoComment pile.
// Do not push an empty string on an empty pile, it is irrelevant.
// If force = true, do not set a new algo comment until this comment is reset.
void NOMAD::MainStep::setAlgoComment(const std::string& algoComment, const bool force)
{
    if (!_forceAlgoComment)
    {
        // Push algo comment to _prevAlgoComment
        if (!_prevAlgoComment.empty() || !_algoComment.empty())
        {
            _prevAlgoComment.push_back(_algoComment);
        }
        _algoComment = algoComment;
    }
    if (force)
    {
        _forceAlgoComment = true;
    }
}


// Pop the previous algo comment from the _prevAlgoComment pile.
void NOMAD::MainStep::resetPreviousAlgoComment(const bool force)
{   
    if (!_forceAlgoComment || force)
    {
        if (_prevAlgoComment.empty())
        {
            _algoComment = "";
        }
        else
        {
            // Remove last element, simulate a "pop".
            _algoComment = std::move(_prevAlgoComment[_prevAlgoComment.size()-1]);
            _prevAlgoComment.erase(_prevAlgoComment.end()-1);
        }
        if (_forceAlgoComment)
        {
            _forceAlgoComment = false;
        }
    }
}


std::shared_ptr<NOMAD::Subproblem> NOMAD::MainStep::getCurrentSubproblem() const
{
    std::shared_ptr<NOMAD::Subproblem> currentSub(nullptr);

    if (_subproblems.size() >= 1)
    {
        currentSub = std::make_shared<NOMAD::Subproblem>(_subproblems[0]);   // This has to be generalized
    }
    else
    {
        // No subproblem defined.
        // Create a default subproblem
        currentSub = std::make_shared<NOMAD::Subproblem>(_pbParams);
    }

    return currentSub;
}


void NOMAD::MainStep::startImp()
{

    if (nullptr == _allParams)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Using Library mode. Parameters must be set prior to running MainStep step.");
    }

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

    createCache();
    updateX0sFromCache();

    // Setup EvaluatorControl
    setNumThreads();

    if (nullptr == _evaluator)
    {
        // Batch mode. Create Evaluator on the go.
        _evaluator = std::unique_ptr<NOMAD::Evaluator>(
                        new NOMAD::Evaluator(_allParams->getEvalParams(), 
                                             NOMAD::EvalType::BB,
                                             getNumThreads(),
                                             NOMAD::EvalXDefined::USE_BB_EVAL));
    }

    auto evaluatorControl = std::make_shared<NOMAD::EvaluatorControl>(std::move(_evaluator),
                                                                      _allParams->getEvaluatorControlParams() );
    NOMAD::EvcInterface::setEvaluatorControl(evaluatorControl);

    // Currently this does nothing.
    NOMAD::EvcInterface::getEvaluatorControl()->start();

    // TODO Create multiple subproblems / groups of variables.
    // Currently, we manage a single subproblem, based on the parameter
    // FIXED_VARIABLE. This Subproblem is created even when FIXED_VARIABLE
    // is not set.
    // The Subproblem class creates a PbParameters for the subproblem. The
    // algorithms will be restricted to these subdimension parameters.
    _subproblems.clear();
    auto subproblem = NOMAD::Subproblem(_allParams->getPbParams());
    _subproblems.push_back(subproblem);

    // Create the Algorithms that we want to solve.
    // Currently available: PhaseOne, Mads, LH, NM, SgtelibModel
    // Note: These algorithm parameters are mutually exclusive:
    // LH_EVAL NM_OPTIMIZATION SGTELIB_MODEL_EVAL.
    // This is caught by checkAndComply().
    _algos.clear();

    auto nbLHEval = _allParams->getRunParams()->getAttributeValue<size_t>("LH_EVAL");
    auto doNMOptimization = _allParams->getRunParams()->getAttributeValue<bool>("NM_OPTIMIZATION");
#ifdef USE_SGTELIB
    bool doSgtelibModelEval = _allParams->getRunParams()->getAttributeValue<bool>("SGTELIB_MODEL_EVAL");
#endif

    if ( nbLHEval > 0 )
    {
        auto lhStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::LHStopType>>();
        
        // All the LH sample points must be evaluated. No opportunism.
        if ( _allParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL") )
            AddOutputInfo("Opportunistic evaluation is disabled for LH when ran as a single algorithm.");
        
        _allParams->setAttributeValue("OPPORTUNISTIC_EVAL",false);
        _allParams->checkAndComply( );
        
        auto lh = std::make_shared<NOMAD::LH>(this,
                                              lhStopReasons ,
                                              _allParams->getRunParams(),
                                              getCurrentSubproblem()->getPbParams());
        _algos.push_back(lh);
    }
    else if ( doNMOptimization )
    {
        auto nmStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::NMStopType>>();
        
        // All the NM points must be evaluated. No opportunism.
        if ( _allParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL") )
            AddOutputInfo("Opportunistic evaluation is disabled for NM when ran as a single algorithm.");
        
        _allParams->setAttributeValue("OPPORTUNISTIC_EVAL",false);
        _allParams->checkAndComply( );
        
        auto nm = std::make_shared<NOMAD::NM>(this,
                                              nmStopReasons ,
                                              _allParams->getRunParams(),
                                              getCurrentSubproblem()->getPbParams());
        _algos.push_back(nm);
    }
#ifdef USE_SGTELIB
    else if (doSgtelibModelEval)
    {
        auto sgtelibModelStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::SgtelibModelStopType>>();

        // All the Sgtelib Model sample points are evaluated. No opportunism.
        _allParams->setAttributeValue("OPPORTUNISTIC_EVAL", false);
        _allParams->checkAndComply();

        std::shared_ptr<NOMAD::Barrier> barrier = nullptr;
        if (NOMAD::CacheBase::getInstance()->size() > 0)
        {
            barrier = std::make_shared<NOMAD::Barrier>();
        }
        auto sgtelibModel = std::make_shared<NOMAD::SgtelibModel>(this,
                                                        sgtelibModelStopReasons,
                                                        barrier,
                                                        _allParams->getRunParams(),
                                                        getCurrentSubproblem()->getPbParams(),
                                                        nullptr);   // no mesh
        _algos.push_back(sgtelibModel);
    }
#endif // USE_SGTELIB
    else
    {
        // The stop reasons for mads
        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

        if (detectPhaseOne())
        {
            // First, run PhaseOne, which has its own Mads.
            // Then, run regular Mads.

            auto PhaseOneStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::PhaseOneStopType>>();
            auto phaseOne = std::make_shared<NOMAD::PhaseOne>(this,
                                                              PhaseOneStopReasons,
                                                              _allParams->getRunParams(),
                                                              getCurrentSubproblem()->getPbParams());
            NOMAD::PhaseOne::setBBOutputTypes(_allParams->getEvalParams()->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"));
            _algos.push_back(phaseOne);

            // TODO At this point, we may want to use setMegaIteration(), e.g. to ensure the MAX_ITERATION is not over.
            // Based on the code for hot restart:
            /*
             // Create a GMesh and an MegaIteration with default values, to be filled
             // by istream is.
             auto barrier = std::make_shared<NOMAD::Barrier>();
             int k = 0;
             auto mesh = std::shared_ptr<NOMAD::MeshBase>(new NOMAD::GMesh(_allParams->getPbParams()));
             NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;
             auto megaIter = std::shared_ptr<NOMAD::MegaIteration>(new NOMAD::MegaIteration(this, k, barrier, mesh, success));
             setMegaIteration(megaIter);
             */

            auto mads = std::make_shared<NOMAD::Mads>(this,
                                                      stopReasons ,
                                                      _allParams->getRunParams(),
                                                      getCurrentSubproblem()->getPbParams());
            _algos.push_back(mads);
        }
        else
        {
            // Default behavior: create Mads and add it to _algos.
            auto mads = std::make_shared<NOMAD::Mads>(this,
                                                      stopReasons ,
                                                      _allParams->getRunParams(),
                                                      getCurrentSubproblem()->getPbParams());
            _algos.push_back(mads);
        }
    }

}


bool NOMAD::MainStep::runImp()
{
    bool ret = false;

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
#pragma omp parallel default(none) shared(algo,ret)
#endif // _OPENMP
        {
            printNumThreads();

            // Start evaluatorControl on all threads
            NOMAD::EvcInterface::getEvaluatorControl()->restart();
            NOMAD::EvcInterface::getEvaluatorControl()->run();

            // Start algorithms on master thread only.
#ifdef _OPENMP
#pragma omp master
#endif // _OPENMP
            {
                // Algo run is done in master thread only.
                ret = algo->run();
                // When algo is done, evaluatorControl is ready to quit all threads.
                NOMAD::EvcInterface::getEvaluatorControl()->stop();
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
#pragma omp master
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
                                     _allParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"));
    }
}


void NOMAD::MainStep::updateX0sFromCache() const
{
    // Update X0s, if needed.
    auto x0s = _allParams->getPbParams()->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
    if ((x0s.empty() || x0s[0].toBeDefined())
        && NOMAD::CacheBase::getInstance()->size() > 0)
    {
        x0s.clear();
        // Use best points in cache for x0.
        std::vector<NOMAD::EvalPoint> evalPointList;
        // Note: We are working in full dimension here, not in subproblem.
        // For this reason, use cache instance directly, not CacheInterface.
        auto fixedVariable = _allParams->getPbParams()->getAttributeValue<NOMAD::Point>("FIXED_VARIABLE");
        NOMAD::CacheBase::getInstance()->findBestFeas(evalPointList,
                                                fixedVariable, getEvalType());
        if (0 == evalPointList.size())
        {
            auto hMax = _allParams->getRunParams()->getAttributeValue<NOMAD::Double>("H_MAX_0");
            NOMAD::CacheBase::getInstance()->findBestInf(evalPointList,
                                                hMax, fixedVariable, getEvalType());
        }
        if (0 == evalPointList.size())
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Cache did not find a best point to initialize X0");
        }
        for (size_t i = 0; i < evalPointList.size(); i++)
        {
            x0s.push_back(evalPointList[i]);
        }

        _allParams->getPbParams()->setAttributeValue("X0", x0s);
        _allParams->getPbParams()->checkAndComply();
    }

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

    // TODO          + "Developer help : " + strExeName + " -d keyword(s) (or 'all')\n"

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
    version += " Beta 1";
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
    std::string filename = "Util/Copyright.hpp";
    if (NOMAD::readAllFile(info, filename))
    {
        NOMAD::OutputQueue::Add(info, NOMAD::OutputLevel::LEVEL_ERROR);
    }
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

        // Do not use the unique_ptr _evaluator because it is NULL in this function
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

        std::cin.clear();
    }

    hotRestartEndHelper();
}
