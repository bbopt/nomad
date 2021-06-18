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

#include <signal.h>

#include "../Algos/Algorithm.hpp"
#include "../Algos/EvcInterface.hpp"
#include "../Algos/Mads/SearchMethodBase.hpp"
#include "../Algos/SubproblemManager.hpp"
#include "../Cache/CacheBase.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Math/RNG.hpp"
#include "../Util/fileutils.hpp"

#ifdef TIME_STATS
#include "../Util/Clock.hpp"
#endif // TIME_STATS

void NOMAD::Algorithm::init()
{
    //_name = "AGenericAlgorithmHasNoName";

    // Verifications that throw Exceptions to the Constructor if not validated.
    verifyParentNotNull();

    if (nullptr == _runParams)
    {
        throw NOMAD::StepException(__FILE__, __LINE__,
                               "A valid RunParameters must be provided to an Algorithm constructor.", this);
    }

    if (nullptr == _pbParams)
    {
        throw NOMAD::StepException(__FILE__, __LINE__,
                               "A valid PbParameters must be provided to the Algorithm constructor.", this);
    }

    if ( nullptr == _stopReasons )
        throw NOMAD::StepException(__FILE__, __LINE__,
                               "Valid stop reasons must be provided to the Algorithm constructor.", this);

    // Check pbParams if needed, ex. if a copy of PbParameters was given to the Algorithm constructor.
    _pbParams->checkAndComply();

    // Instantiate generic algorithm termination
    _termination    = std::make_unique<NOMAD::Termination>( this );

    // Update SubproblemManager
    NOMAD::Point fullFixedVariable = isRootAlgo() ? _pbParams->getAttributeValue<NOMAD::Point>("FIXED_VARIABLE")
                                   : NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(_parentStep);

    NOMAD::Subproblem subproblem(_pbParams, fullFixedVariable);
    NOMAD::SubproblemManager::getInstance()->addSubproblem(this, subproblem);
    _pbParams = subproblem.getPbParams();
    _pbParams->checkAndComply();

    /** Step::userInterrupt() will be called if CTRL-C is pressed.
     * Currently, the main thread will wait for all evaluations to be complete.
     * \todo Propage interruption to all threads, for all parallel evaluations of blackbox.
     */
    signal(SIGINT, userInterrupt);
    signal(SIGSEGV, debugSegFault);

}


NOMAD::Algorithm::~Algorithm()
{
    NOMAD::SubproblemManager::getInstance()->removeSubproblem(this);
}


void NOMAD::Algorithm::startImp()
{
#ifdef TIME_STATS
    if (isRootAlgo())
    {
        _startTime = NOMAD::Clock::getCPUTime();
    }
#endif // TIME_STATS

    // All stop reasons are reset.
    _stopReasons->setStarted();

    // SuccessType is reset
    _algoSuccessful = false;
    _algoBestSuccess = NOMAD::SuccessType::NOT_EVALUATED;

    if (isRootAlgo())
    {
        // Update hot restart info
        readInformationForHotRestart();
        NOMAD::CacheBase::getInstance()->setStopWaiting(false);
    }

    // By default reset the lap counter for BbEval and set the lap maxBbEval to INF
    NOMAD::EvcInterface::getEvaluatorControl()->resetLapBbEval();
    NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( NOMAD::INF_SIZE_T );
    NOMAD::EvcInterface::getEvaluatorControl()->resetModelEval();

    if (nullptr == _megaIteration)
    {
        // Default behavior - not hot restart.
        // Clear cache hits.
        // Initialization.
        // Eval X0s.

        if (isRootAlgo())
        {
            // Ensure we do not count cache hits which may have been read in the cache.
            NOMAD::CacheBase::resetNbCacheHits();
        }

        // Perform algo initialization only when available.
        if (nullptr != _initialization)
        {
            _initialization->start();
            _initialization->run();
            _initialization->end();
        }

    }
    else
    {
        // Hot restart situation.
        // We will not need Initialization.
        auto barrier = _megaIteration->getBarrier();

        // Update X0s
        // Use best points.
        auto bestPoints = barrier->getAllPoints();

        NOMAD::ArrayOfPoint x0s;
        if (!bestPoints.empty())
        {
            std::transform(bestPoints.begin(), bestPoints.end(), std::back_inserter(x0s),
                       [](NOMAD::EvalPoint evalPoint) -> NOMAD::EvalPoint { return evalPoint; });

        }
        _pbParams->setAttributeValue<NOMAD::ArrayOfPoint>("X0", x0s);
        _pbParams->checkAndComply();
    }
}


void NOMAD::Algorithm::endImp()
{
    if ( _endDisplay )
    {
        displayBestSolutions();
#ifdef TIME_STATS
        if (isRootAlgo())
        {
            _totalRealAlgoTime = NOMAD::Clock::getTimeSinceStart();
            _totalCPUAlgoTime += NOMAD::Clock::getCPUTime() - _startTime;
        }
#endif // TIME_STATS

        displayEvalCounts();
    }

    // Update the SearchMethod success type with best success found.
    if ( _algoSuccessful )
    {
        // The parent can be a SearchMethod (NM-Mads Search) or not (that is NM is a standalone optimization)
        auto searchMethodConst = dynamic_cast<const NOMAD::SearchMethodBase*>(_parentStep);

        if (searchMethodConst != nullptr)
        {
            auto searchMethod = const_cast<NOMAD::SearchMethodBase*>(searchMethodConst);
            searchMethod->setSuccessType(_algoBestSuccess);
        }

    }

    // By default reset the lap counter for BbEval and set the lap maxBbEval to INF
    NOMAD::EvcInterface::getEvaluatorControl()->resetLapBbEval();
    NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( NOMAD::INF_SIZE_T );

    if (isRootAlgo())
    {
        saveInformationForHotRestart();
        NOMAD::CacheBase::getInstance()->setStopWaiting(true);
    }
}


void NOMAD::Algorithm::hotRestartOnUserInterrupt()
{
#ifdef TIME_STATS
    if (isRootAlgo())
    {
        _totalCPUAlgoTime += NOMAD::Clock::getCPUTime() - _startTime;
    }
#endif // TIME_STATS
    hotRestartBeginHelper();

    hotRestartEndHelper();
#ifdef TIME_STATS
    if (isRootAlgo())
    {
        _startTime = NOMAD::Clock::getCPUTime();
    }
#endif // TIME_STATS
}


void NOMAD::Algorithm::saveInformationForHotRestart() const
{
    // If we want to stop completely and then be able to restart
    // from where we were, we need to save some information on file.
    //
    // Issue 372: Maybe we need to write current parameters. If we write them,
    // ignore initial values, only take latest values down the Parameter tree.
    // For now, using initial parameters.

    // Cache file is treated independently from hot restart file.
    // As long as the cache file name is set, it is written.
    // This is the behavior of NOMAD 3.
    std::string cacheFile = NOMAD::CacheBase::getInstance()->getFileName();
    if (!cacheFile.empty())
    {
        NOMAD::CacheBase::getInstance()->write();
    }
    if ( _runParams->getAttributeValue<bool>("HOT_RESTART_WRITE_FILES"))
    {
        std::cout << "Save information for hot restart." << std::endl;
        std::cout << "Write hot restart file." << std::endl;
        NOMAD::write(*this, _runParams->getAttributeValue<std::string>("HOT_RESTART_FILE"));
    }
}


void NOMAD::Algorithm::displayBestSolutions() const
{
    std::vector<NOMAD::EvalPoint> evalPointList;
    // Display best feasible solutions.
    std::string sFeas;
    // Output level is very high if there are no parent algorithm
    // Output level is info if this algorithm is a sub part of another algorithm.
    NOMAD::OutputLevel outputLevel = isSubAlgo() ? NOMAD::OutputLevel::LEVEL_INFO
                                                 : NOMAD::OutputLevel::LEVEL_VERY_HIGH;
    auto solFormat = NOMAD::OutputQueue::getInstance()->getSolFormat();
    auto computeType = NOMAD::EvcInterface::getEvaluatorControl()->getComputeType();
    auto surrogateAsBB = NOMAD::EvcInterface::getEvaluatorControl()->getSurrogateOptimization();
    if (isRootAlgo())
    {
        solFormat.set(-1);
    }
    NOMAD::OutputInfo displaySolFeas(getName(), sFeas, outputLevel);
    auto fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this);

    sFeas = "Best feasible solution";
    auto barrier = getMegaIterationBarrier();
    if (nullptr != barrier)
    {
        evalPointList = barrier->getAllXFeas();
        NOMAD::convertPointListToFull(evalPointList, fixedVariable);
    }
    size_t nbBestFeas = evalPointList.size();

    if (0 == nbBestFeas)
    {
        sFeas += ":     Undefined.";
        displaySolFeas.addMsg(sFeas);
    }
    else if (1 == nbBestFeas)
    {
        sFeas += ":     ";
        displaySolFeas.addMsg(sFeas + evalPointList[0].display(computeType,
                                                        solFormat,
                                                        NOMAD::DISPLAY_PRECISION_FULL,
                                                        surrogateAsBB));
    }
    else
    {
        sFeas += "s:    ";
        displaySolFeas.addMsg(sFeas + evalPointList[0].display(computeType,
                                                        solFormat,
                                                        NOMAD::DISPLAY_PRECISION_FULL,
                                                        surrogateAsBB));
    }


    const size_t maxSolCount = 8;
    if (nbBestFeas > 1)
    {
        std::vector<NOMAD::EvalPoint>::const_iterator it;
        size_t solCount = 0;
        for (it = evalPointList.begin(); it != evalPointList.end(); ++it)
        {
            solCount++;
            if (evalPointList.begin() == it)
            {
                continue;   // First element already added
            }
            sFeas = "                            ";
            displaySolFeas.addMsg(sFeas + it->display(computeType, solFormat,
                                                      NOMAD::DISPLAY_PRECISION_FULL,
                                                      surrogateAsBB));
            if (solCount >= maxSolCount)
            {
                // We printed enough solutions already.
                displaySolFeas.addMsg("... A total of " + std::to_string(evalPointList.size()) + " feasible solutions were found.");
                break;
            }
        }
    }

    NOMAD::OutputQueue::Add(std::move(displaySolFeas));

    evalPointList.clear();


    // Display best infeasible solutions.
    std::string sInf;
    NOMAD::OutputInfo displaySolInf(getName(), sInf, outputLevel);
    sInf = "Best infeasible solution";
    if (nullptr != barrier)
    {
        evalPointList = barrier->getAllXInf();
        NOMAD::convertPointListToFull(evalPointList, fixedVariable);
    }
    size_t nbBestInf = evalPointList.size();

    if (0 == nbBestInf)
    {
        sInf += ":   Undefined.";
        displaySolInf.addMsg(sInf);
    }
    else if (1 == nbBestInf)
    {
        sInf += ":   ";
        displaySolInf.addMsg(sInf + evalPointList[0].display(computeType,
                                                        solFormat,
                                                        NOMAD::DISPLAY_PRECISION_FULL,
                                                        surrogateAsBB));
    }
    else
    {
        sInf += "s:  ";
        displaySolInf.addMsg(sInf + evalPointList[0].display(computeType,
                                                        solFormat,
                                                        NOMAD::DISPLAY_PRECISION_FULL,
                                                        surrogateAsBB));
    }

    if (nbBestInf > 1)
    {
        size_t solCount = 0;
        std::vector<NOMAD::EvalPoint>::const_iterator it;
        for (it = evalPointList.begin(); it != evalPointList.end(); ++it)
        {
            solCount++;
            if (evalPointList.begin() == it)
            {
                continue;   // First element already added
            }
            displaySolInf.addMsg("                            " + it->display(computeType,
                                                                        solFormat,
                                                                        NOMAD::DISPLAY_PRECISION_FULL,
                                                                        surrogateAsBB));
            if (solCount >= maxSolCount)
            {
                // We printed enough solutions already.
                displaySolInf.addMsg("... A total of " + std::to_string(evalPointList.size()) + " infeasible solutions were found.");
                break;
            }
        }
    }

    NOMAD::OutputQueue::Add(std::move(displaySolInf));
}


void NOMAD::Algorithm::displayEvalCounts() const
{
    // Display evaluation information

    // Used to display or not certain values
    bool isSub = isSubAlgo();

    // Actual numbers
    size_t bbEval       = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
    size_t lapBbEval    = NOMAD::EvcInterface::getEvaluatorControl()->getLapBbEval();
    size_t nbEval       = NOMAD::EvcInterface::getEvaluatorControl()->getNbEval();
    size_t modelEval    = NOMAD::EvcInterface::getEvaluatorControl()->getModelEval();
    size_t totalModelEval = NOMAD::EvcInterface::getEvaluatorControl()->getTotalModelEval();
    size_t nbCacheHits  = NOMAD::CacheBase::getNbCacheHits();
    int nbEvalNoCount   = static_cast<int>(nbEval - bbEval - nbCacheHits);

    // What needs to be shown, according to the counts and to the value of isSub
    bool showNbEvalNoCount  = (nbEvalNoCount > 0);
    bool showModelEval      = isSub && (modelEval > 0);
    bool showTotalModelEval = (totalModelEval > 0);
    bool showNbCacheHits    = (nbCacheHits > 0);
    bool showNbEval         = (nbEval > bbEval);
    bool showLapBbEval      = isSub && (bbEval > lapBbEval && lapBbEval > 0);

    // Output levels will be modulated depending on the counts and on the Algorithm level.
    NOMAD::OutputLevel outputLevelHigh = isSub ? NOMAD::OutputLevel::LEVEL_INFO
                                               : NOMAD::OutputLevel::LEVEL_HIGH;
    NOMAD::OutputLevel outputLevelNormal = isSub ? NOMAD::OutputLevel::LEVEL_INFO
                                                 : NOMAD::OutputLevel::LEVEL_NORMAL;

    // Padding for nice presentation
    std::string sFeedBbEval, sFeedLapBbEval, sFeedNbEvalNoCount, sFeedModelEval,
                sFeedTotalModelEval, sFeedCacheHits, sFeedNbEval;

    // Conditional values: showNbEval, showNbEvalNoCount, showLapBbEval
    if (showLapBbEval)  // Longest title
    {
        sFeedBbEval += "                 ";
        //sFeedLapBbEval += "";
        sFeedNbEvalNoCount += "   ";
        sFeedModelEval += "                    ";
        sFeedTotalModelEval += "              ";
        sFeedCacheHits += "                           ";
        sFeedNbEval += "          ";
    }
    else if (showNbEvalNoCount) // Second longest
    {
        sFeedBbEval += "              ";
        //sFeedLapBbEval += "";
        //sFeedNbEvalNoCount += "";
        sFeedModelEval += "                 ";
        sFeedTotalModelEval += "           ";
        sFeedCacheHits += "                        ";
        sFeedNbEval += "       ";
    }
    else if (showNbEval)    // 3rd longest title
    {
        sFeedBbEval += "       ";
        //sFeedLapBbEval += "";
        //sFeedNbEvalNoCount += "";
        sFeedModelEval += "          ";
        sFeedTotalModelEval += "    ";
        sFeedCacheHits += "                 ";
        //sFeedNbEval += "";
    }
    else if (showTotalModelEval)
    {
        sFeedBbEval += "   ";
        //sFeedLapBbEval += "";
        //sFeedNbEvalNoCount += "";
        //sFeedModelEval += "          ";
        //sFeedTotalModelEval += "    ";
        //sFeedCacheHits += "                 ";
        //sFeedNbEval += "";
    }

    std::string sBbEval         = "Blackbox evaluations: " + sFeedBbEval + NOMAD::itos(bbEval);
    std::string sLapBbEval      = "Sub-optimization blackbox evaluations: " + sFeedLapBbEval + NOMAD::itos(lapBbEval);
    std::string sNbEvalNoCount  = "Blackbox evaluation (not counting): " + sFeedNbEvalNoCount + NOMAD::itos(nbEvalNoCount);
    std::string sModelEval      = "Model evaluations: " + sFeedModelEval + NOMAD::itos(modelEval);
    std::string sTotalModelEval = "Total model evaluations: " + sFeedTotalModelEval + NOMAD::itos(totalModelEval);
    std::string sCacheHits      = "Cache hits: " + sFeedCacheHits + NOMAD::itos(nbCacheHits);
    std::string sNbEval         = "Total number of evaluations: " + sFeedNbEval + NOMAD::itos(nbEval);


#ifdef TIME_STATS
    std::string sTotalRealTime  = "Total real time (round s):    " + std::to_string(_totalRealAlgoTime);
    std::string sTotalCPUTime   = "Total CPU time (s):           " + std::to_string(_totalCPUAlgoTime);
 #endif // TIME_STATS

    AddOutputInfo("", outputLevelHigh); // skip line
    // Always show number of blackbox evaluations
    AddOutputInfo(sBbEval, outputLevelHigh);
    // The other values are conditional to the show* booleans
    if (showLapBbEval)
    {
        AddOutputInfo(sLapBbEval, outputLevelNormal);
    }
    if (showNbEvalNoCount)
    {
        AddOutputInfo(sNbEvalNoCount, outputLevelNormal);
    }
    if (showModelEval)
    {
        AddOutputInfo(sModelEval, outputLevelNormal);
    }
    if (showTotalModelEval)
    {
        AddOutputInfo(sTotalModelEval, outputLevelNormal);
    }
    if (showNbCacheHits)
    {
        AddOutputInfo(sCacheHits, outputLevelNormal);
    }
    if (showNbEval)
    {
        AddOutputInfo(sNbEval, outputLevelNormal);
    }

#ifdef TIME_STATS
    {
        AddOutputInfo(sTotalRealTime, outputLevelNormal);
        AddOutputInfo(sTotalCPUTime, outputLevelNormal);
    }
#endif // TIME_STATS
}


bool NOMAD::Algorithm::isSubAlgo() const
{
    bool isSub = false;

    auto parentAlgo = getParentOfType<NOMAD::Algorithm*>();
    if (nullptr != parentAlgo)
    {
        isSub = true;
    }

    return isSub;
}


bool NOMAD::Algorithm::terminate(const size_t iteration)
{
    return _termination->terminate(iteration);
}


void NOMAD::Algorithm::display ( std::ostream& os ) const
{

    os << "MEGA_ITERATION " << std::endl;
    os << *_megaIteration << std::endl;
    os << "NB_EVAL " << NOMAD::EvcInterface::getEvaluatorControl()->getNbEval() << std::endl;
    os << "NB_BB_EVAL " << NOMAD::EvcInterface::getEvaluatorControl()->getBbEval() << std::endl;
    uint32_t x, y, z;
    NOMAD::RNG::getPrivateSeed(x, y, z);
    os << "RNG " << x << " " << y << " " << z << std::endl;

}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::Algorithm & mads)
{
    mads.display(os);
    return os;
}


void NOMAD::Algorithm::read(std::istream& is)
{
    // Read line by line
    std::string name;

    int nbEval = 0, nbBbEval = 0;
    uint32_t x, y, z;

    while (is >> name && is.good() && !is.eof())
    {
        if ("MEGA_ITERATION" == name)
        {
            is >> *_megaIteration;
        }
        else if ("NB_EVAL" == name)
        {
            is >> nbEval;
        }
        else if ("NB_BB_EVAL" == name)
        {
            is >> nbBbEval;
        }
        else if ("RNG" == name)
        {
            is >> x >> y >> z;
            NOMAD::RNG::setPrivateSeed(x, y, z);
        }
        else
        {
            // Put back name to istream. Maybe there is a simpler way.
            for (unsigned i = 0; i < name.size(); i++)
            {
                is.unget();
            }
            break;
        }
    }

    NOMAD::EvcInterface::getEvaluatorControl()->setBbEval(nbBbEval);
    NOMAD::EvcInterface::getEvaluatorControl()->setNbEval(nbEval);

}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::Algorithm& algo)
{
    algo.read(is);
    return is;

}
