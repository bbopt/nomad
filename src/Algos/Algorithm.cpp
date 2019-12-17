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

#include <signal.h>

#include "../Algos/EvcInterface.hpp"
#include "../Algos/Algorithm.hpp"

#include "../Math/RNG.hpp"

#include "../Param/AllParameters.hpp"  // Used for hot restart only.

#include "../Util/fileutils.hpp"


void NOMAD::Algorithm::init()
{

    _name = "AGenericAlgorithmHasNoName";
    
    // Verifications that throw Exceptions to the Constructor if not validated.
    verifyParentNotNull();

    if (nullptr == _runParams)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "A valid RunParameters must be provided to an Algorithm constructor.");
    }

    if (nullptr == _pbParams)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "A valid PbParameters must be provided to the Algorithm constructor.");
    }

    if ( nullptr == _stopReasons )
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Valid stop reasons must be provided to the Algorithm constructor.");

    // Instanciate generic algorithm termination
    _termination    = std::make_unique<NOMAD::Termination>( this );


    /** Step::userInterrupt() will be called if CTRL-C is pressed.
     * Currently, the master thread will wait for all evaluations to be complete.
     * \todo Propage interruption to all threads, for all parallel evaluations of blackbox.
     */
    signal(SIGINT, userInterrupt);

}


NOMAD::Algorithm::~Algorithm()
{
}


void NOMAD::Algorithm::startImp()
{

    // Comment to appear at the end of stats lines
    // By default, nothing is added
    NOMAD::MainStep::setAlgoComment("");

    // All stop reasons are reset.
    _stopReasons->setStarted();

    if (isMainAlgo())
    {
        // Update hot restart info
        readInformationForHotRestart();
    }

    // By default reset the lap counter for BbEval and set the lap maxBbEval to INF
    NOMAD::EvcInterface::getEvaluatorControl()->resetLapBbEval();
    NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( NOMAD::INF_SIZE_T );
    NOMAD::EvcInterface::getEvaluatorControl()->resetSgteEval();

    if (nullptr == _megaIteration)
    {
        // Default behavior - not hot restart.
        // Clear cache hits.
        // Initialization.
        // Eval X0s.

        // Ensure we do not count cache hits which may have been read in the cache.
        NOMAD::CacheBase::resetNbCacheHits();

        _initialization->start();
        _initialization->run();
        _initialization->end();

    }
    else
    {
        // Hot restart situation.
        // We will not need Initialization.
        auto barrier = _megaIteration->getBarrier();

        // Update X0s
        // Use best feasible, or best infeasible points.
        auto bestFeasPoints = barrier->getAllXFeas();
        auto bestInfPoints  = barrier->getAllXInf();

        NOMAD::ArrayOfPoint x0s;
        if (!bestFeasPoints.empty())
        {
            std::transform(bestFeasPoints.begin(), bestFeasPoints.end(), std::back_inserter(x0s),
                       [](NOMAD::EvalPointPtr evalPointPtr) -> NOMAD::EvalPoint { return *evalPointPtr; });

        }
        else if (!bestInfPoints.empty())
        {
            std::transform(bestInfPoints.begin(), bestInfPoints.end(), std::back_inserter(x0s),
                       [](NOMAD::EvalPointPtr evalPointPtr) -> NOMAD::EvalPoint { return *evalPointPtr; });
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
        displayEvalCounts();
    }

    // By default reset the lap counter for BbEval and set the lap maxBbEval to INF
    NOMAD::EvcInterface::getEvaluatorControl()->resetLapBbEval();
    NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( NOMAD::INF_SIZE_T );

    if (isMainAlgo())
    {
        saveInformationForHotRestart();
    }

    // Reset stats comment 
    NOMAD::MainStep::resetPreviousAlgoComment();
}


void NOMAD::Algorithm::hotRestartOnUserInterrupt()
{
    hotRestartBeginHelper();

    // TODO: To investigate: Reset mesh because parameters have changed.
    // - Maybe already done as a side effect of calling parent step and
    //   resetting parameters.
    // - Maybe should be done at another level, Iteration or MegaIteration.
    // - Useful code:
    /*
    std::stringstream ss;
    auto mesh = getIterationMesh();
    ss << *mesh;
    // Reset pointer
    mesh.reset();
    mesh = std::make_shared<NOMAD::GMesh>(_pbParams);
    // Get old mesh values
    ss >> *mesh;
    */

    hotRestartEndHelper();
}


void NOMAD::Algorithm::saveInformationForHotRestart() const
{
    // If we want to stop completely and then be able to restart
    // from where we were, we need to save some information on file.
    //
    // TODO: Maybe we need to write current parameters. If we write them,
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
    NOMAD::OutputInfo displaySolFeas(_name, sFeas, outputLevel);

    sFeas = "Best feasible solution";
    size_t nbBestFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointList,
                                    getSubFixedVariable(), getEvalType());
    if (0 == nbBestFeas)
    {
        sFeas += ":     Undefined.";
        displaySolFeas.addMsg(sFeas);
    }
    else if (1 == nbBestFeas)
    {
        sFeas += ":     ";
        displaySolFeas.addMsgAndSol(sFeas, *evalPointList.begin());
    }
    else
    {
        sFeas += "s:    ";
        displaySolFeas.addMsgAndSol(sFeas, *evalPointList.begin());
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
            displaySolFeas.addMsgAndSol("                            ",*it);
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
    auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
    if (nullptr != _megaIteration)
    {
        hMax = _megaIteration->getBarrier()->getHMax();
    }
    NOMAD::OutputInfo displaySolInf(_name, sInf, outputLevel);
    sInf = "Best infeasible solution";
    size_t nbBestInf = NOMAD::CacheBase::getInstance()->findBestInf(evalPointList,
                                hMax, getSubFixedVariable(), getEvalType());
    if (0 == nbBestInf)
    {
        sInf += ":   Undefined.";
        displaySolInf.addMsg(sInf);
    }
    else if (1 == nbBestInf)
    {
        sInf += ":   ";
        displaySolInf.addMsgAndSol(sInf, *evalPointList.begin());
    }
    else
    {
        sInf += "s:  ";
        displaySolInf.addMsgAndSol(sInf, *evalPointList.begin());
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
            displaySolInf.addMsgAndSol("                            ",(*it));
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
    size_t sgteEval     = NOMAD::EvcInterface::getEvaluatorControl()->getSgteEval();
    size_t totalSgteEval = NOMAD::EvcInterface::getEvaluatorControl()->getTotalSgteEval();
    size_t nbCacheHits  = NOMAD::CacheBase::getNbCacheHits();
    int nbEvalNoCount   = static_cast<int>(nbEval - bbEval - nbCacheHits);

    // What needs to be shown, according to the counts and to the value of isSub
    bool showNbEvalNoCount  = (nbEvalNoCount > 0);
    bool showSgteEval       = isSub && (sgteEval > 0);
    bool showTotalSgteEval  = (totalSgteEval > 0);
    bool showNbCacheHits    = (nbCacheHits > 0);
    bool showNbEval         = (nbEval > bbEval);
    bool showLapBbEval      = isSub && (bbEval > lapBbEval && lapBbEval > 0);

    // Padding for nice presentation
    std::string sFeedBbEval, sFeedLapBbEval, sFeedNbEvalNoCount, sFeedSgteEval,
                sFeedTotalSgteEval, sFeedCacheHits, sFeedNbEval;

    // Conditional values: showNbEval, showNbEvalNoCount, showLapBbEval
    if (showLapBbEval)  // Longest title
    {
        sFeedBbEval += "                 ";
        //sFeedLapBbEval += "";
        sFeedNbEvalNoCount += "   ";
        sFeedSgteEval += "                     ";
        sFeedTotalSgteEval += "               ";
        sFeedCacheHits += "                           ";
        sFeedNbEval += "          ";
    }
    else if (showNbEvalNoCount) // Second longest
    {
        sFeedBbEval += "              ";
        //sFeedLapBbEval += "";
        //sFeedNbEvalNoCount += "";
        sFeedSgteEval += "                  ";
        sFeedTotalSgteEval += "            ";
        sFeedCacheHits += "                        ";
        sFeedNbEval += "       ";
    }
    else if (showNbEval)    // 3rd longest title
    {
        sFeedBbEval += "       ";
        //sFeedLapBbEval += "";
        //sFeedNbEvalNoCount += "";
        sFeedSgteEval += "           ";
        sFeedTotalSgteEval += "     ";
        sFeedCacheHits += "                 ";
        //sFeedNbEval += "";
    }

    std::string sBbEval         = "Blackbox evaluations: " + sFeedBbEval + NOMAD::itos(bbEval);
    std::string sLapBbEval      = "Sub-optimization blackbox evaluations: " + sFeedLapBbEval + NOMAD::itos(lapBbEval);
    std::string sNbEvalNoCount  = "Blackbox evaluation (not counting): " + sFeedNbEvalNoCount + NOMAD::itos(nbEvalNoCount);
    std::string sSgteEval       = "Sgte evaluations: " + sFeedSgteEval + NOMAD::itos(sgteEval);
    std::string sTotalSgteEval  = "Total sgte evaluations: " + sFeedTotalSgteEval + NOMAD::itos(totalSgteEval);
    std::string sCacheHits      = "Cache hits: " + sFeedCacheHits + NOMAD::itos(nbCacheHits);
    std::string sNbEval         = "Total number of evaluations: " + sFeedNbEval + NOMAD::itos(nbEval);
    
    // Output levels will be modulated depending on the counts and on the Algorithm level.
    NOMAD::OutputLevel outputLevelHigh = isSub ? NOMAD::OutputLevel::LEVEL_INFO
                                               : NOMAD::OutputLevel::LEVEL_HIGH;
    NOMAD::OutputLevel outputLevelNormal = isSub ? NOMAD::OutputLevel::LEVEL_INFO
                                                 : NOMAD::OutputLevel::LEVEL_NORMAL;

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
    if (showSgteEval)
    {
        AddOutputInfo(sSgteEval, outputLevelNormal);
    }
    if (showTotalSgteEval)
    {
        AddOutputInfo(sTotalSgteEval, outputLevelNormal);
    }
    if (showNbCacheHits)
    {
        AddOutputInfo(sCacheHits, outputLevelNormal);
    }
    if (showNbEval)
    {
        AddOutputInfo(sNbEval, outputLevelNormal);
    }
}


bool NOMAD::Algorithm::isSubAlgo() const
{
    bool isSub = false;

    // Get Parent algorithm. No need to cast: We only want to know if it exists.
    const NOMAD::Step* parentAlgo = getParentOfType<NOMAD::Algorithm*>();
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
