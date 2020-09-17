
#include <algorithm>    // For std::merge and std::unique

#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/Search.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Output/OutputQueue.hpp"

#ifdef TIME_STATS
#include "../../Algos/EvcInterface.hpp"
#include "../../Util/Clock.hpp"

// Initialize static variables
double NOMAD::MadsIteration::_iterTime = 0.0;
double NOMAD::MadsIteration::_searchTime = 0.0;
double NOMAD::MadsIteration::_searchEvalTime = 0.0;
double NOMAD::MadsIteration::_pollTime = 0.0;
double NOMAD::MadsIteration::_pollEvalTime = 0.0;
#endif // TIME_STATS


void NOMAD::MadsIteration::init()
{
    _name = NOMAD::Iteration::getName();
}


void NOMAD::MadsIteration::startImp()
{
#ifdef TIME_STATS
    _iterStartTime = NOMAD::Clock::getCPUTime();
#endif // TIME_STATS
}


bool NOMAD::MadsIteration::runImp()
{
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    bool iterationSuccess = false;
    NOMAD::SuccessType bestSuccessYet = NOMAD::SuccessType::NOT_EVALUATED;

    // Parameter Update is handled at the upper level - MegaIteration.

    // 1. Search
    if ( ! _stopReasons->checkTerminate() )
    {
#ifdef TIME_STATS
        double searchStartTime = NOMAD::Clock::getCPUTime();
        double searchEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS
        NOMAD::Search search(this );
        search.start();
        iterationSuccess = search.run();

        NOMAD::SuccessType success = search.getSuccessType();
        if (success > bestSuccessYet)
        {
            bestSuccessYet = success;
        }
        search.end();
#ifdef TIME_STATS
        _searchTime += NOMAD::Clock::getCPUTime() - searchStartTime;
        _searchEvalTime += NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime() - searchEvalStartTime;
#endif // TIME_STATS

    }

    if ( ! _stopReasons->checkTerminate() )
    {
        if (iterationSuccess)
        {
            OUTPUT_INFO_START
            AddOutputInfo("Search Successful. Enlarge Delta frame size.");
            OUTPUT_INFO_END
        }
        else
        {
#ifdef TIME_STATS
            double pollStartTime = NOMAD::Clock::getCPUTime();
            double pollEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS
            // 2. Poll
            NOMAD::Poll poll( this );
            poll.start();
            // Iteration is a success if either a better xFeas or
            // a better xInf (partial success or dominating) xInf was found.
            // See Algorithm 12.2 from DFBO.
            iterationSuccess = poll.run();

            NOMAD::SuccessType success = poll.getSuccessType();
            if (success > bestSuccessYet)
            {
                bestSuccessYet = success;
            }
            poll.end();
#ifdef TIME_STATS
            _pollTime += NOMAD::Clock::getCPUTime() - pollStartTime;
            _pollEvalTime += NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime() - pollEvalStartTime;
#endif // TIME_STATS
        }
    }

    setSuccessType(bestSuccessYet);

    // End of the iteration: iterationSuccess is true if we have a success (partial or full).
    return iterationSuccess;
}


#ifdef TIME_STATS
void NOMAD::MadsIteration::endImp()
{
    _iterTime += NOMAD::Clock::getCPUTime() - _iterStartTime;
}
#endif // TIME_STATS


bool NOMAD::MadsIteration::isMainIteration() const
{
    // This MadsIteration is the main iteration if it has the same mesh and k
    // as its parent MadsMegaIteration, and if the poll center is the first point of the MadsMegaIteration's barrier.
    bool ret = false;

    auto megaIter = getParentOfType<NOMAD::MadsMegaIteration*>();
    if (nullptr != megaIter)
    {
        ret = (megaIter->getMesh() == _mesh && megaIter->getK() == _k);

        if (ret)
        {
            auto firstBarrierPoint = megaIter->getBarrier()->getFirstXFeas();
            if (nullptr == firstBarrierPoint)
            {
                firstBarrierPoint = megaIter->getBarrier()->getFirstXInf();
            }
            ret = (*_frameCenter == *firstBarrierPoint);
        }
    }

    return ret;
}

