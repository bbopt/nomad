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

#include "../Cache/CacheBase.hpp"
#include "../Eval/EvaluatorControl.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Util/AllStopReasons.hpp"
#include "../Util/Clock.hpp"
#include <unistd.h> // For usleep

/*------------------------*/
/* Class EvaluatorControl */
/*------------------------*/
// Initialize EvaluatorControl class.
// To be called by the Constructor.
void NOMAD::EvaluatorControl::init(std::shared_ptr<NOMAD::Evaluator> evaluator,
                                   const std::shared_ptr<NOMAD::EvaluatorControlParameters>& evalContParams)
{
#ifdef _OPENMP
    omp_init_lock(&_evalQueueLock);
#endif // _OPENMP

    auto evalStopReason = std::make_shared<NOMAD::StopReason<NOMAD::EvalStopType>>();
    addMainThread(NOMAD::getThreadNum(), evalStopReason, evaluator, evalContParams);
}


// Terminate EvaluatorControl class.
// To be called by the Destructor.
void NOMAD::EvaluatorControl::destroy()
{
    if (!_evalPointQueue.empty())
    {
        // Show warnings and debug info.
        // Do not scare the user if display degree is medium or low.
        OUTPUT_INFO_START
        std::cerr << "Warning: deleting EvaluatorControl with EvalPoints remaining." << std::endl;
        OUTPUT_INFO_END
        bool showDebug = false;
        OUTPUT_DEBUG_START
        showDebug = true;
        OUTPUT_DEBUG_END
        clearQueue(-1, false /* waitRunning */, showDebug);
    }

    for (int mainThreadNum : _mainThreads)
    {
        if (remainsEvaluatedPoints(mainThreadNum))
        {
            OUTPUT_INFO_START
            std::cerr << "Warning: deleting EvaluatorControl with evaluated points remaining." << std::endl;
            OUTPUT_INFO_END
            OUTPUT_DEBUG_START
            // retrieveAllEvaluatedPoints will also clear _evaluatedPoints.
            for (auto evalPoint : retrieveAllEvaluatedPoints(mainThreadNum))
            {
                std::string s = "Delete evaluated point: ";
                s += evalPoint.display();
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            }
        }
        OUTPUT_DEBUG_END
    }

#ifdef _OPENMP
    omp_destroy_lock(&_evalQueueLock);
#endif // _OPENMP
}


std::shared_ptr<NOMAD::Evaluator> NOMAD::EvaluatorControl::setEvaluator(std::shared_ptr<NOMAD::Evaluator> evaluator)
{
    return getMainThreadInfo().setEvaluator(evaluator);
}


void NOMAD::EvaluatorControl::addMainThread(const int threadNum,
                                            const std::shared_ptr<NOMAD::StopReason<NOMAD::EvalStopType>>,
                                            const std::shared_ptr<NOMAD::Evaluator>& evaluator,
                                            const std::shared_ptr<NOMAD::EvaluatorControlParameters>& evalContParams)
{
    if (!isMainThread(threadNum))
    {
        OUTPUT_DEBUG_START
        std::string s = "Add main thread: " + NOMAD::itos(threadNum);
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        _mainThreads.insert(threadNum);
        // Since EvcMainThreadInfo has atomic members, we have to use map::emplace and the only way I found to 
        // make it work was using this convulated formulation.
        _mainThreadInfo.emplace(std::piecewise_construct, std::forward_as_tuple(threadNum), std::forward_as_tuple(evaluator, evalContParams));
    }
}


NOMAD::EvcMainThreadInfo& NOMAD::EvaluatorControl::getMainThreadInfo(const int threadNum) const
{
    int mainThreadNum = threadNum;
    if (-1 == mainThreadNum)
    {
        mainThreadNum = NOMAD::getThreadNum();
    }
    if (!isMainThread(mainThreadNum))
    {
        std::string s = "Thread " + NOMAD::itos(threadNum) + " is not a main thread";
        NOMAD::Exception(__FILE__,__LINE__,s);
    }

    return _mainThreadInfo.at(mainThreadNum);
}


size_t NOMAD::EvaluatorControl::getSgteEval(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getSgteEval();
}


void NOMAD::EvaluatorControl::incSgteEval(const int threadNum, const size_t countEval)
{
    getMainThreadInfo(threadNum).incSgteEval(countEval);
}


void NOMAD::EvaluatorControl::resetSgteEval()
{
    getMainThreadInfo().resetSgteEval();
}


size_t NOMAD::EvaluatorControl::getBbEvalInSubproblem() const
{
    return getMainThreadInfo().getBbEvalInSubproblem();
}


void NOMAD::EvaluatorControl::incBbEvalInSubproblem(const size_t countEval)
{
    getMainThreadInfo().incBbEvalInSubproblem(countEval);
}


void NOMAD::EvaluatorControl::resetBbEvalInSubproblem()
{
    getMainThreadInfo().resetBbEvalInSubproblem();
}


// getNbEval() and setNbEval interface with Mads.
// The number of cache hits must be added to getNbEval(), and removed when
// setting _nbEvalSentToEvaluator.
size_t NOMAD::EvaluatorControl::getNbEval() const
{
    return _nbEvalSentToEvaluator + NOMAD::CacheBase::getNbCacheHits();
}


void NOMAD::EvaluatorControl::setNbEval(const size_t nbEval)
{
    if (nbEval < NOMAD::CacheBase::getNbCacheHits())
    {
        std::cerr << "Warning: trying to set EvaluatorControl NbEval to negative value: " << nbEval << " - " << NOMAD::CacheBase::getNbCacheHits() << std::endl;
    }
    else
    {
        _nbEvalSentToEvaluator = nbEval - NOMAD::CacheBase::getNbCacheHits();
    }
}


size_t NOMAD::EvaluatorControl::getLapMaxBbEval() const
{
    return getMainThreadInfo().getLapMaxBbEval();
}


void NOMAD::EvaluatorControl::setLapMaxBbEval(const size_t maxBbEval)
{
    getMainThreadInfo().setLapMaxBbEval(maxBbEval);
}


void NOMAD::EvaluatorControl::resetLapBbEval()
{
    getMainThreadInfo().resetLapBbEval();
}


size_t NOMAD::EvaluatorControl::getLapBbEval(const int mainThreadNum) const
{
    return getMainThreadInfo(mainThreadNum).getLapBbEval();
}


void NOMAD::EvaluatorControl::incLapBbEval(const size_t countEval)
{
    getMainThreadInfo().incLapBbEval(countEval);
}


size_t NOMAD::EvaluatorControl::getQueueSize(const int threadNum) const
{
    if (-1 == threadNum)
    {
        return _evalPointQueue.size();
    }

    return std::count_if(_evalPointQueue.begin(), _evalPointQueue.end(),
                         [threadNum](const std::shared_ptr<EvalQueuePoint>& evalQueuePoint){return threadNum == evalQueuePoint->getThreadAlgo();});
}


bool NOMAD::EvaluatorControl::getDoneWithEval(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getDoneWithEval();
}


void NOMAD::EvaluatorControl::setDoneWithEval(const int threadNum, const bool doneWithEval)
{
    getMainThreadInfo(threadNum).setDoneWithEval(doneWithEval);
}


void NOMAD::EvaluatorControl::setBarrier(const std::shared_ptr<NOMAD::Barrier>& barrier)
{
    getMainThreadInfo().setBarrier(barrier);
}


const std::shared_ptr<NOMAD::Barrier>& NOMAD::EvaluatorControl::getBarrier(const int mainThreadNum) const
{
    return getMainThreadInfo(mainThreadNum).getBarrier();
}


void NOMAD::EvaluatorControl::setComputeSuccessTypeFunction(const ComputeSuccessFunction& computeSuccessFunction)
{
    getMainThreadInfo().setComputeSuccessTypeFunction(computeSuccessFunction);
}


void NOMAD::EvaluatorControl::setStopReason(const int threadNum, const EvalStopType& s)
{
    getMainThreadInfo(threadNum).setStopReason(s);
}


const NOMAD::StopReason<NOMAD::EvalStopType>& NOMAD::EvaluatorControl::getStopReason(const int mainThreadNum) const
{
    return getMainThreadInfo(mainThreadNum).getStopReason();
}


std::string NOMAD::EvaluatorControl::getStopReasonAsString(const int mainThreadNum) const
{
    return getMainThreadInfo(mainThreadNum).getStopReasonAsString();
}


bool NOMAD::EvaluatorControl::testIf(const EvalStopType& s) const
{
    return getMainThreadInfo().testIf(s);
}


bool NOMAD::EvaluatorControl::checkEvalTerminate() const
{
    return getMainThreadInfo().checkEvalTerminate();
}


std::vector<NOMAD::EvalPoint> NOMAD::EvaluatorControl::retrieveAllEvaluatedPoints(const int threadNum)
{
    return getMainThreadInfo(threadNum).retrieveAllEvaluatedPoints();
}


void NOMAD::EvaluatorControl::addEvaluatedPoint(const int threadNum, const EvalPoint& evaluatedPoint)
{
    getMainThreadInfo(threadNum).addEvaluatedPoint(evaluatedPoint);
}


bool NOMAD::EvaluatorControl::remainsEvaluatedPoints(const int threadNum) const
{
    return getMainThreadInfo(threadNum).remainsEvaluatedPoints();
}


void NOMAD::EvaluatorControl::clearEvaluatedPoints(const int threadNum)
{
    getMainThreadInfo(threadNum).clearEvaluatedPoints();
}


const NOMAD::SuccessType& NOMAD::EvaluatorControl::getSuccessType(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getSuccessType();
}


void NOMAD::EvaluatorControl::setSuccessType(const int threadNum, const NOMAD::SuccessType& success)
{
    getMainThreadInfo(threadNum).setSuccessType(success);
}


size_t NOMAD::EvaluatorControl::getCurrentlyRunning(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getCurrentlyRunning();
}


void NOMAD::EvaluatorControl::incCurrentlyRunning(const int threadNum, const size_t n)
{
    getMainThreadInfo(threadNum).incCurrentlyRunning(n);
}


void NOMAD::EvaluatorControl::decCurrentlyRunning(const int threadNum, const size_t n)
{
    getMainThreadInfo(threadNum).decCurrentlyRunning(n);
}


NOMAD::Double NOMAD::EvaluatorControl::getHMax(const int threadNum) const
{
    NOMAD::Double hMax = NOMAD::INF;
    auto barrier = getBarrier(threadNum);
    if (nullptr != barrier)
    {
        hMax = barrier->getHMax();
    }

    return hMax;
}


std::shared_ptr<NOMAD::EvalParameters> NOMAD::EvaluatorControl::getEvalParams(const int threadNum) const
{
    return getMainThreadInfo().getEvalParams();
}


bool NOMAD::EvaluatorControl::getOpportunisticEval(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getOpportunisticEval();
}


void NOMAD::EvaluatorControl::setOpportunisticEval(const bool opportunisticEval)
{
    getMainThreadInfo().setOpportunisticEval(opportunisticEval);
}


bool NOMAD::EvaluatorControl::getUseCache(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getUseCache();
}


void NOMAD::EvaluatorControl::setUseCache(const bool useCache)
{
    getMainThreadInfo().setUseCache(useCache);
}


NOMAD::EvalType NOMAD::EvaluatorControl::getEvalType(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getEvalType();
}


size_t NOMAD::EvaluatorControl::getMaxBbEvalInSubproblem(const int threadNum) const
{
    return getMainThreadInfo(threadNum).getMaxBbEvalInSubproblem();
}


void NOMAD::EvaluatorControl::lockQueue()
{
#ifdef _OPENMP
    // Sanity check before locking the queue.
    // Verify we are in a main thread.
    const int threadNum = NOMAD::getThreadNum();
    if (!isMainThread(threadNum))
    {
        std::string err = "Error: EvaluatorControl::lockQueue called from thread ";
        err += std::to_string(threadNum);
        err += ", which is not a main thread.";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    // Note: The queue could be already locked, ex. by popEvalPoint.
    omp_set_lock(&_evalQueueLock);
#endif // _OPENMP
}


void NOMAD::EvaluatorControl::unlockQueue(bool doSort)
{
#ifdef _OPENMP
    // Sanity checks before unlocking the queue
    // 1- Verify we are in a main thread.
    const int threadNum = NOMAD::getThreadNum();
    if (!isMainThread(threadNum))
    {
        std::string err = "Error: EvaluatorControl::unlockQueue called from thread ";
        err += std::to_string(threadNum);
        err += ", which is not a main thread.";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    // 2- Verify the queue was already locked. The lock should have been set by lockQueue().
    if (omp_test_lock(&_evalQueueLock))
    {
        // Queue was not locked. Queue is now locked.
        std::string err = "Error: tring to unlock a queue that was not locked.";
        omp_unset_lock(&_evalQueueLock);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    // Now, unlock the queue.
    omp_unset_lock(&_evalQueueLock);
#endif // _OPENMP

    // The EvalQueuePoints are added randomly.
    // Sort the queue, using default sort, if doSort is true (default).
    // In non-opportunistic context, it is useless to sort.
    if (doSort && getOpportunisticEval() && getQueueSize() > 1)
    {
        sort(_comp);
    }
}


// Add an EvalPoint to the Queue
bool NOMAD::EvaluatorControl::addToQueue(const NOMAD::EvalQueuePointPtr &evalQueuePoint)
{
    bool pointInserted = false;

    if (!evalQueuePoint->ArrayOfDouble::isComplete())
    {
        std::string err("EvaluatorControl: addToQueue: Adding an undefined Point for evaluation: ");
        err += evalQueuePoint->getX()->NOMAD::Point::display();
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
#ifdef _OPENMP
    // Thread-safety necessary here.
    // Ensure we will keep adding points until we ask to eval the queue, using run().
    // I.e. Ensure the queue is already locked.
    if (omp_test_lock(&_evalQueueLock))
    {
        std::string err = "Error: tring to add an element to a queue that was not locked.";
        // Unlock queue before throwing exception.
        // If we are in this section, it means that the queue was locked by
        // the call to omp_test_lock().
        omp_unset_lock(&_evalQueueLock);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
#endif // _OPENMP

    // Insert at the beginning of the queue. The points will be sorted when the
    // queue is unlocked.
    // If points are not sorted, they remain in reverse lexicographical order.
    NOMAD::EvalPoint foundEvalPoint;
    auto evalType = evalQueuePoint->getEvalType();
    if (std::find_if(_evalPointQueue.begin(), _evalPointQueue.end(), [evalQueuePoint](NOMAD::EvalQueuePointPtr eqp){ return (*eqp == *evalQueuePoint); }) != _evalPointQueue.end())
    {
        // Point is already in queue, do not insert it again.
    }
    else if (NOMAD::CacheBase::getInstance()->find(*evalQueuePoint, foundEvalPoint)
             && NOMAD::EvalStatusType::EVAL_IN_PROGRESS == foundEvalPoint.getEvalStatus(evalType))
    {
         OUTPUT_DEBUG_START
         std::string s = "Evaluation is already in progress for point: " + foundEvalPoint.displayAll();
         NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
         OUTPUT_DEBUG_END
    }
    else
    {
        pointInserted = (_evalPointQueue.end() != _evalPointQueue.insert(_evalPointQueue.begin(), evalQueuePoint));
    }

    return pointInserted;
}


// Get the top EvalPoint from the Queue and pop it
// Return true if it worked, false if it failed.
bool NOMAD::EvaluatorControl::popEvalPoint(NOMAD::EvalQueuePointPtr &evalQueuePoint,
                                           int mainThreadNum)
{
    bool success = false;
#ifdef _OPENMP
    omp_set_lock(&_evalQueueLock);
#endif // _OPENMP
    if (!_evalPointQueue.empty())
    {
        if (-1 == mainThreadNum)
        {
            // Remove last element, simulate a "pop".
            evalQueuePoint = std::move(_evalPointQueue[_evalPointQueue.size()-1]);
            _evalPointQueue.erase(_evalPointQueue.end()-1);
            mainThreadNum = evalQueuePoint->getThreadAlgo();
            success = true;
        }
        else
        {
            // Find the last element with the given thread number and pop it
            // Done in a naive way; this part could be made nicer.
            for (int i = static_cast<int>(_evalPointQueue.size())-1; i >= 0; i--)
            {
                if (mainThreadNum == _evalPointQueue[i]->getThreadAlgo())
                {
                    evalQueuePoint = std::move(_evalPointQueue[i]);
                    _evalPointQueue.erase(_evalPointQueue.begin()+i);
                    success = true;
                    break;
                }
            }
        }
        // Count this evaluation as running as soon as it is popped from the queue.
        if (success)
        {
            incCurrentlyRunning(mainThreadNum, 1);
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_evalQueueLock);
#endif // _OPENMP

    return success;
}


bool NOMAD::EvaluatorControl::popBlock(NOMAD::BlockForEval &block)
{
    bool success = false;
    bool popWorks = true;
    size_t bbBlockSize = NOMAD::INF_SIZE_T;
    size_t sgteBlockSize = NOMAD::INF_SIZE_T;
    size_t blockSize = 1;
    bool gotBlockSize = false;


    // NOTE: Some racing conditions may happen with value BB_MAX_BLOCK_SIZE.
    // As a workaround, try to get the attribute until no exception is thrown.
    while (!gotBlockSize)
    {
        try
        {
            bbBlockSize = _evalContGlobalParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
            sgteBlockSize = _evalContGlobalParams->getAttributeValue<size_t>("SGTE_MAX_BLOCK_SIZE");
            gotBlockSize = true;
        }
        catch (NOMAD::ParameterToBeChecked &e)
        {
            // While will loop - Retry
        }
    }


    // Need to pop a block of points which will all be evaluated using the
    // same evaluator.
    // A block of points contains points that were all generated by the same main thread.
    int mainThreadNum = -1;
    bool firstPop = true;
    while (block.size() < blockSize && popWorks)
    {
        NOMAD::EvalQueuePointPtr evalQueuePoint;
        popWorks = popEvalPoint(evalQueuePoint, mainThreadNum);
        if (popWorks)
        {
            mainThreadNum = evalQueuePoint->getThreadAlgo();
            block.push_back(std::move(evalQueuePoint));
            success = true;
            if (firstPop)
            {
                // Update blockSize depending on eval type
                auto blockEvalType = getEvalType(mainThreadNum);
                switch (blockEvalType)
                {
                    case NOMAD::EvalType::BB:
                        blockSize = bbBlockSize;
                        break;
                    case NOMAD::EvalType::SGTE:
                        blockSize = sgteBlockSize;
                        break;
                    default:
                        std::cerr << "EvaluatorControl::popBlock: Unknown eval type " << blockEvalType << std::endl;
                }
                firstPop = false;
            }
        }
    }

    return success;
}


void NOMAD::EvaluatorControl::sort(NOMAD::ComparePriority comp)
{
#ifdef _OPENMP
    omp_set_lock(&_evalQueueLock);
#endif // _OPENMP

    std::string s;
    OUTPUT_DEBUG_START
    s = "Evaluation queue before sort:";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    for (auto evalPoint : _evalPointQueue)
    {
        s = "\t" + evalPoint->display();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    OUTPUT_DEBUG_END

    std::sort(_evalPointQueue.begin(), _evalPointQueue.end(), comp);

    OUTPUT_DEBUG_START
    s = "Evaluation queue after sort:";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    for (auto evalPoint : _evalPointQueue)
    {
        s = "\t" + evalPoint->display();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    OUTPUT_DEBUG_END

#ifdef _OPENMP
    omp_unset_lock(&_evalQueueLock);
#endif // _OPENMP
}


void NOMAD::EvaluatorControl::clearQueue(int mainThreadNum, const bool waitRunning, const bool showDebug)
{
    if (-1 == mainThreadNum)
    {
        mainThreadNum = NOMAD::getThreadNum();
    }

    // Wait for any currently running evaluations to be done.
    // Otherwise, the queue would be empty but old evaluations
    // could pop up unexpectedtly.
    if (waitRunning && mainThreadNum >= 0)
    {
        bool warningShown = false;
        while (getCurrentlyRunning(mainThreadNum) > 0)
        {
            OUTPUT_INFO_START
            if (!warningShown)
            {
                std::string s = "Clear queue for main thread " + NOMAD::itos(mainThreadNum);
                s += ": Waiting for " + NOMAD::itos(getCurrentlyRunning(mainThreadNum));
                s += " evaluations to complete.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
                warningShown = true;
            }
            OUTPUT_INFO_END
            usleep(10);
        }
    }

#ifdef _OPENMP
    omp_set_lock(&_evalQueueLock);
#endif // _OPENMP

    auto it = std::remove_if(_evalPointQueue.begin(),
                             _evalPointQueue.end(),
                             [mainThreadNum, showDebug](const NOMAD::EvalQueuePointPtr& evalQueuePoint)
                             {
                                 OUTPUT_DEBUG_START
                                 if (showDebug)
                                 {
                                     std::string s = "Delete point from queue: ";
                                     s += evalQueuePoint->display();
                                     NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                                 }
                                 OUTPUT_DEBUG_END
                                 return (mainThreadNum < 0 || mainThreadNum == evalQueuePoint->getThreadAlgo());
                             });
    _evalPointQueue.erase(it, _evalPointQueue.end());

#ifdef _OPENMP
    omp_unset_lock(&_evalQueueLock);
#endif // _OPENMP
}


// Evaluate all points in the queue, or stop under some conditions.
// If strategy is opportunistic (parameter OPPORTUNISTIC_EVAL), stop
// as soon as a successful point is found, and flush the queue.
//
// Points must already be in the cache.
NOMAD::SuccessType NOMAD::EvaluatorControl::run()
{
    const int threadNum = NOMAD::getThreadNum();
    const bool inMainThread = isMainThread(threadNum);
    // Main threads only:
    //  - Reset success
    //  - Set Barrier and flag OpportunisticEval
    //  - Print info
    if (inMainThread)
    {
        // At this point, the threads other than the current main thread might have already
        // started evaluating.

        setSuccessType(threadNum, NOMAD::SuccessType::UNSUCCESSFUL);
        // Transfer stop reason from AllStopReasons to EvcMainThreadInfo.
        setStopReason(threadNum, NOMAD::AllStopReasons::getEvalStopReason().get());

        // An empty eval queue must be accounted for
        if (_evalPointQueue.empty())
        {
            setStopReason(threadNum, NOMAD::EvalStopType::EMPTY_LIST_OF_POINTS);
        }

        // Update stop reason.
        if (   checkEvalTerminate()
            || NOMAD::AllStopReasons::checkBaseTerminate())
        {
            OUTPUT_DEBUG_START
            std::string sStopReason = "EvaluatorControl stop reason: ";
            sStopReason += (NOMAD::AllStopReasons::checkBaseTerminate())
                ? NOMAD::AllStopReasons::getBaseStopReasonAsString() + " (Base) "
                :  getStopReasonAsString(threadNum) + " (Eval) ";

            NOMAD::OutputQueue::Add(sStopReason, NOMAD::OutputLevel::LEVEL_DEBUG);
            NOMAD::OutputQueue::Flush();
            OUTPUT_DEBUG_END
        }
        else
        {
            // Ready to run.
            // Reset EvalStopType to STARTED.
            setStopReason(threadNum, NOMAD::EvalStopType::STARTED);
        }

        OUTPUT_DEBUG_START
        std::string s = "Start evaluation.";
#ifdef _OPENMP
        s += " Main thread = ";
        s += std::to_string(threadNum);
        s += ",";
#endif // _OPENMP
        s += " Opportunism = ";
        s += NOMAD::boolToString(getOpportunisticEval(threadNum));
        s += ", UseCache = ";
        s += NOMAD::boolToString(getUseCache(threadNum));
        s += ", Barrier";
        auto barrier = getBarrier(threadNum);
        s += ((nullptr == barrier) ? " = NULL" : ":\n" + barrier->display(4)); // Display a maximum of 4 xFeas and 4 xInf
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
    }

    // Queue runs forever on non-main threads.
    // On main threads, queue runs until stopMainEval() is true.
    bool conditionForStop = false;

    // For debug display
    time_t lastDisplayed = 0;

    // conditionForStop is true if we are in a main thread and stopMainEval() returns true.
    // conditionForStop is true in any thread if reachedMaxEval() returns true; otherwise, it is always false.
    while (!conditionForStop && !_allDoneWithEval)
    {
        // Check for stop conditions
        if (inMainThread)
        {
            conditionForStop = stopMainEval();
        }
        // If we reached max eval, we also stop (valid for all threads).
        conditionForStop = conditionForStop || reachedMaxEval();

        NOMAD::BlockForEval block;
        if (!conditionForStop && popBlock(block))
        {
            // mainThreadNum is the same for all points in the block.
            const int mainThreadNum = block[0]->getThreadAlgo();
            if (evalBlock(block))
            {
                // Update SuccessType
                // success is a member of EvcMainThreadInfo and so it can be shared between secundary threads.
#ifdef _OPENMP
                #pragma omp critical(updateSuccessType)
#endif // _OPENMP
                {
                    for (auto it = block.begin(); it < block.end(); it++)
                    {
                        NOMAD::EvalQueuePointPtr evalQueuePoint = (*it);
                        // Update success type for return
                        if (evalQueuePoint->getSuccess() > getSuccessType(mainThreadNum))
                        {
                            setSuccessType(mainThreadNum, evalQueuePoint->getSuccess());
                            evalQueuePoint->setRelativeSuccess(true);
                        }
                    }
                }   // End critical(updateSuccessType)

                // No need for a critical region here, stopReason would not be set back to NO_STOP
                // by another thread. It could only be set to another stop reason which
                // would also be valid.
                if (getOpportunisticEval(mainThreadNum) && getSuccessType(mainThreadNum) >= NOMAD::SuccessType::PARTIAL_SUCCESS)
                {
                    setStopReason(mainThreadNum, NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS);
                }

                AddStatsInfo(block);
            }
            decCurrentlyRunning(mainThreadNum, block.size());
        }
        else if (!_allDoneWithEval)
        {
            displayDebugWaitingInfo(lastDisplayed);
        }
        else // Queue is empty and we are doneWithEval
        {
            if (inMainThread)
            {
                OUTPUT_DEBUG_START
                std::string s = "Queue is empty and we are done with evaluations.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
            }
            break;
        }

    }   // End of while loop: Exit for this main thread.
        // Other threads keep on looping.


    const bool clearEvalQueue = _evalContGlobalParams->getAttributeValue<bool>("CLEAR_EVAL_QUEUE");
    if (inMainThread)
    {
        // Special case:
        // stopReason is ALL_POINTS_EVALUATED because there are no points left
        // in the queue.
        // But some points are still being evaluated.
        // We must wait for all points to really be evaluated.
        bool warningShown = false;
        // Note that when all points are evaluated, getSuccessType() has the correct
        // value, even if it was modified by those last points being evaluated.
        while ( (clearEvalQueue || testIf(NOMAD::EvalStopType::ALL_POINTS_EVALUATED))
                && getCurrentlyRunning() > 0)
        {
            OUTPUT_INFO_START
            if (!warningShown)
            {
                std::string s = "Waiting for " + NOMAD::itos(getCurrentlyRunning());
                s += " evaluations to complete.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
                warningShown = true;
            }
            OUTPUT_INFO_END
            usleep(10);

            // Update stopReason in case we found a success
            if (getOpportunisticEval(threadNum)
                && getSuccessType(threadNum) >= NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                setStopReason(threadNum, NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS);
            }

            // Update stopReason in case we reached the different max eval criterions
            reachedMaxEval();
            reachedMaxStepEval(threadNum);
        }
    }

    if (inMainThread)
    {
        if (clearEvalQueue)
        {
            // Remove remaining points from queue, to start fresh next time.
            // Otherwise, we keep on evaluating the points remaining in the queue.
            OUTPUT_INFO_START
            NOMAD::OutputQueue::Add("Evaluation is done. Clear queue of " + NOMAD::itos(getQueueSize()) + " points.");
            OUTPUT_INFO_END
            clearQueue(threadNum, true /* waitRunning */, false /* showDebug */);
        }

        // Transfer main thread stop reason from EvcMainThreadInfo to static in
        // AllStopReasons, if it is not already set.
        if (NOMAD::AllStopReasons::testIf(NOMAD::EvalStopType::STARTED))
        {
            NOMAD::AllStopReasons::set(getStopReason(threadNum).get());
        }
        OUTPUT_DEBUG_START
        NOMAD::OutputQueue::Add("EvaluatorControl stop reason: " + NOMAD::AllStopReasons::getEvalStopReasonAsString(), NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        NOMAD::OutputQueue::Flush();
    }

    return (inMainThread) ? getSuccessType(threadNum) : NOMAD::SuccessType::UNSUCCESSFUL;
}


void NOMAD::EvaluatorControl::stop()
{
    std::string s;
    const int mainThreadNum = getThreadNum();
    setDoneWithEval(mainThreadNum, true);
    OUTPUT_DEBUG_START
    s = "Stop evaluation queue for main thread " + std::to_string(mainThreadNum);
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    // Go through all main thread info to see if _allDoneWithEval must be set.
    bool allDone = true;
    for (int mainThreadNum : _mainThreads)
    {
        if (!getDoneWithEval(mainThreadNum))
        {
            // Not done with queue because thread mainThreadNum is not done.
            allDone = false;
            break;
        }
    }

    if (allDone)
    {
        OUTPUT_DEBUG_START
        s = "All main threads are done. Done with evaluation queue.";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        _allDoneWithEval = true;
    }
}


void NOMAD::EvaluatorControl::restart()
{
    _allDoneWithEval = false;
    for (int mainThreadNum : _mainThreads)
    {
        setDoneWithEval(mainThreadNum, false);
    }
}


bool NOMAD::EvaluatorControl::stopMainEval() const
{
    // Check for stop conditions for a main thread.
    // This function should not be called from other threads. But we
    // do not verify the thread number.
    const int mainThreadNum = NOMAD::getThreadNum();

    // Inspect stopReasons
    bool doStopEval = getStopReason(mainThreadNum).checkTerminate();
    doStopEval = doStopEval || NOMAD::AllStopReasons::checkEvalTerminate();

    // Inspect for opportunistic success set in a previous point
    doStopEval = doStopEval || testIf(NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS);

    // Update stopReason if there is no more point to evaluate for this thread
    if (_evalPointQueue.empty() && (!doStopEval || testIf(NOMAD::EvalStopType::EMPTY_LIST_OF_POINTS)))
    {
        getMainThreadInfo().setStopReason(NOMAD::EvalStopType::ALL_POINTS_EVALUATED);
        doStopEval = true;
    }

    // Update on max eval, valid for all threads
    doStopEval = doStopEval || reachedMaxStepEval(mainThreadNum) || reachedMaxEval();

    // Inspect on base stop reason
    bool doStopBase = NOMAD::AllStopReasons::checkBaseTerminate();

    OUTPUT_DEBUG_START
    std::string s = "stopMainEval: return true because: ";
    if (doStopEval)
    {
        s += NOMAD::AllStopReasons::getEvalStopReasonAsString();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }

    if (doStopBase)
    {
        s += NOMAD::AllStopReasons::getBaseStopReasonAsString();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    NOMAD::OutputQueue::Flush();
    OUTPUT_DEBUG_END

    return ( doStopBase || doStopEval ) ;
}


// Eval at most a total maxBbEval, maxEval, or maxBlockEval points in the queue.
// If this condition is true, we assume that it will remain true for the
// rest of the optimization.
// If the condition is true temporary, for example lap or sgte evals,
// use stopMainEval().
bool NOMAD::EvaluatorControl::reachedMaxEval() const
{
    // There is a possiblity that these values are modified
    // in master thread, and that getAttributeValue() gets called before
    // checkAndComply(). In that case, catch the exception and ignore
    // the test.
    bool ret = false;

    if (   NOMAD::AllStopReasons::testIf(NOMAD::EvalStopType::MAX_BB_EVAL_REACHED)
        || NOMAD::AllStopReasons::testIf(NOMAD::EvalStopType::MAX_EVAL_REACHED)
        || NOMAD::AllStopReasons::testIf(NOMAD::EvalStopType::MAX_BLOCK_EVAL_REACHED))
    {
        // Conditions already reached.
        return true;
    }

    size_t maxBbEval    = NOMAD::INF_SIZE_T;
    size_t maxEval      = NOMAD::INF_SIZE_T;
    size_t maxBlockEval = NOMAD::INF_SIZE_T;

    try
    {
        maxBbEval    = _evalContGlobalParams->getAttributeValue<size_t>("MAX_BB_EVAL");
        maxEval      = _evalContGlobalParams->getAttributeValue<size_t>("MAX_EVAL");
        maxBlockEval = _evalContGlobalParams->getAttributeValue<size_t>("MAX_BLOCK_EVAL");
    }
    catch (NOMAD::ParameterToBeChecked &e)
    {
        // Early out. Ignore the test: return false.
        return false;
    }

    std::string s = "Reached stop criterion: ";
    if (maxBbEval < NOMAD::INF_SIZE_T && _bbEval >= maxBbEval)
    {
        // Reached maxBbEval.
        // Note that we have (total number of threads -1) other threads currently
        // running, so the total number of bb evaluations will be over
        // maxBbEval.
        NOMAD::AllStopReasons::set(NOMAD::EvalStopType::MAX_BB_EVAL_REACHED);
        s += NOMAD::AllStopReasons::getEvalStopReasonAsString() + " " + NOMAD::itos(_bbEval);
        ret = true;
    }
    else if (maxEval < NOMAD::INF_SIZE_T && getNbEval() >= maxEval)
    {
        // Reached maxEval.
        NOMAD::AllStopReasons::set(NOMAD::EvalStopType::MAX_EVAL_REACHED);
        s += NOMAD::AllStopReasons::getEvalStopReasonAsString() + " " + NOMAD::itos(getNbEval());
        ret = true;
    }
    else if (maxBlockEval < NOMAD::INF_SIZE_T && _blockEval >= maxBlockEval)
    {
        // Reached maxBlockEval.
        NOMAD::AllStopReasons::set(NOMAD::EvalStopType::MAX_BLOCK_EVAL_REACHED);
        s += NOMAD::AllStopReasons::getEvalStopReasonAsString() + " " + NOMAD::itos(_blockEval);
        ret = true;
    }

#ifdef _OPENMP
    #pragma omp single nowait
#endif // _OPENMP
    {
        if (ret)
        {
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_NORMAL);
        }
    }

    return ret;
}


// Have we reached max eval for a sub step: Number of laps, number of Sgte evals?
bool NOMAD::EvaluatorControl::reachedMaxStepEval(const int mainThreadNum) const
{
    bool ret = false;

    if (   testIf(NOMAD::EvalStopType::MAX_SGTE_EVAL_REACHED)
        || testIf(NOMAD::EvalStopType::LAP_MAX_BB_EVAL_REACHED)
        || testIf(NOMAD::EvalStopType::SUBPROBLEM_MAX_BB_EVAL_REACHED))
    {
        // Conditions already reached.
        return true;
    }

    // Note: Default for NOMAD 3 is 10000. Default for NOMAD 4 is 100.
    // Look if/how we can increment the number of evaluations without
    // increasing the time too much.
    const size_t maxSgteEval = _evalContGlobalParams->getAttributeValue<size_t>("MAX_SGTE_EVAL");
    const size_t maxLapBbEval = getMainThreadInfo(mainThreadNum).getLapMaxBbEval();
    const size_t maxBbEvalInSub = getMaxBbEvalInSubproblem(mainThreadNum);
    std::string s = "Reached sub step stop criterion: ";
    if (maxSgteEval < NOMAD::INF_SIZE_T && getSgteEval(mainThreadNum) >= maxSgteEval)
    {
        // Reached maxSgteEval, max number of eval in Sgte context.
        getMainThreadInfo(mainThreadNum).setStopReason(NOMAD::EvalStopType::MAX_SGTE_EVAL_REACHED);
        s += getStopReasonAsString(mainThreadNum) + " " + NOMAD::itos(getSgteEval(mainThreadNum));
        ret = true;
    }
    else if (maxLapBbEval < NOMAD::INF_SIZE_T && getLapBbEval(mainThreadNum) >= maxLapBbEval)
    {
        // Reached lapMaxEval.
        getMainThreadInfo(mainThreadNum).setStopReason(NOMAD::EvalStopType::LAP_MAX_BB_EVAL_REACHED);
        s += getStopReasonAsString(mainThreadNum) + " " + NOMAD::itos(getLapBbEval(mainThreadNum));
        ret = true;
    }
    else if (maxBbEvalInSub < NOMAD::INF_SIZE_T && getBbEvalInSubproblem() > maxBbEvalInSub)
    {
        // Reached max eval in supbroblem.
        getMainThreadInfo(mainThreadNum).setStopReason(NOMAD::EvalStopType::SUBPROBLEM_MAX_BB_EVAL_REACHED);
        s += getStopReasonAsString(mainThreadNum) + " " + NOMAD::itos(getBbEvalInSubproblem());
        ret = true;
    }

#ifdef _OPENMP
    #pragma omp single nowait
#endif // _OPENMP
    {
        if (ret)
        {
            OUTPUT_DEBUG_START
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
        }
    }

    return ret;
}


void NOMAD::EvaluatorControl::displayDebugWaitingInfo(time_t &lastDisplayed) const
{
#ifdef _OPENMP
    // Check lastDisplayed value againt now, so this message
    // is not continuously displayed.
    time_t now;
    time(&now);
    if (difftime(now, lastDisplayed) >= 1)
    {
        OUTPUT_DEBUG_START
        std::string s = "Thread: " + NOMAD::itos(NOMAD::getThreadNum());
        s += " Waiting for points.";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        lastDisplayed = now;
    }
#endif // _OPENMP
}


void NOMAD::EvaluatorControl::AddStatsInfo(const NOMAD::BlockForEval& block) const
{
    OUTPUT_STATS_START
    for (auto it = block.begin(); it < block.end(); it++)
    {
        NOMAD::EvalQueuePointPtr evalQueuePoint = (*it);
        const int mainThreadNum = evalQueuePoint->getThreadAlgo();

        // SGTE optimizations generate a lot of stats output. Do not show them.
        // Only show BB optimizations.
        if (NOMAD::EvalType::BB != evalQueuePoint->getEvalType())
        {
            return;
        }
        // Evaluation info for output
        NOMAD::StatsInfoUPtr stats(new NOMAD::StatsInfo());
        size_t n = evalQueuePoint->size();

        // Values that are unknown at this point
        NOMAD::ArrayOfDouble unknownMeshIndex(n);
        const int threadNum = NOMAD::getThreadNum();

        // As of February 2019, values for bbEval (BBE) and nbEval (EVAL) make
        // enough sense in threaded environment to be printed out. If at some point
        // they do not, for example if too many iterations end up with the same
        // value, or if the user is dissatisfied with this output, we may update it
        // in a stricter way.

        stats->setObj(evalQueuePoint->getF(NOMAD::EvalType::BB));
        stats->setConsH(evalQueuePoint->getH(NOMAD::EvalType::BB));
        stats->setHMax(getHMax(mainThreadNum));
        stats->setBBE(_bbEval);
        stats->setLap(getLapBbEval(mainThreadNum));
        stats->setSgte(getSgteEval(mainThreadNum));
        stats->setTotalSgte(_totalSgteEval);
        stats->setBlkEva(_blockEval);
        stats->setBlkSize(block.size());
        stats->setBBO(evalQueuePoint->getBBO(NOMAD::EvalType::BB));
        stats->setEval(_nbEvalSentToEvaluator);
        stats->setCacheHits(NOMAD::CacheBase::getNbCacheHits());
        stats->setCacheSize(NOMAD::CacheBase::getInstance()->size());
        stats->setIterNum(evalQueuePoint->getK());
        stats->setTime(NOMAD::Clock::getTimeSinceStart());
        stats->setMeshIndex(unknownMeshIndex);
        stats->setMeshSize(evalQueuePoint->getMeshSize());
        stats->setFrameSize(evalQueuePoint->getFrameSize());
        stats->setSol(*(evalQueuePoint->getX()));
        stats->setThreadAlgo(mainThreadNum);
        stats->setThreadNum(threadNum);
        stats->setRelativeSuccess(evalQueuePoint->getRelativeSuccess());
        stats->setComment(evalQueuePoint->getComment());
        stats->setGenStep(evalQueuePoint->getGenStep());

        std::string s = "Evaluated point: " + evalQueuePoint->display();
        NOMAD::OutputInfo outputInfo("EvaluatorControl", s, NOMAD::OutputLevel::LEVEL_STATS);
        outputInfo.setStatsInfo(std::move(stats));
        NOMAD::OutputQueue::Add(std::move(outputInfo));
    }
    OUTPUT_STATS_END
}


// Eval a block (vector) of EvalQueuePointPtr
bool NOMAD::EvaluatorControl::evalBlock(NOMAD::BlockForEval& blockForEval)
{
    if (blockForEval.empty())
    {
        return false;
    }

    // All EvalPoints in blockForEval have the same mainThreadNum, and are to be evaluated
    // with the same evaluator, using the same EvalType.
    // TODO: idea: make BlockForEval a class with its own field mainThreadNum
    const int mainThreadNum = blockForEval[0]->getThreadAlgo();
    NOMAD::EvalType evalType = getEvalType(mainThreadNum);
    const NOMAD::Double hMax = getHMax(mainThreadNum);

    // Create a block of EvalPoints (Block) from the given block of EvalQueuePoints (BlockForEval),
    // to give to the Evaluator.
    NOMAD::Block block;
    for (auto it = blockForEval.begin(); it < blockForEval.end(); it++)
    {
        block.push_back(*it);
    }

    // Block evaluation
    std::vector<bool> evalOk = evalBlockOfPoints(block, mainThreadNum, hMax);

    for (size_t i = 0; i < blockForEval.size(); i++)
    {
        NOMAD::EvalPoint evalPoint = *block[i];
        // Put points in evaluatedPoints list for threadAlgo,
        // so the user can get back evaluation info (especially if cache is not used)
        addEvaluatedPoint(mainThreadNum, evalPoint);
        // Update Eval in BlockForEval. Workaround.
        NOMAD::Eval* eval = evalPoint.getEval(evalType);
        if (nullptr != eval)
        {
            blockForEval[i]->setEval(*eval, evalType);
        }

        computeSuccess(blockForEval[i], evalOk[i], hMax);

    }

    size_t nbEvalOk = std::count(evalOk.begin(), evalOk.end(), true);

    return (nbEvalOk > 0);
}


// Helper to update EvalQueuePointPtr
void NOMAD::EvaluatorControl::computeSuccess(NOMAD::EvalQueuePointPtr evalQueuePoint,
                                             const bool evalOk,
                                             const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::UNSUCCESSFUL;
    NOMAD::EvalType evalType = evalQueuePoint->getEvalType();
    auto mainThreadNum = evalQueuePoint->getThreadAlgo();

    if (evalOk)
    {
        NOMAD::EvalPointPtr xFeas, xInf;
        auto barrier = getBarrier(mainThreadNum);
        if (nullptr != barrier)
        {
            // Use best xFeas and xInf for comparison.
            // Only the Eval part is used.
            xFeas = barrier->getRefBestFeas();
            xInf  = barrier->getRefBestInf();
        }

        auto computeSuccessType = getMainThreadInfo(mainThreadNum).getComputeSuccessType();
        if (evalQueuePoint->isFeasible(evalType))
        {
            // Feasible - Compare with xFeas
            success = computeSuccessType(evalQueuePoint, xFeas);
        }
        else
        {
            // Infeasible - Compare with xInf
            success = computeSuccessType(evalQueuePoint, xInf, hMax);
        }
    }

    evalQueuePoint->setSuccess(success);

    OUTPUT_DEBUG_START
    std::string s = NOMAD::evalTypeToString(evalType) + " Evaluation done for ";
    s += evalQueuePoint->displayAll();
    s += ". Success found: " + NOMAD::enumStr(evalQueuePoint->getSuccess());
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END
}


// Eval a single EvalPoint, calling EvcMainThreadInfo's Evaluator.
// Underlying EvalPoint must already be in the cache.
// Updates only the Eval members of the EvalPoint.
bool NOMAD::EvaluatorControl::evalSinglePoint(NOMAD::EvalPoint &evalPoint,
                                              const int mainThreadNum,
                                              const NOMAD::Double &hMax)
{
    bool evalOk = false;

    // Create a block of one point and evaluate it.
    NOMAD::Block block;
    std::shared_ptr<NOMAD::EvalPoint> epp = std::make_shared<NOMAD::EvalPoint>(evalPoint);
    block.push_back(epp);
    std::vector<bool> vectorEvalOk = evalBlockOfPoints(block, mainThreadNum, hMax);
    size_t nbEvalOk = std::count(vectorEvalOk.begin(), vectorEvalOk.end(), true);
    evalOk = (nbEvalOk > 0);

    // Copy back values to evalPoint. Yes, that is ugly.
    // Disclaimer: We expect this method to be called only from Initialization Step.
    evalPoint = *epp;

    return evalOk;
}


// Eval a block of EvalPoints, calling Evaluator.
// Underlying EvalPoints must already be in the cache.
// Updates only the Eval members of the EvalPoints.
// Returns a vector of bools, of the same size as block.
std::vector<bool> NOMAD::EvaluatorControl::evalBlockOfPoints(
                                    NOMAD::Block &block,
                                    const int mainThreadNum,
                                    const NOMAD::Double &hMax)
{
    auto evalType = getEvalType(mainThreadNum);

    std::vector<bool> evalOk(block.size(), false);
    std::vector<bool> countEval(block.size(), false);

    for (auto it = block.begin(); it != block.end(); it++)
    {
        NOMAD::EvalPoint evalPoint = *(*it).get();
        if (!evalPoint.ArrayOfDouble::isComplete())
        {
            std::string err("EvaluatorControl: Eval Single Point: Trying to evaluate an undefined Point ");
            err += evalPoint.getX()->NOMAD::Point::display();
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }
    }

    // NOTE: Put everything in a single OutputInfo structure and then
    // add it to the OutputQueue all at once when the thread part is done.
    std::string s_thread_info;
    OUTPUT_DEBUG_START
    s_thread_info = "Thread = " + NOMAD::itos(NOMAD::getThreadNum());
    OUTPUT_DEBUG_END
    // Use "EvaluatorControl" as originator.
    NOMAD::OutputInfo evalInfo("EvaluatorControl", s_thread_info, NOMAD::OutputLevel::LEVEL_DEBUG);

    // Evaluation of the block
    try
    {
        for (size_t i = 0; i < block.size(); i++)
        {
            std::shared_ptr<NOMAD::EvalPoint> evalPoint = block[i];
            updateEvalStatusBeforeEval(*evalPoint, mainThreadNum);
        }
        OUTPUT_INFO_START
        std::string startMsg = "Start evaluation of block of " + NOMAD::itos(block.size()) + " points.";
        evalInfo.addMsg(startMsg);
        OUTPUT_INFO_END
#ifdef TIME_STATS
        double evalStartTime = NOMAD::Clock::getCPUTime();
#endif // TIME_STATS
        evalOk = getMainThreadInfo(mainThreadNum).getEvaluator()->eval_block(block, hMax, countEval);
#ifdef TIME_STATS
#ifdef _OPENMP
        #pragma omp critical(computeEvalTime)
#endif // _OPENMP
        {
            _evalTime += NOMAD::Clock::getCPUTime() - evalStartTime;
        }
#endif // TIME_STATS

        // Note: _bbEval, and EvcMainThreadInfo members _lapBbEval, _sgteEval, _subBbEval, are atomic.
        // Update these counters with evaluations that count.
        auto nbCountEval = std::count(countEval.begin(), countEval.end(), true);
        if (EvalType::SGTE == evalType)
        {
            incSgteEval(mainThreadNum, nbCountEval);
            _totalSgteEval += nbCountEval;
        }
        else
        {
            _bbEval += nbCountEval;
            getMainThreadInfo(mainThreadNum).incLapBbEval(nbCountEval);

            // All bb evals count for _nbEvalSentToEvaluator.
            _nbEvalSentToEvaluator += block.size();
        }

        for (size_t index = 0; index < block.size(); index++)
        {
            std::shared_ptr<NOMAD::EvalPoint> evalPoint = block[index];
            // Update eval status if needed.
            updateEvalStatusAfterEval(*evalPoint, mainThreadNum, evalOk[index]);
            // increment counter
            evalPoint->incNumberEval();
        }
    }
    catch (std::exception &e)
    {
        std::string err("EvaluatorControl: Eval Block of Points: eval_x returned an exception: ");
        err += e.what();
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (size_t index = 0; index < block.size(); index++)
    {
        NOMAD::EvalPoint evalPoint = *block[index].get();

        if (evalOk[index] && !evalPoint.getF(evalType).isDefined())
        {
            std::string modifMsg = "Warning: EvaluatorControl: Point ";
            modifMsg += evalPoint.display() + ": Eval ok but f not defined. Setting evalOk to false.";
            OUTPUT_INFO_START
            evalInfo.addMsg(modifMsg);
            OUTPUT_INFO_END
            std::cerr << modifMsg << std::endl;

            evalOk[index] = false;
            evalPoint.setEvalStatus(NOMAD::EvalStatusType::EVAL_FAILED, evalType);
        }

        if (evalOk[index] && nullptr == evalPoint.getEval(evalType))
        {
            // Error: evalOk is true, but Eval is NULL.
            std::string err = "EvaluatorControl: Eval Single Point: no Eval on EvalPoint that was just evaluated. " + evalPoint.display();
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        if (getUseCache(mainThreadNum))
        {
            // Update cache
            NOMAD::CacheBase::getInstance()->update(evalPoint, evalType);
        }
    }

    if (NOMAD::EvalType::BB == evalType)
    {
        // One more block evaluated (count only blocks of BB evaluations).
        _blockEval++;
    }

    OUTPUT_INFO_START
    NOMAD::OutputQueue::Add(std::move(evalInfo));
    OUTPUT_INFO_END

    return evalOk;
}


void NOMAD::EvaluatorControl::updateEvalStatusBeforeEval(NOMAD::EvalPoint &evalPoint,
                                                         const int mainThreadNum)
{
    std::string err;
    // Find the EvalPoint in the cache and set its eval status to IN_PROGRESS.
    NOMAD::EvalPoint foundEvalPoint;
    if (getUseCache(mainThreadNum))
    {
        if (!NOMAD::CacheBase::getInstance()->find(evalPoint, foundEvalPoint))
        {
            err = "NOMAD::EvaluatorControl: updateEvalStatusBeforeEval: EvalPoint not found: ";
            err += evalPoint.display();
            //throw NOMAD::Exception(__FILE__, __LINE__, err);
            // Don't throw an exception, for now - make it a Warning.
            err = "Warning: " + err;
            NOMAD::OutputQueue::Add(err);
            return;
        }
    }
    else
    {
        foundEvalPoint = evalPoint;
    }

    NOMAD::EvalType evalType = getEvalType(mainThreadNum);
    NOMAD::EvalStatusType evalStatus = foundEvalPoint.getEvalStatus(evalType);
    if (evalStatus == NOMAD::EvalStatusType::EVAL_FAILED
        || evalStatus == NOMAD::EvalStatusType::EVAL_ERROR
        || evalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
        || evalStatus == NOMAD::EvalStatusType::EVAL_CONS_H_OVER
        || evalStatus == NOMAD::EvalStatusType::EVAL_OK)
    {
        if (NOMAD::EvalType::BB == evalType)
        {
            err = "Warning: Point " + foundEvalPoint.display() + " will be re-evaluated.";
            NOMAD::OutputQueue::Add(err, NOMAD::OutputLevel::LEVEL_WARNING);
        }
    }
    else if (evalStatus == NOMAD::EvalStatusType::EVAL_IN_PROGRESS)
    {
        // Should not happen
        err = "NOMAD::EvaluatorControl: updateEvalStatusBeforeEval: ";
        err += "Evaluation of EvalPoint ";
        err += evalPoint.getX()->NOMAD::Point::display();
        err += " is already in progress";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    else if (evalStatus == NOMAD::EvalStatusType::EVAL_NOT_STARTED
             || evalStatus == NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED)
    {
        // All good
    }
    else
    {
        // Failsafe
        err = "Unknown eval status: " + NOMAD::enumStr(evalStatus);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    evalPoint.setEvalStatus(NOMAD::EvalStatusType::EVAL_IN_PROGRESS, evalType);
    if (getUseCache(mainThreadNum))
    {
        NOMAD::CacheBase::getInstance()->update(evalPoint, evalType);
    }
}


void NOMAD::EvaluatorControl::updateEvalStatusAfterEval(NOMAD::EvalPoint &evalPoint,
                                                        const int mainThreadNum,
                                                        bool evalOk)
{
    NOMAD::EvalType evalType = getEvalType(mainThreadNum);
    NOMAD::EvalStatusType evalStatus = evalPoint.getEvalStatus(evalType);
    if (evalStatus == NOMAD::EvalStatusType::EVAL_FAILED
        || evalStatus == NOMAD::EvalStatusType::EVAL_ERROR
        || evalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
        || evalStatus == NOMAD::EvalStatusType::EVAL_CONS_H_OVER
        || evalStatus == NOMAD::EvalStatusType::EVAL_OK)
    {
        // Nothing to do
    }
    else if (evalStatus == NOMAD::EvalStatusType::EVAL_IN_PROGRESS)
    {
        // Evaluator did not update evalStatus. Lean on evalOk.
        evalPoint.setEvalStatus(evalOk ? NOMAD::EvalStatusType::EVAL_OK : NOMAD::EvalStatusType::EVAL_FAILED, evalType);
    }
    else if (evalStatus == NOMAD::EvalStatusType::EVAL_NOT_STARTED
             || evalStatus == NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED)
    {
        std::string err = "Eval status after evaluation is: " + NOMAD::enumStr(evalStatus);
        err += ". Cannot be handled.";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    else
    {
        std::string err = "Unknown eval status: " + NOMAD::enumStr(evalStatus);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}
