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
#include "../Eval/EvaluatorControl.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Util/Clock.hpp"
#include <unistd.h> // For usleep

// Initialization of static stop reasons used in the libNomadEval library
NOMAD::StopReason<NOMAD::BaseStopType> NOMAD::AllStopReasons::_baseStopReason = NOMAD::StopReason<NOMAD::BaseStopType>();
NOMAD::StopReason<NOMAD::EvalStopType> NOMAD::AllStopReasons::_evalStopReason = NOMAD::StopReason<NOMAD::EvalStopType>();

/*------------------------*/
/* Class EvaluatorControl */
/*------------------------*/
// Initialize EvaluatorControl class.
// To be called by the Constructor.
void NOMAD::EvaluatorControl::init()
{
#ifdef _OPENMP
    omp_init_lock(&_evalQueueLock);
#endif // _OPENMP

    // Set opportunism.
    // The parameter will be re-read in run(), because it may change.
    _opportunisticEval = _evalContParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL");
}


// Terminate EvaluatorControl class.
// To be called by the Destructor.
void NOMAD::EvaluatorControl::destroy()
{
    if (!_evalPointQueue.empty())
    {
        // Show warnings and debug info.
        // Do not scare the user if display degree is medium or low.
        int displayDegree = NOMAD::OutputQueue::getInstance()->getDisplayDegree();
        if (displayDegree >= 3)
        {
            std::cerr << "Warning: deleting EvaluatorControl with EvalPoints remaining." << std::endl;
        }
        bool showDebug = (displayDegree >= 4);
        clearQueue(false /* waitRunning */, showDebug);
    }

#ifdef _OPENMP
    omp_destroy_lock(&_evalQueueLock);
#endif // _OPENMP
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


std::shared_ptr<NOMAD::EvalParameters> NOMAD::EvaluatorControl::getEvalParams() const
{
    std::shared_ptr<NOMAD::EvalParameters> evalParams = nullptr;
    if (_evaluator)
    {   
        evalParams = _evaluator->getEvalParams();
    }

    return evalParams;
}


void NOMAD::EvaluatorControl::lockQueue()
{
#ifdef _OPENMP
    // Sanity check before locking the queue.
    // 1- Verify we are in master thread.
    if (0 != omp_get_thread_num())
    {
        std::string err = "Error: EvaluatorControl::lockQueue called from non-master thread ";
        err += std::to_string(omp_get_thread_num());
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
    // 1- Verify we are in master thread.
    if (0 != omp_get_thread_num())
    {
        std::string err = "Error: EvaluatorControl::unlockQueue called from non-master thread ";
        err += std::to_string(omp_get_thread_num());
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

#ifndef USE_PRIORITY_QUEUE
    // When using a vector instead of priority, the EvalQueuePoints are added randomly.
    // Sort them, using default sort, if doSort is true (default).
    // In non-opportunistic context, it is useless to sort.
    if (doSort && _opportunisticEval)
    {
        sort(_comp);
    }
#endif // USE_PRIORITY_QUEUE
}


// Add an EvalPoint to the Queue
void NOMAD::EvaluatorControl::addToQueue(const NOMAD::EvalQueuePointPtr &evalQueuePoint)
{
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

#ifdef USE_PRIORITY_QUEUE
    _evalPointQueue.push(evalQueuePoint);
#else

    // CT Lexicographic order instead of inverse lexicographic order
    _evalPointQueue.insert ( _evalPointQueue.begin() , evalQueuePoint );
    //_evalPointQueue.push_back(evalQueuePoint);
#endif // USE_PRIORITY_QUEUE

}


// Get the top EvalPoint from the Queue and pop it
// Return true if it worked, false if it failed.
bool NOMAD::EvaluatorControl::popEvalPoint(NOMAD::EvalQueuePointPtr &evalQueuePoint)
{
    bool success = false;
#ifdef _OPENMP
    omp_set_lock(&_evalQueueLock);
#endif // _OPENMP
    if (!_evalPointQueue.empty())
    {
#ifdef USE_PRIORITY_QUEUE
        evalQueuePoint = std::move(_evalPointQueue.top());
        _evalPointQueue.pop();
#else
        // Remove last element, simulate a "pop".
        evalQueuePoint = std::move(_evalPointQueue[_evalPointQueue.size()-1]);
        _evalPointQueue.erase(_evalPointQueue.end()-1);
#endif // USE_PRIORITY_QUEUE
        success = true;
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
    size_t blockSize = NOMAD::INF_SIZE_T;
    bool gotBlockSize = false;
    

    // NOTE: Some racing conditions may happen with value BB_MAX_BLOCK_SIZE.
    // As a workaround, try to get the attribute until no exception is thrown.
    while (!gotBlockSize)
    {
        try
        {
            blockSize = _evalContParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
            gotBlockSize = true;
        }
        catch (NOMAD::Exception &e)
        {
            // While will loop - Retry
        }
    }

    while (block.size() < blockSize && popWorks)
    {
        NOMAD::EvalQueuePointPtr evalQueuePoint;
        popWorks = popEvalPoint(evalQueuePoint);
        if (popWorks)
        {
            block.push_back(std::move(evalQueuePoint));
            success = true;
        }
    }

    return success;

}


#ifndef USE_PRIORITY_QUEUE
void NOMAD::EvaluatorControl::sort(NOMAD::ComparePriority comp)
{
#ifdef _OPENMP
    omp_set_lock(&_evalQueueLock);
#endif // _OPENMP
    std::sort(_evalPointQueue.begin(), _evalPointQueue.end(), comp);
#ifdef _OPENMP
    omp_unset_lock(&_evalQueueLock);
#endif // _OPENMP
}
#endif // USE_PRIORITY_QUEUE


void NOMAD::EvaluatorControl::clearQueue(const bool waitRunning, const bool showDebug)
{
    // Wait for any currently running evaluations to be done.
    // Otherwise, the queue would be empty but old evaluations
    // could pop up unexpectedtly.
    if (waitRunning)
    {
        while (_currentlyRunning > 0)
        {
            std::string s = "Waiting for " + NOMAD::itos(_currentlyRunning);
            s += " evaluations to complete.";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
            usleep(10000);
        }
    }

#ifdef _OPENMP
    omp_set_lock(&_evalQueueLock);
#endif // _OPENMP

#ifdef USE_PRIORITY_QUEUE
    while (!_evalPointQueue.empty())
    {
        if (showDebug)
        {
            std::string s = "Delete point from queue: ";
            auto evalQueuePoint = std::move(_evalPointQueue.top());
            s += evalQueuePoint.tostring();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        }
        _evalPointQueue.pop();
    }
#else
    // Queue is a vector: Show all (if showDebug is true), then delete all.
    if (showDebug)
    {
        for (auto evalQueuePoint : _evalPointQueue)
        {
            std::string s = "Delete point from queue: ";
            s += evalQueuePoint->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        }
    }
    _evalPointQueue.clear();
#endif // USE_PRIORITY_QUEUE

#ifdef _OPENMP
    omp_unset_lock(&_evalQueueLock);
#endif // _OPENMP
}


// Evaluate all points in the queue, or stop under some conditions.
// If strategy is opportunistic (opportunisticEval), stop
// as soon as a successful point is found, and flush the queue.
//
// Points must already be in the cache.
NOMAD::SuccessType NOMAD::EvaluatorControl::run()
{
    // Master thread only:
    //  - Reset success
    //  - Set Barrier and flag OpportunisticEval
    //  - Print info
#ifdef _OPENMP
    #pragma omp master
#endif // _OPENMP
    {
        // At this point, the threads other than master might have already
        // started evaluating. So the number of points in the queue is already moot.
        //NOMAD::OutputQueue::Add("Eval " + NOMAD::itos(_evalPointQueue.size()) + " points.", NOMAD::OutputLevel::LEVEL_DEBUG);

        _success = NOMAD::SuccessType::UNSUCCESSFUL;

        // An empty eval queue must be accounted for
        if (_evalPointQueue.empty())
        {
            NOMAD::AllStopReasons::set(NOMAD::EvalStopType::EMPTY_LIST_OF_POINTS);
        }

        // Update stop reason.
        if (   NOMAD::AllStopReasons::checkEvalTerminate()
            || NOMAD::AllStopReasons::checkBaseTerminate())
        {
            std::string sStopReason = "EvaluatorControl stop reason: ";
            sStopReason += (NOMAD::AllStopReasons::checkBaseTerminate())
                ? NOMAD::AllStopReasons::getBaseStopReasonAsString() + " (Base) "
                : NOMAD::AllStopReasons::getEvalStopReasonAsString() + " (Eval) ";

            NOMAD::OutputQueue::Add(sStopReason, NOMAD::OutputLevel::LEVEL_DEBUG);
            NOMAD::OutputQueue::Flush();
        }
        else
        {
            // Ready to run.
            // Reset EvalStopType to STARTED.
            NOMAD::AllStopReasons::set(NOMAD::EvalStopType::STARTED);
        }

        _opportunisticEval = _evalContParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL");

        std::string s = "Start evaluation. Opportunism = ";
        s += NOMAD::boolToString(_opportunisticEval);
        s += ", Barrier =";
        s += ((nullptr == _barrier) ? " NULL" : "\n" + _barrier->display(4)); // Display a maximum of 4 xFeas and 4 xInf
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);

    }

    // Queue runs forever on non-master threads.
    // On master, queue runs until stopMainEval() is true.
    bool conditionForStop = false;

    // For debug display
    time_t lastDisplayed = 0;

    // conditionForStop is true if we are in master thread and stopMainEval() returns true.
    // conditionForStop is true in any thread if reachedMaxEval() returns true; otherwise, it is always false.
    while (!conditionForStop && !_doneWithEval)
    {
        // Check for stop conditions
#ifdef _OPENMP
        #pragma omp master
#endif // _OPENMP
        {
            conditionForStop = stopMainEval();
        }
        // If we reached max eval, we also stop (valid for all threads).
        conditionForStop = conditionForStop || reachedMaxEval();

        NOMAD::BlockForEval block;
        if (!conditionForStop && popBlock(block))
        {
            _currentlyRunning += block.size();
            bool evalOk = evalBlock(block);
            _currentlyRunning -= block.size();
            if (evalOk)
            {
                // Update SuccessType
                // _success is a member so it can be shared between threads.
#ifdef _OPENMP
                #pragma omp critical(updateSuccessType)
#endif // _OPENMP
                {
                    for (auto it = block.begin(); it < block.end(); it++)
                    {
                        NOMAD::EvalQueuePointPtr evalQueuePoint = (*it);
                        // Update success type for return
                        if (evalQueuePoint->getSuccess() > _success)
                        {
                            _success = evalQueuePoint->getSuccess();
                            evalQueuePoint->setRelativeSuccess(true);
                        }
                    }
                }   // End critical(updateSuccessType)

                // No need for critical here, stopReason would not be set back to NO_STOP
                // by another thread. It could only be set to another stop reason which
                // would also be valid.
                if (_opportunisticEval && _success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
                {
                    NOMAD::AllStopReasons::set(NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS);
                }

                AddStatsInfo(block);
            }
        }
        else if (!_doneWithEval)
        {
            displayDebugWaitingInfo(lastDisplayed);
        }
        else // Queue is empty and we are doneWithEval
        {
#ifdef _OPENMP
            #pragma omp master
#endif // _OPENMP
            {
                std::string s = "Queue is empty and we are done with evaluations.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            }
            break;
        }

    }   // End of while loop: Exit for master thread.
        // Other threads keep on looping.


#ifdef _OPENMP
#pragma omp master
#endif
    {
        // Special case:
        // stopReason is ALL_POINTS_EVALUATED because there are no points left
        // in the queue.
        // But some points are still being evaluated.
        // We must wait for all points to really be evaluated.
        // Note that when all points are evaluated, _success has the correct
        // value, even if it was modified by those last points being evaluated.
        while (( NOMAD::AllStopReasons::testIf ( NOMAD::EvalStopType::ALL_POINTS_EVALUATED )
                || _evalContParams->getAttributeValue<bool>("CLEAR_EVAL_QUEUE"))
               && _currentlyRunning > 0)
        {
            std::string s = "Waiting for " + NOMAD::itos(_currentlyRunning);
            s += " evaluations to complete.";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
            usleep(10000);

            // Update stopReason in case we found a success
            if (_opportunisticEval && _success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                NOMAD::AllStopReasons::set( NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS );
            }

            // Update stopReason in case we reached the different max eval criterions
            reachedMaxEval();
            reachedMaxStepEval();
        }
    }
    
#ifdef _OPENMP
    #pragma omp master
#endif // _OPENMP
    {
        if (_evalContParams->getAttributeValue<bool>("CLEAR_EVAL_QUEUE"))
        {
            // Remove remaining points from queue, to start fresh next time.
            // Otherwise, we keep on evaluationg the points remaining in the queue.
            NOMAD::OutputQueue::Add("Evaluation is done. Clear queue of " + NOMAD::itos(_evalPointQueue.size()) + " points.");
            clearQueue(true /* waitRunning */, false /* showDebug */);
        }

        

        NOMAD::OutputQueue::Add("EvaluatorControl stop reason: " + NOMAD::AllStopReasons::getEvalStopReasonAsString() , NOMAD::OutputLevel::LEVEL_DEBUG);
        NOMAD::OutputQueue::Flush();
        
    }

    return _success;
}


void NOMAD::EvaluatorControl::stop()
{
    _doneWithEval = true;
}


bool NOMAD::EvaluatorControl::stopMainEval()
{
    // Check for stop conditions for main thread.
    // This function should not be called from other threads. But we
    // do not verify the thread number.

    // Inspect stopReasons
    bool doStopEval = NOMAD::AllStopReasons::checkEvalTerminate() ;

    // Inspect for opportunistic success set in a previous point
    doStopEval = doStopEval || ( NOMAD::AllStopReasons::testIf( NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS) );

    // Update stopReason if queue is empty
    if (_evalPointQueue.empty() && (!doStopEval || NOMAD::AllStopReasons::testIf(NOMAD::EvalStopType::EMPTY_LIST_OF_POINTS)) )
    {
        NOMAD::AllStopReasons::set(NOMAD::EvalStopType::ALL_POINTS_EVALUATED);
        doStopEval = true;
    }

    // Update on max eval, valid for all threads
    doStopEval = doStopEval || reachedMaxStepEval() || reachedMaxEval();
    std::string s = "stopMainEval: return true because: ";
    if (doStopEval)
    {
        s += NOMAD::AllStopReasons::getEvalStopReasonAsString() ;
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }


    // Inspect on base stop reason
    bool doStopBase = NOMAD::AllStopReasons::checkBaseTerminate() ;
    if (doStopBase)
    {
        s += NOMAD::AllStopReasons::getBaseStopReasonAsString() ;
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }

    NOMAD::OutputQueue::Flush();

    return ( doStopBase || doStopEval ) ;
}


// Eval at most a total maxBbEval, maxEval, or maxBlockEval points in the queue.
// If this condition is true, we assume that it will remain true for the
// rest of the optimization.
// If the condition is true temporary, for example lap or sgte evals,
// use stopMainEval().
bool NOMAD::EvaluatorControl::reachedMaxEval() const
{
    bool ret = false;
    // There is a possiblity that these values are modified
    // in master thread, and that getAttributeValue() gets called before
    // checkAndComply(). In that case, catch the exception and ignore
    // the test.
    size_t maxBbEval    = NOMAD::INF_SIZE_T;
    size_t maxEval      = NOMAD::INF_SIZE_T;
    size_t maxBlockEval = NOMAD::INF_SIZE_T;

    try
    {
        maxBbEval       = _evalContParams->getAttributeValue<size_t>("MAX_BB_EVAL");
        maxEval         = _evalContParams->getAttributeValue<size_t>("MAX_EVAL");
        maxBlockEval    = _evalContParams->getAttributeValue<size_t>("MAX_BLOCK_EVAL");
    }
    catch (NOMAD::Exception &e)
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
    #pragma omp master
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
bool NOMAD::EvaluatorControl::reachedMaxStepEval() const
{
    bool ret = false;
    size_t maxSgteEval  = NOMAD::INF_SIZE_T;
    try
    {
        // TODO: Default for NOMAD 3 is 10000. Default for NOMAD 4 is 100.
        // Look if/how we can increment the number of evaluations without 
        // increasing the time too much.
        maxSgteEval = _evalContParams->getAttributeValue<size_t>("SGTELIB_MODEL_EVAL_NB");
    }
    catch (NOMAD::Exception &e)
    {
        // ignore value for maxSgteEval if we could not get it
    }

    std::string s = "Reached sub step stop criterion: ";
    if ((EvalType::SGTE == _evaluator->getEvalType()) && maxSgteEval < NOMAD::INF_SIZE_T && _sgteEval >= maxSgteEval)
    {
        // Reached maxSgteEval, max number of eval in Sgte context.
        NOMAD::AllStopReasons::set(NOMAD::EvalStopType::MAX_SGTE_EVAL_REACHED);
        s += NOMAD::AllStopReasons::getEvalStopReasonAsString() + " " + NOMAD::itos(_sgteEval);
        ret = true;
    }
    else if (_lapMaxBbEval < NOMAD::INF_SIZE_T && _lapBbEval >= _lapMaxBbEval)
    {
        // Reached lapMaxEval.
        NOMAD::AllStopReasons::set(NOMAD::EvalStopType::LAP_MAX_BB_EVAL_REACHED);
        s += NOMAD::AllStopReasons::getEvalStopReasonAsString() + " " + NOMAD::itos(_lapBbEval);
        ret = true;
    }

#ifdef _OPENMP
    #pragma omp master
#endif // _OPENMP
    {
        if (ret)
        {
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
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
        if (NOMAD::OutputQueue::getInstance()->getDisplayDegree() >= 4)
        {
            std::string s = "Thread: " + NOMAD::itos(omp_get_thread_num());
            s += " Waiting for points.";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        }
        lastDisplayed = now;
    }
#endif // _OPENMP
}


void NOMAD::EvaluatorControl::AddStatsInfo(const NOMAD::BlockForEval& block) const
{
    // SGTE optimizations generate a lot of stats output. Do not show them.
    // Only show BB optimizations.
    if (NOMAD::EvalType::BB != getEvalType())
    {
        return;
    }
    for (auto it = block.begin(); it < block.end(); it++)
    {
        NOMAD::EvalQueuePointPtr evalQueuePoint = (*it);

        // Evaluation info for output
        NOMAD::StatsInfoUPtr stats(new NOMAD::StatsInfo());
        size_t n = evalQueuePoint->size();

        // Values that are unknown at this point
        NOMAD::ArrayOfDouble unknownMeshIndex(n);
        int threadNum = 0;
#ifdef _OPENMP
        threadNum = omp_get_thread_num();
#endif // _OPENMP

        // As of February 2019, values for bbEval (BBE) and nbEval (EVAL) make
        // enough sense in threaded environment to be printed out. If at some point
        // they do not, for example if too many iterations end up with the same
        // value, or if the user is dissatisfied with this output, we may update it
        // in a stricter way.

        stats->setObj(evalQueuePoint->getF(NOMAD::EvalType::BB));
        stats->setConsH(evalQueuePoint->getH(NOMAD::EvalType::BB));
        stats->setHMax(getHMax());
        stats->setBBE(_bbEval);
        stats->setLap(_lapBbEval);
        stats->setSgte(_sgteEval);
        stats->setTotalSgte(_totalSgteEval);
        stats->setBlkEva(_blockEval);
        stats->setBlkSize(block.size());
        stats->setBBO(evalQueuePoint->getBBO(NOMAD::EvalType::BB));
        stats->setEval(_nbEvalSentToEvaluator);
        stats->setCacheHits(NOMAD::CacheBase::getNbCacheHits());
        stats->setTime(NOMAD::Clock::getTimeSinceStart());
        stats->setMeshIndex(unknownMeshIndex);
        stats->setMeshSize(evalQueuePoint->getMeshSize());
        stats->setFrameSize(evalQueuePoint->getFrameSize());
        stats->setSol(*(evalQueuePoint->getX()));
        stats->setThreadNum(threadNum);
        stats->setRelativeSuccess(evalQueuePoint->getRelativeSuccess());
        stats->setComment(evalQueuePoint->getComment());
        stats->setGenStep(evalQueuePoint->getGenStep());

        std::string s = "Evaluated point: " + evalQueuePoint->display();
        NOMAD::OutputInfo outputInfo("EvaluatorControl", s, NOMAD::OutputLevel::LEVEL_STATS);
        outputInfo.setStatsInfo(std::move(stats));
        NOMAD::OutputQueue::Add(std::move(outputInfo));
    }
}


// Eval a block (vector) of EvalQueuePointPtr
bool NOMAD::EvaluatorControl::evalBlock(NOMAD::BlockForEval& blockForEval)
{
    const NOMAD::Double hMax = getHMax();
    NOMAD::Block block;
    for (auto it = blockForEval.begin(); it < blockForEval.end(); it++)
    {
        block.push_back(*it);
    }

    std::vector<bool> evalOk = evalBlockOfPoints(block, hMax);

    for (size_t i = 0; i < blockForEval.size(); i++)
    {
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

    if (evalOk)
    {
        NOMAD::EvalPointPtr xFeas, xInf;
        if (nullptr != _barrier)
        {
            // Use first xFeas and xInf - their Eval must be equivalent, and it
            // is the only part that is used for comparison.
            xFeas = _barrier->getFirstXFeas();
            xInf  = _barrier->getFirstXInf();
        }

        NOMAD::ComputeSuccessType computeSuccessType;
        if (evalQueuePoint->isFeasible(getEvalType()))
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

    std::string s = NOMAD::evalTypeToString(getEvalType()) + " Evaluation done for ";
    s += evalQueuePoint->displayAll();
    s += ". Success found: " + NOMAD::enumStr(evalQueuePoint->getSuccess());
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
}


// Eval a single EvalPoint, calling _evaluator.
// Underlying EvalPoint must already be in the cache.
// Updates only the Eval members of the EvalPoint.
bool NOMAD::EvaluatorControl::evalSinglePoint(NOMAD::EvalPoint &evalPoint, const NOMAD::Double &hMax)
{
    bool evalOk = false;

    // Create a block of one point and evaluate it.
    NOMAD::Block block;
    std::shared_ptr<NOMAD::EvalPoint> epp = std::make_shared<NOMAD::EvalPoint>(evalPoint);
    block.push_back(epp);
    std::vector<bool> vectorEvalOk = evalBlockOfPoints(block, hMax);
    size_t nbEvalOk = std::count(vectorEvalOk.begin(), vectorEvalOk.end(), true);
    evalOk = (nbEvalOk > 0);

    // Copy back values to evalPoint. Yes, that is ugly.
    // Disclaimer: We expect this method to be called only from Initialization Step.
    evalPoint = *epp;

    return evalOk;
}


// Eval a block of EvalPoints, calling _evaluator.
// Underlying EvalPoints must already be in the cache.
// Updates only the Eval members of the EvalPoints.
// Returns a vector of bools, of the same size as block.
std::vector<bool> NOMAD::EvaluatorControl::evalBlockOfPoints(
                                    NOMAD::Block &block,
                                    const NOMAD::Double &hMax)
{
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
    std::string s_thread_num = "0";
#ifdef _OPENMP
    s_thread_num = NOMAD::itos(omp_get_thread_num());
#endif // _OPENMP
    std::string s_thread_info = "Thread = " + s_thread_num;
    // Use "EvaluatorControl" as originator.
    NOMAD::OutputInfo evalInfo("EvaluatorControl", s_thread_info, NOMAD::OutputLevel::LEVEL_DEBUG);

    // Evaluation of the block
    try
    {
        for (auto it = block.begin(); it != block.end(); it++)
        {
            std::shared_ptr<NOMAD::EvalPoint> evalPoint = (*it);
            updateEvalStatusBeforeEval(*evalPoint);
        }
        std::string startMsg = "Start evaluation of block of " + NOMAD::itos(block.size()) + " points.";
        evalInfo.addMsg(startMsg);
        evalOk = _evaluator->eval_block(block, hMax, countEval);

        // Note: _bbEval, _lapBbEval, _sgteEval, and _nbEvalSentToEvaluator are atomic.
        // Add evals that count to _bbEval, _lapBbEval and _sgteEval.
        auto nbCountEval = std::count(countEval.begin(), countEval.end(), true);
        if ((EvalType::SGTE == _evaluator->getEvalType()))
        {
            _sgteEval += nbCountEval;
            _totalSgteEval += nbCountEval;
        }
        else
        {
            _bbEval += nbCountEval;
            _lapBbEval += nbCountEval;

            // All bb evals count for _nbEvalSentToEvaluator.
            _nbEvalSentToEvaluator += block.size();
        }

        for (size_t index = 0; index < block.size(); index++)
        {
            std::shared_ptr<NOMAD::EvalPoint> evalPoint = block[index];
            // Update eval status if needed.
            updateEvalStatusAfterEval(*evalPoint, evalOk[index]);
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

        if (evalOk[index] && !evalPoint.getF(getEvalType()).isDefined())
        {
            std::string modifMsg = "Warning: EvaluatorControl: Point ";
            modifMsg += evalPoint.display() + ": Eval ok but f not defined. Setting evalOk to false.";
            evalInfo.addMsg(modifMsg);
            std::cerr << modifMsg << std::endl;

            evalOk[index] = false;
            evalPoint.setEvalStatus(NOMAD::EvalStatusType::EVAL_FAILED, getEvalType());
        }

        if (evalOk[index] && nullptr == evalPoint.getEval(getEvalType()))
        {
            // Error: evalOk is true, but Eval is NULL.
            std::string err = "EvaluatorControl: Eval Single Point: no Eval on EvalPoint that was just evaluated. " + evalPoint.display();
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        // Update cache
        NOMAD::CacheBase::getInstance()->update(evalPoint, getEvalType());
    }

    // One more block evaluated.
    _blockEval++;

    NOMAD::OutputQueue::Add(std::move(evalInfo));

    return evalOk;
}


void NOMAD::EvaluatorControl::updateEvalStatusBeforeEval(NOMAD::EvalPoint &evalPoint)
{
    const NOMAD::EvalType evalType = getEvalType();
    std::string err;
    // Find the EvalPoint in the cache and set its eval status to IN_PROGRESS.
    NOMAD::EvalPoint foundEvalPoint;
    size_t nbFound = NOMAD::CacheBase::getInstance()->find(evalPoint, foundEvalPoint);
    if (nbFound == 0)
    {
        err = "NOMAD::EvaluatorControl: updateEvalStatusBeforeEval: EvalPoint not found: ";
        err += evalPoint.display();
        //throw NOMAD::Exception(__FILE__, __LINE__, err);
        // Don't throw an exception, for now - make it a Warning.
        err = "Warning: " + err;
        NOMAD::OutputQueue::Add(err);
        return;
    }

    NOMAD::EvalStatusType evalStatus = foundEvalPoint.getEvalStatus(evalType);
    if (evalStatus == NOMAD::EvalStatusType::EVAL_FAILED
        || evalStatus == NOMAD::EvalStatusType::EVAL_ERROR
        || evalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
        || evalStatus == NOMAD::EvalStatusType::EVAL_CONS_H_OVER
        || evalStatus == NOMAD::EvalStatusType::EVAL_OK)
    {
        err = "Warning: Point " + foundEvalPoint.display() + " will be re-evaluated.";
        NOMAD::OutputQueue::Add(err, NOMAD::OutputLevel::LEVEL_WARNING);
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
    NOMAD::CacheBase::getInstance()->update(evalPoint, evalType);
}


void NOMAD::EvaluatorControl::updateEvalStatusAfterEval(NOMAD::EvalPoint &evalPoint, bool evalOk)
{
    const NOMAD::EvalType evalType = getEvalType();
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
