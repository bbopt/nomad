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
 \file   EvaluatorControl.hpp
 \brief  Management for evaluation.
 \author Viviane Rochon Montplaisir
 \date   November 2017
 \see    EvaluatorControl.cpp
 */

#ifndef __NOMAD_4_0_EVALUATORCONTROL__
#define __NOMAD_4_0_EVALUATORCONTROL__

#include "../Eval/Barrier.hpp"
#include "../Eval/ComparePriority.hpp"
#include "../Eval/EvalQueuePoint.hpp"
#include "../Eval/EvcMainThreadInfo.hpp"
#include "../Param/EvaluatorControlGlobalParameters.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP
#include <time.h>

#include "../nomad_nsbegin.hpp"


/// \brief Class to control the evaluation of points using a queue.
/**
 * \todo Complete the description.
 */
class EvaluatorControl {
private:
    std::shared_ptr<EvaluatorControlGlobalParameters> _evalContGlobalParams;  ///< The parameters controlling the behavior of the class

    std::set<int> _mainThreads;     // Thread numbers of main threads

    mutable std::map<int, EvcMainThreadInfo> _mainThreadInfo;  ///< Info about main threads

    /// The queue of points to be evaluated.
      /**
       \todo Have means to reorder the queue and change the comparison function during the run. \n
     */
    /**
     * The queue is implemented as a vector. Points are added at the end of
      the queue. Points are sorted using ComparePriority which is called
      in unlockQueue() when the user is done adding points to the queue. \n
     * Sorting the queue can also be done by providing another function.
     */
    std::vector<EvalQueuePointPtr> _evalPointQueue;

    std::shared_ptr<ComparePriorityMethod>  _userCompMethod;    ///< User-implemented comparison method to sort points before evaluation

#ifdef _OPENMP
    /// To lock the queue
    mutable omp_lock_t _evalQueueLock;
#endif // _OPENMP

    /// The number of blackbox evaluations performed
    /**
     \remark
     * Atomic for thread-safety. \n
     * Atomic is also supported by OpenMP, but here
     * I find it clearer to declare it atomic from the start.
     */
    std::atomic<size_t> _bbEval;

    /// The number of blackbox evaluations that are not ok
    /**
     \remark
     * Atomic for thread-safety. \n
     * Atomic is also supported by OpenMP, but here
     * I find it clearer to declare it atomic from the start.
     */
    std::atomic<size_t> _bbEvalNotOk;

    /// The number of feasible blackbox evaluations obtained
    /**
     \remark
     * Atomic for thread-safety. \n
     * Atomic is also supported by OpenMP, but here
     * I find it clearer to declare it atomic from the start.
     */
    std::atomic<size_t> _feasBBEval;

    /// The number of infeasible blackbox evaluations obtained
    /**
     \remark
     * Atomic for thread-safety. \n
     * Atomic is also supported by OpenMP, but here
     * I find it clearer to declare it atomic from the start.
     */
    std::atomic<size_t> _infBBEval;

    /// The number of static surrogate evaluations performed
    std::atomic<size_t> _surrogateEval;

    /// The total number of quad or modellib model evaluations. Used for stats exclusively.
    std::atomic<size_t> _totalModelEval;

    /// The number of block evaluations performed
    /**
     \remark Atomic for thread-safety.
     */
    std::atomic<size_t> _blockEval;

    /// The index of the last successfull evaluation block
    /**
     \remark Atomic for thread-safety.
     */
    std::atomic<size_t> _indexSuccBlockEval;

    /// The index of the best feasible evaluation
    /**
     \remark Atomic for thread-safety.
     */
    std::atomic<size_t> _indexBestFeasEval;

    /// The index of the best infeasible evaluation
    /**
     \remark Atomic for thread-safety.
     */
    std::atomic<size_t> _indexBestInfeasEval;

    /**
     The number of evaluations performed through \c this:
        - Including evaluations for which countEval returned \c false.
        - Not including cache hits.
     */
    std::atomic<size_t> _nbEvalSentToEvaluator;

    /**
     The number of full success evaluations performed through \c this.
     */
    std::atomic<size_t> _nbRelativeSuccess;

    /**
     The number of PhaseOne success evaluations performed through \c this.
     */
    std::atomic<size_t> _nbPhaseOneSuccess;

    /**
     The VNS Mads search neighborhood parameter
     */
    std::atomic<size_t> _VNSMadsNeighParameter;

    bool _allDoneWithEval;     ///< All evaluations done. The queue can be destroyed.

#ifdef TIME_STATS
    double _evalTime;  ///< Total time spent running evaluations
#endif // TIME_STATS

public:
    /// Constructor
    /**
     \param evaluator       The blackbox evaluator -- \b IN.
     \param evalContGlobalParams  The parameters controlling how the class works -- \b IN.
     \param evalContParams  The parameters for main threads -- \b IN.
     */
    explicit EvaluatorControl(std::shared_ptr<Evaluator> evaluator,
                              const std::shared_ptr<EvaluatorControlGlobalParameters>& evalContGlobalParams,
                              const std::shared_ptr<EvaluatorControlParameters>& evalContParams)
      : _evalContGlobalParams(evalContGlobalParams),
        _mainThreads(),
        _mainThreadInfo(),
        _evalPointQueue(),
        _userCompMethod(nullptr),
#ifdef _OPENMP
        _evalQueueLock(),
#endif // _OPENMP
        _bbEval(0),
        _bbEvalNotOk(0),
        _feasBBEval(0),
        _infBBEval(0),
        _surrogateEval(0),
        _totalModelEval(0),
        _blockEval(0),
        _indexSuccBlockEval(0),
        _indexBestFeasEval(0),
        _indexBestInfeasEval(0),
        _nbEvalSentToEvaluator(0),
        _nbRelativeSuccess(0),
        _nbPhaseOneSuccess(0),
        _VNSMadsNeighParameter(0),
        _allDoneWithEval(false)
#ifdef TIME_STATS
        ,_evalTime(0.0)
#endif // TIME_STATS
    {
        init(evaluator, evalContParams);
    }

    /// Destructor.
    virtual ~EvaluatorControl()
    {
        destroy();
    }

    /*---------------*/
    /* Get/Set       */
    /*---------------*/
    void addMainThread(const int threadNum,
                       const std::shared_ptr<StopReason<EvalMainThreadStopType>> evalMainThreadStopReason,
                       const std::shared_ptr<Evaluator>& evaluator,
                       const std::shared_ptr<EvaluatorControlParameters>& evalContParams);
    bool isMainThread(const int threadNum) const { return (_mainThreads.end() != _mainThreads.find(threadNum)); }

    const std::set<int>& getMainThreads() const { return _mainThreads; }
    int getNbMainThreads() const { return int(_mainThreads.size()); }

    std::shared_ptr<Evaluator> setEvaluator(std::shared_ptr<Evaluator> evaluator);

    /// Get the number of blackbox evaluations.
    size_t getBbEval() const { return _bbEval; }

    /// Set blackbox evaluations number.
    void setBbEval(const size_t bbEval) { _bbEval = bbEval; }

    /// Get the number of blackbox evaluations that are not ok.
    size_t getBbEvalNotOk() const { return _bbEvalNotOk; }

    /// Get the number of feasible blackbox evaluations.
    size_t getFeasBbEval() const { return _feasBBEval; }

    /// Get the number of infeasible blackbox evaluations.
    size_t getInfeasBbEval() const { return _infBBEval; }

    /// Get the number of surrogate evaluations.
    size_t getSurrogateEval() const { return _surrogateEval; }

    size_t getModelEval(const int mainThreadNum = -1) const;
    void resetModelEval(const int mainThreadNum = -1);
    size_t getTotalModelEval() const { return _totalModelEval; }

    size_t getBbEvalInSubproblem(const int mainThreadNum = -1) const;
    void resetBbEvalInSubproblem(const int mainThreadNum = -1);

    /// Get the number of block evaluations.
    size_t getBlockEval() const { return _blockEval; }

    /// Get the index  of block evaluations.
    size_t getIndexSuccBlockEval() const { return _indexSuccBlockEval; }

    /** Get the total number of evaluations.
     * Total number of evaluations, including:
        - blackbox evaluations (EvaluatorControl::_bbEval),
        - blackbox evaluations for which countEval returned \c false
        - cache hits.
       Not including quad and sgtelib model evaluations (EvaluatorControl::_modelEval).

     \note Member EvaluatorControl:: holds only the first two values.
     Cache hits are added in this method.
     */
    size_t getNbEval() const;

    /// Set the number of evaluations
    /**
     Similarly to EvaluatorControl:getNbEval(), the number of cache hits is removed before
    setting EvaluatorControl::_nbEvalSentToEvaluator.
     */
    void setNbEval(const size_t nbEval);

    size_t getLapMaxBbEval(const int mainThreadNum) const;
    void setLapMaxBbEval(const size_t maxBbEval);
    void resetLapBbEval();
    size_t getLapBbEval(const int threadNum = -1) const;

    size_t getVNSMadsNeighParameter() const { return _VNSMadsNeighParameter; }
    void resetVNSMadsNeighParameter() { _VNSMadsNeighParameter = 0; }
    void incrementVNSMadsNeighParameter() { _VNSMadsNeighParameter++ ; }

    size_t getNbPhaseOneSuccess() const {return  _nbPhaseOneSuccess; }
    size_t getNbRelativeSuccess() const {return  _nbRelativeSuccess; }
    size_t getIndexFeasEval() const { return _indexBestFeasEval; }
    size_t getIndexInfeasEval() const { return _indexBestInfeasEval; }

#ifdef TIME_STATS
    /// Get the total time spend in Evaluator
    double getEvalTime() const { return _evalTime; }
#endif // TIME_STATS

    size_t getQueueSize(const int mainThreadNum = -1) const;

    bool getDoneWithEval(const int mainThreadNum) const;
    void setDoneWithEval(const int mainThreadNum, const bool doneWithEval);

    void setBarrier(const std::shared_ptr<Barrier>& barrier);
    const std::shared_ptr<Barrier>& getBarrier(const int threadNum = -1) const;

    void setBestIncumbent(const int mainThreadNum, const std::shared_ptr<EvalPoint>& bestIncumbent);
    const std::shared_ptr<EvalPoint>& getBestIncumbent(const int mainThreadNum) const;

    void setUserCompMethod(const std::shared_ptr<ComparePriorityMethod>& compMethod) { _userCompMethod = compMethod; }

    void setComputeType(const ComputeType& computeType);
    const ComputeType& getComputeType(const int mainThreadNum = -1) const;

    void setLastSuccessfulFeasDir(const std::shared_ptr<Direction>& feasDir);
    void setLastSuccessfulInfDir(const std::shared_ptr<Direction>& infDir);
    const std::shared_ptr<Direction>& getLastSuccessfulFeasDir() const;
    const std::shared_ptr<Direction>& getLastSuccessfulInfDir() const;

    void setStopReason(const int mainThreadNum, const EvalMainThreadStopType& s);
    const StopReason<EvalMainThreadStopType>& getStopReason(const int mainThreadNum) const;
    std::string getStopReasonAsString(const int mainThreadNum) const;
    bool testIf(const EvalMainThreadStopType& s) const;
    bool checkEvalTerminate(const int mainThreadNum) const;

    /// Get all points that were just evaluated. This is especially useful when cache is not used.
    /**
     \note List of evaluated points is cleared
     */
    std::vector<EvalPoint> retrieveAllEvaluatedPoints(const int threadNum = -1);
    void addEvaluatedPoint(const int threadNum, const EvalPoint& evaluatedPoint);
    bool remainsEvaluatedPoints(const int threadNum) const;
    void clearEvaluatedPoints(const int threadNum);

    const SuccessType& getSuccessType(const int threadNum) const;
    void setSuccessType(const int threadNum, const SuccessType& success);

    /// Get the max infeasibility to keep a point in barrier
    Double getHMax(const int threadNum) const;

    /// Get the global parameters for \c *this
    const std::shared_ptr<EvaluatorControlGlobalParameters>& getEvaluatorControlGlobalParams() const { return _evalContGlobalParams; }

    /// Get the Evaluator's parameters.
    std::shared_ptr<EvalParameters> getEvalParams(const int threadNum = -1) const;

    /// Get or Set the value of some parameters
    bool getOpportunisticEval(const int mainThreadNum = -1) const;
    void setOpportunisticEval(const bool opportunistic);
    bool getUseCache(const int mainThreadNum = -1) const;
    void setUseCache(const bool usecache);
    EvalType getEvalType(const int mainThreadNum = -1) const;
    size_t getMaxBbEvalInSubproblem(const int mainThreadNum = -1) const;
    void setMaxBbEvalInSubproblem(const size_t maxBbEval);
    bool getSurrogateOptimization(const int mainThreadNum = -1) const;
    void setSurrogateOptimization(const bool surrogateOptimization);

    /*---------------*/
    /* Other methods */
    /*---------------*/

    /// Lock the queue.
    /**
     Typically, this is done before adding points.
     */
    void lockQueue();

    /// Unlock the queue.
    /**
     We are done adding points.
     If doSort is \c true, the points in the queue are sorted using ComparePriority.
     */
    void unlockQueue(const bool doSort = true, const size_t keepN = INF_SIZE_T, const StepType& removeStepType = StepType::UNDEFINED);

    /// Add a single point to the queue
    /**
     \return \c true if point was inserted in queue
    **/
    bool addToQueue(const EvalQueuePointPtr& evalQueuePoint);

    /// Get the top point from the queue and pop it.
    /**
    \param evalQueuePoint   The eval point popped from the queue -- \b OUT.
    \param evaluator        Evaluator with which this point must be evaluated -- \b IN / OUT
    \param hMax             hMax for evaluation of this point -- \b IN / OUT
    \return                 \c true if it worked, \c false otherwise.
    */
    bool popEvalPoint(EvalQueuePointPtr &evalQueuePoint, Evaluator*& evaluator, Double& hMax);

    /// Pop eval points from the queue to fill a block of size BB_MAX_BLOCK_SIZE.
    /**
    \param block   The eval queue point block created -- \b OUTr.
    \return        \c true if the block has at least one point, \c false otherwise.
    */
    bool popBlock(BlockForEval &block);

    /// Clear queue.
    /**
     * \param mainThreadNum   Clear points generated by this main thread only, If -1, clear all queue.
     * \param showDebug   If \c true, print points erased from the queue.
     * \return Number of points erased
    */
    size_t clearQueue(const int mainThreadNum, const bool showDebug = false);

    /// Start evaluation.
    void start() {}

    /// Continuous evaluation - running on all threads simultaneously.
    /**
     * Stop reasons may be controled by parameters MAX_BB_EVAL, MAX_EVAL, EVAL_OPPORTUNISTIC. \n
     * If strategy is opportunistic, stop as soon as a successful point is found. \n
     \return    The success type of the evaluations.
     */
    SuccessType run();

    /// Stop evaluation
    void stop();

    /// Restart
    void restart();

    /// Evaluates a block of points
    /**
     Updates the Eval members of the evaluation points.
     Also updates the fields specific to EvalQueuePoints.

    \param block   The block of points to evaluate -- \b IN/OUT.
    \return        \c true if at least one evaluation worked (evalOk), \c false otherwise.
     */
    bool evalBlock(BlockForEval& block);

    /// Evaluates a single point.
    /**
     * Creates a block with a single point and evaluates it using EvaluatorControl::evalBlockOfPoints(). \n
     * Updates only the Eval members.
     \param evalPoint   The point to evaluate -- \b IN.
     \param mainThreadNum   Thread number of the main thread that generated this point
     \param hMax        The max infeasibility for keeping points in barrier -- \b IN.
     \return            \c true if evaluation worked (evalOk), \c false otherwise.
     */
    bool evalSinglePoint(EvalPoint &evalPoint,
                         const int mainThreadNum,
                         const Double &hMax = INF);

    /// Evaluates a block of points.
    /**
     For each point, \c true if evaluation worked \c false otherwise.
     \remark Updates only the Eval members.
     \param block   The block of points to evaluate -- \b IN/OUT.
     \param evaluator   Evaluator to be used for all these points
     \param hMax    The max infeasibilyt to keep a point in barrier -- \b IN.
     \return        A vector of booleans, of the same size as block.
     */
    std::vector<bool> evalBlockOfPoints(Block &block,
                                        const Evaluator& evaluator,
                                        const Double &hMax = INF);

    /// Updates eval status.
    /**
     Find the point in the cache and update its evalStatus
      to IN_PROGRESS, knowing that the evaluation is about to start.
     \param evalPoint   The evaluation point -- \b IN/OUT.
     \return \c true if the point must be evaluated, \c false otherwise.
     */
    bool updateEvalStatusBeforeEval(EvalPoint &evalPoint);

    /// Updates eval status.
    /**
     Update point's evalStatus, knowing that the evaluation has just ended.
     \param evalPoint       The evalPoint -- \b IN/OUT.
     \param evalOk          Status of evaluation -- \b IN.
     */
    void updateEvalStatusAfterEval(EvalPoint &evalPoint,
                                   bool evalOk);

    /// Did we reach one of the evaluation parameters: MAX_EVAL, MAX_BB_EVAL, MAX_BLOCK_EVAL ?
    bool reachedMaxEval() const;

    /// Did we reach Max LAP, MODEL_MAX_EVAL (temporary max evals)?
    /**
     * \param mainThreadNum Number of the main thread to update if a stop condition is reached. -1 means use current thread.
     */
    bool reachedMaxStepEval(const int mainThreadNum = -1) const;

    /// For debugging purposes. Show the contents of the evaluation queue.
    void debugDisplayQueue() const;

private:

    /// Helper for constructor
    void init(std::shared_ptr<Evaluator> evaluator,
              const std::shared_ptr<EvaluatorControlParameters>& evalContParams);
    /// Helper for destructor
    void destroy();

    /// Get the EvcMainThreadInfo associated with this thread number.
    /**
     * Get EvcMainThreadInfo associated with this thread number.
     \param mainThreadNum Main thread number. If -1, use current thread number, assuming it is a main thread. -- \b IN.
     \return EvcMainThreadInfo associated with this thread number.
     */
    EvcMainThreadInfo& getMainThreadInfo(const int mainThreadNum = -1) const;

    /// Helper function called during evaluation.
    /**
     Update EvalQueuePointPtr after evaluation.

     \param evalQueuePoint The queue point of interest -- \b IN/OUT.
     \param evalOk         Flag to specific if evaluation was OK -- \b IN.
     \param hMax           The max infeasibility to keep a point in barrier -- \b IN.
     */
    void computeSuccess(EvalQueuePointPtr evalQueuePoint,
                        const bool evalOk,
                        const Double& hMax);

    /// Helper for EvalType: BB OR (SURROGATE and EVAL_AS_SURROGATE)
    bool evalTypeAsBB(const EvalType& evalType, const int mainThreadNum) const;
    /// Helper for EvalType: BB OR SURROGATE. MODEL does not count for some counters.
    bool evalTypeCounts(const EvalType& evalType) const;

    /// Sort the queue with respect to the selected sort strategy. Called by unlockQueue().
    void sort();

    /// Helper for unlockQueue(), to validate if a point may be removed from the queue after sort.
    bool canErase(const EvalQueuePointPtr &evalQueuePoint,
                  const int threadNum,
                  const StepType& removeStepType) const;

    /// Did we reach a stop condition (for a main thread)?
    /**
     \param mainThreadNum    The main thread number to verify --\b IN.
     */
    bool stopMainEval(const int mainThreadNum) const;

    /// Debug trace when a thread is waiting.
    void displayDebugWaitingInfo(time_t &lastDisplayed) const;

    /// Stats Output
    void addStatsInfo(const BlockForEval& block) const;

    /// History and Solution file output
    void addDirectToFileInfo(EvalQueuePointPtr evalQueuePoint) const;
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_EVALUATORCONTROL__
