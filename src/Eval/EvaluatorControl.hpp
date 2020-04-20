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
/**
 \file   EvaluatorControl.hpp
 \brief  Management for evaluation.
 \author Viviane Rochon Montplaisir
 \date   November 2017
 \see    EvaluatorControl.cpp
 */

#ifndef __NOMAD400_EVALUATORCONTROL__
#define __NOMAD400_EVALUATORCONTROL__

#include "../Cache/CacheBase.hpp"
#include "../Eval/Barrier.hpp"
#include "../Eval/EvalQueuePoint.hpp"
#include "../Eval/Evaluator.hpp"
#include "../Param/EvaluatorControlParameters.hpp"

#include "../Algos/AllStopReasons.hpp"

#include <atomic>       // For atomic
#include <time.h>
#include <vector>

#include "../nomad_nsbegin.hpp"


/// \brief Class to control the evaluation of points using a queue.
/**
 * \todo Complete the description.
 */
class EvaluatorControl {
private:
    
    std::unique_ptr<Evaluator> _evaluator;///< The Evaluator for evaluating points.

    std::shared_ptr<EvaluatorControlParameters> _evalContParams;  ///< The parameters controlling the behavior of the class
    
    /// The queue of points to be evaluated.
      /**
       \todo Have means to reorder the queue and change the comparison function during the run. \n
     */
    /**
     * The queue is implemented as a vector. Points are added at the end of
      the queue. Points are sorted using the _comp() which is called
      in unlockQueue() when the user is done adding points to the queue. \n
     * Sorting the queue can also be done by providing another function.
     */
    std::vector<EvalQueuePointPtr> _evalPointQueue;
    ComparePriority _comp;

#ifdef _OPENMP
    /// To lock the queue
    omp_lock_t _evalQueueLock;
#endif // _OPENMP

    /**
     * Barrier used for the current runs. \n
     * Modified only by master thread. Used by all threads.
     */
    std::shared_ptr<Barrier> _barrier;
    bool _opportunisticEval; ///< Is opportunistic ?
    bool _updateCache; ///< Update cache before and after optimization
    std::vector<EvalPoint> _evaluatedPoints; ///< Where evaluated points are put temporarily

    SuccessType _success; ///< Success type of the last run

    std::atomic<size_t> _currentlyRunning; ///< Count number of evaluations currently running

    /// The number of blackbox evaluations performed
    /**
     \remark
     * Atomic for thread-safety. \n
     * Atomic is also supported by OpenMP, but here
     * I find it clearer to declare it atomic from the start.
     */
    std::atomic<size_t> _bbEval;

    /// The number of blackbox evaluations performed by a given sub algorithm (reset at Algorithm start).
    std::atomic<size_t> _lapBbEval;

    /// The maximum number of blackbox evaluations that can be performed by a sub algorithm.
    /**
     * For now we have a single counter used for a sub-algorithm.\n
     * The main counter EvaluatorControl::_bbEval is used for the main algo. \n
     * \remark Optionnaly set during Algorithm creation.
     */
    std::atomic<size_t> _lapMaxBbEval;

    /// The number of sgte evaluations performed (since last reset)
    std::atomic<size_t> _sgteEval;

    /// The total number of sgte evaluations. Used for stats exclusively.
    std::atomic<size_t> _totalSgteEval;

    /// The number of block evaluations performed
    /**
     \remark Atomic for thread-safety.
     */
    std::atomic<size_t> _blockEval;

    /**
     The number of evaluations performed through \c this:
        - Including evaluations for which countEval returned \c false.
        - Not including cache hits.
     */
    std::atomic<size_t> _nbEvalSentToEvaluator;
    
    bool _doneWithEval;     ///< All evaluations done. The queue can be destroyed.

#ifdef TIME_STATS
    double _evalTime;  ///< Total time spent running evaluations
#endif // TIME_STATS

public:
    
    /// Constructor
    /**
     \param evaluator       The blackbox evaluator -- \b IN.
     \param evalContParams  The parameters controlling how the class works -- \b IN.
     \param comp            The priority comparison function -- \b IN.
     */
    explicit EvaluatorControl(std::unique_ptr<Evaluator> evaluator,
                              const std::shared_ptr< EvaluatorControlParameters> &evalContParams,
                              ComparePriority comp = ComparePriority() )
      : _evaluator(std::move(evaluator)),
        _evalContParams(evalContParams),
        _evalPointQueue(),
        _comp(comp),
#ifdef _OPENMP
        _evalQueueLock(),
#endif // _OPENMP
        _barrier(),
        _opportunisticEval(false),
        _updateCache(true),
        _evaluatedPoints(),
        _success(SuccessType::UNSUCCESSFUL),
        _currentlyRunning(0),
        _bbEval(0),
        _lapBbEval(0),
        _lapMaxBbEval(INF_SIZE_T),
        _sgteEval(0),
        _totalSgteEval(0),
        _blockEval(0),
        _nbEvalSentToEvaluator(0),
        _doneWithEval(false)
#ifdef TIME_STATS
        ,_evalTime(0.0)
#endif // TIME_STATS
    {
        init();
    }

    /// Destructor.
    virtual ~EvaluatorControl()
    {
        destroy();
    }

    /*---------------*/
    /* Get/Set       */
    /*---------------*/
    /// Outside of the class, the only way to get the Evaluator is to acquire it.
    std::unique_ptr<Evaluator> getEvaluatorUPtr() { return std::move(_evaluator); }
    void setEvaluator(std::unique_ptr<Evaluator> evaluator)
    {
        _evaluator = std::move(evaluator);
    }
    
    /// Get the number of blackbox evaluations.
    size_t getBbEval() const { return _bbEval; }
    
    /// Set blackbox evaluations number.
    void setBbEval(const size_t bbEval) { _bbEval = bbEval; }

    size_t getSgteEval() const { return _sgteEval; }
    size_t getTotalSgteEval() const { return _totalSgteEval; }
    void resetSgteEval() { _sgteEval = 0; }

    /// Get the number of block evaluations.
    size_t getBlockEval() const { return _blockEval; }
    
    /** Get the total number of evaluations.
     * Total number of evaluations, including:
        - blackbox evaluations (EvaluatorControl::_bbEval),
        - blackbox evaluations for which countEval returned \c false
        - cache hits.
       Not including sgte evaluations (EvaluatorControl::_sgteEval).
     
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
    
    void setLapMaxBbEval(const size_t maxBbEval) { _lapMaxBbEval = maxBbEval; }
    void resetLapBbEval ( void ) { _lapBbEval = 0; }
    size_t getLapBbEval ( void ) { return _lapBbEval; }

#ifdef TIME_STATS
    /// Get the total time spend in Evaluator
    double getEvalTime() const { return _evalTime; }
#endif // TIME_STATS
    
    size_t getQueueSize() const { return _evalPointQueue.size(); }

    void setBarrier(const std::shared_ptr<Barrier> barrier) { _barrier = barrier; }
    const std::shared_ptr<Barrier> getBarrier() const { return _barrier; }
    
    bool getOpportunisticEval() const;
    bool getUseCache() const;

    /// Get all points that were just evaluated. This is especially useful when cache is not used.
    /**
     \note              _evaluatedPoints is cleared
     */
    std::vector<EvalPoint> getAllEvaluatedPoints();


    /// Get the max infeasibility to keep a point in barrier
    Double getHMax() const
    {
        return (nullptr == _barrier) ? INF : _barrier->getHMax();
    }

    /// Get the parameters for \c *this
    std::shared_ptr<EvaluatorControlParameters> getEvaluatorControlParams() const
    {
        return _evalContParams;
    }

    /// Get the Evaluator's parameters.
    std::shared_ptr<EvalParameters> getEvalParams() const;

    /*---------------*/
    /* Other methods */
    /*---------------*/
    
    /// Lock the queue.
    /**
     Typically, this is done before adding points.
     It is also used to modify parameters.
     */
    void lockQueue();
  
    /// Unlock the queue.
    /**
     We are done adding points, or modify and checking the parameters.
     If doSort is \c true, the points in the queue are sorted using _comp.
     */
    void unlockQueue(const bool doSort = true);
  
    /// Add a single point to the queue
    void addToQueue(const EvalQueuePointPtr& evalQueuePoint);
    
    /// Get the top point from the queue and pop it.
    /** 
    \param evalQueuePoint   The eval point popped from the queue -- \b OUT.
    \return                 \c true if it worked, \c false otherwise.
    */
    bool popEvalPoint(EvalQueuePointPtr &evalQueuePoint);

    /// Pop eval points from the queue to fill a block of size BB_MAX_BLOCK_SIZE.
    /** 
    \param block   The eval queue point block created -- \b OUTr.
    \return        \c true if the block has at least one point, \c false otherwise.
    */
    bool popBlock(BlockForEval &block);

    /// Sort the queue with respect to the comparison function comp.
    void sort(ComparePriority comp);
  
    /// Use the default comparison function _comp.
    void sort() { sort(_comp); }

    /// Clear queue.
    /**
     * If waitRunning is \c true, wait for all currently running blackboxes to end
     before clearing. \n
     * If showDebug is \c true, print points remaining in the queue.
    */
    void clearQueue(const bool waitRunning = false, const bool showDebug = false);

    /// Start evaluation.
    void start() {}
  
    /// Continuous evaluation - running on all threads simultaneously.
    /**
     * Stop reasons may be controled by parameters MAX_BB_EVAL, MAX_EVAL, OPPORTUNISTIC_EVAL. \n
     * If strategy is opportunistic, stop as soon as a successful point is found. \n
     \return    The success type of the evaluations.
     */
    SuccessType run();
  
    /// Stop evaluation
    void stop() ;

    /// Restart
    void restart() { _doneWithEval = false; }

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
     \param hMax        The max infeasibility for keeping points in barrier -- \b IN.
     \return            \c true if evaluation worked (evalOk), \c false otherwise.
     */
    bool evalSinglePoint(EvalPoint &evalPoint, const Double &hMax = INF);

    /// Evaluates a block of points.
    /**
     For each point, \c true if evaluation worked \c false otherwise.
     \remark Updates only the Eval members.
     
     \param block   The block of points to evaluate -- \b IN/OUT.
     \param hMax    The max infeasibilyt to keep a point in barrier -- \b IN.
     \return        A vector of booleans, of the same size as block.
     */
    std::vector<bool> evalBlockOfPoints(Block &block,
                                        const Double &hMax = INF);
    
    /// Updates eval status.
    /**
     Find the point in the cache and update its evalStatus
      to IN_PROGRESS, knowing that the evaluation is about to start.
     
     \param evalPoint   The evaluation point -- \b IN/OUT.
     */
    void updateEvalStatusBeforeEval(EvalPoint &evalPoint);

    /// Updates eval status.
    /**
     Update point's evalStatus, knowing that the evaluation has just ended.
     
     \param evalPoint    The evalPoint -- \b IN/OUT.
     \param evalOk       Status of evaluation -- \b IN.
     */
    void updateEvalStatusAfterEval(EvalPoint &evalPoint, bool evalOk);

    /// Did we reach one of the evaluation parameters: MAX_EVAL, MAX_BB_EVAL, MAX_BLOCK_EVAL ?
    bool reachedMaxEval() const;

    /// Did we reach Max LAP, MAX_SGTE_EVAL (temporary max evals)?
    bool reachedMaxStepEval() const;


private:

    /// Helper for constructor
    void init();
    
    /// Helper for destructor
    void destroy();
 
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

    /// Did we reach a stop condition (for main thread)?
    bool stopMainEval();

    /// Debug trace when a thread is waiting.
    void displayDebugWaitingInfo(time_t &lastDisplayed) const;

    /// Stats Output
    void AddStatsInfo(const BlockForEval& block) const;

    const EvalType& getEvalType() const { return _evaluator->getEvalType(); }
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVALUATORCONTROL__
