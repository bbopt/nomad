/**
 \file   EvcMainThreadInfo.hpp
 \brief  Information about main threads
 \author Viviane Rochon Montplaisir
 \date   June 2020
 \see    EvcMainThreadInfo.cpp
 */

#ifndef __NOMAD400_EVCMAINTHREADINFO__
#define __NOMAD400_EVCMAINTHREADINFO__

#include <atomic>   // for atomic

#include "../Eval/Barrier.hpp"
#include "../Eval/ComputeSuccessType.hpp"
#include "../Eval/Evaluator.hpp"
#include "../Param/EvaluatorControlParameters.hpp"
#include "../Util/AllStopReasons.hpp"

#include "../nomad_nsbegin.hpp"


/// \brief Class for info related to one of the main threads.
class EvcMainThreadInfo
{
private:
    std::shared_ptr<Evaluator> _evaluator;              ///< The Evaluator for either blackbox or surrogate evaluations.
    std::shared_ptr<EvaluatorControlParameters> _evalContParams;  ///< The parameters controlling the behavior of EvaluatorControl for this main thread
    bool                            _doneWithEval;      ///< All evaluations done for this main thread
    std::shared_ptr<Barrier>        _barrier;
    std::vector<EvalPoint>          _evaluatedPoints;   ///< Where evaluated points are put temporarily
    SuccessType                     _success;           ///< Success type of the last run
    std::atomic<size_t>             _currentlyRunning;  ///< Count number of evaluations currently running.
    size_t                          _lapMaxBbEval;      ///< The maximum number of blackbox evaluations that can be performed by a sub algorithm.
    std::atomic<size_t>             _lapBbEval;         ///< The number of blackbox evaluations performed by a given sub algorithm (reset at Algorithm start).
    std::atomic<size_t>             _sgteEval;          ///< The number of sgte evaluations performed (since last reset)
    std::atomic<size_t>             _subBbEval;         ///< The number of bb eval for a subproblem (e.g. SSD-Mads context)
    ComputeSuccessType              _computeSuccessType;    ///< Function used to compute SuccessType. Depends on Evaluator's EvalType.
    StopReason<EvalStopType>        _stopReason;

public:
    /// Constructor
    /**
     \param evaluator       The Evaluator for either blackbox or surrogate evaluations-- \b IN.
     \param evalContParams  The parameters controlling how the EvaluatorControl behaves for this main thread-- \b IN.
     */
    explicit EvcMainThreadInfo(std::shared_ptr<NOMAD::Evaluator> evaluator,
                               const std::shared_ptr< EvaluatorControlParameters> &evalContParams)
      : _evaluator(evaluator),
        _evalContParams(evalContParams),
        _doneWithEval(false),
        _barrier(),
        _evaluatedPoints(),
        _success(SuccessType::UNSUCCESSFUL),
        _currentlyRunning(0),
        _lapMaxBbEval(INF_SIZE_T),
        _lapBbEval(0),
        _sgteEval(0),
        _subBbEval(0),
        _computeSuccessType(ComputeSuccessType(ComputeSuccessType::defaultComputeSuccessType)),
        _stopReason()
    {
        init();
    }

    /// Set Evaluator and return old Evaluator.
    /**
     \param evaluator       The Evaluator for either blackbox or surrogate evaluations-- \b IN.
     \return                The previous Evaluator.
     */
    std::shared_ptr<NOMAD::Evaluator> setEvaluator(std::shared_ptr<NOMAD::Evaluator> evaluator);
    const Evaluator* getEvaluator() { return _evaluator.get(); }
    std::shared_ptr<EvalParameters> getEvalParams() const;
    EvalType getEvalType() const;

    bool getDoneWithEval() const { return _doneWithEval; }
    void setDoneWithEval(const bool doneWithEval) { _doneWithEval = doneWithEval; }

    // Get and set parameters
    bool getOpportunisticEval() const;
    void setOpportunisticEval(const bool opportunisticEval);
    bool getUseCache() const;
    void setUseCache(const bool useCache);
    size_t getMaxBbEvalInSubproblem() const;

    // Get and set counters
    void setLapMaxBbEval(const size_t maxBbEval) { _lapMaxBbEval = maxBbEval; }
    size_t getLapBbEval() const { return _lapBbEval; }
    size_t getLapMaxBbEval() const { return _lapMaxBbEval; }
    void incLapBbEval(const size_t countEval) { _lapBbEval += countEval; }
    void resetLapBbEval() { _lapBbEval = 0; }
    size_t getSgteEval() const { return _sgteEval; }
    void resetSgteEval() { _sgteEval = 0; }
    void incSgteEval(const size_t countEval) { _sgteEval += countEval; }
    size_t getBbEvalInSubproblem() const { return _subBbEval; }
    void incBbEvalInSubproblem(const size_t countEval) { _subBbEval += countEval; }
    void resetBbEvalInSubproblem() { _subBbEval = 0; }

    const std::shared_ptr<Barrier>& getBarrier() const { return _barrier; }
    void setBarrier(const std::shared_ptr<Barrier>& barrier) { _barrier = barrier; }

    std::vector<EvalPoint> retrieveAllEvaluatedPoints();
    void addEvaluatedPoint(const EvalPoint& evaluatedPoint);
    bool remainsEvaluatedPoints() const { return !_evaluatedPoints.empty(); }
    void clearEvaluatedPoints() { _evaluatedPoints.clear(); }
    size_t getNbEvaluatedPoints() const { return _evaluatedPoints.size(); }

    const SuccessType& getSuccessType() const { return _success; }
    void setSuccessType(const SuccessType& success);

    size_t getCurrentlyRunning() const { return _currentlyRunning; }
    void incCurrentlyRunning(const size_t n);
    void decCurrentlyRunning(const size_t n);

    void setComputeSuccessTypeFunction(const ComputeSuccessFunction &computeSuccessFunction) { _computeSuccessType.setComputeSuccessTypeFunction(computeSuccessFunction); }
    ComputeSuccessType getComputeSuccessType() const { return _computeSuccessType; }

    void setStopReason(const EvalStopType& s);
    const StopReason<EvalStopType>& getStopReason() const { return _stopReason; }
    std::string getStopReasonAsString() const { return _stopReason.getStopReasonAsString(); }
    bool testIf(const EvalStopType& s) const { return s == _stopReason.get(); }
    bool checkEvalTerminate() const { return _stopReason.checkTerminate(); }

private:
    /// Helper for constructor
    void init();
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVCMAINTHREADINFO__
