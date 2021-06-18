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
 \file   EvcMainThreadInfo.hpp
 \brief  Information about main threads
 \author Viviane Rochon Montplaisir
 \date   June 2020
 \see    EvcMainThreadInfo.cpp
 */

#ifndef __NOMAD_4_0_EVCMAINTHREADINFO__
#define __NOMAD_4_0_EVCMAINTHREADINFO__

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
    std::shared_ptr<Evaluator>      _evaluator;         ///< The Evaluator for either blackbox or model evaluations.
    const std::unique_ptr<EvaluatorControlParameters> _evalContParams;  ///< The parameters controlling the behavior of EvaluatorControl for this main thread
    std::atomic<size_t>             _nbPointsInQueue;   ///< Number of points in the evaluation queue for this main thread
    bool                            _doneWithEval;      ///< All evaluations done for this main thread
    std::shared_ptr<Barrier>        _barrier;
    std::shared_ptr<EvalPoint>      _bestIncumbent;     ///< Temporary value useful for display only
    std::vector<EvalPoint>          _evaluatedPoints;   ///< Where evaluated points are put temporarily
    SuccessType                     _success;           ///< Success type of the last run
    std::atomic<size_t>             _currentlyRunning;  ///< Count number of evaluations currently running.
    size_t                          _lapMaxBbEval;      ///< The maximum number of blackbox evaluations that can be performed by a sub algorithm.
    std::atomic<size_t>             _lapBbEval;         ///< The number of blackbox evaluations performed by a given sub algorithm (reset at Algorithm start).
    std::atomic<size_t>             _modelEval;         ///< The number of quad or sgtelib model evaluations performed (since last reset)
    std::atomic<size_t>             _subBbEval;         ///< The number of bb eval for a subproblem (e.g. SSD-Mads context)
    ComputeType                     _computeType;       ///< How to compute f and h
    std::shared_ptr<Direction>      _lastSuccessfulFeasDir; ///< Direction of last success for feasible points. May be used to sort points before evaluation.
    std::shared_ptr<Direction>      _lastSuccessfulInfDir; ///< Direction of last success for infeasible points. May be used to sort points before evaluation.
    StopReason<EvalMainThreadStopType> _stopReason;

public:
    /// Constructor
    /**
     \param evaluator       The Evaluator for either blackbox or model evaluations-- \b IN.
     \param evalContParams  The parameters controlling how the EvaluatorControl behaves for this main thread-- \b IN.
     */
    explicit EvcMainThreadInfo(std::shared_ptr<Evaluator> evaluator,
                               std::unique_ptr<EvaluatorControlParameters> evalContParams)
      : _evaluator(evaluator),
        _evalContParams(std::move(evalContParams)),
        _nbPointsInQueue(0),
        _doneWithEval(false),
        _barrier(),
        _bestIncumbent(),
        _evaluatedPoints(),
        _success(SuccessType::UNSUCCESSFUL),
        _currentlyRunning(0),
        _lapMaxBbEval(INF_SIZE_T),
        _lapBbEval(0),
        _modelEval(0),
        _subBbEval(0),
        _computeType(ComputeType::STANDARD),
        _lastSuccessfulFeasDir(nullptr),
        _lastSuccessfulInfDir(nullptr),
        _stopReason()
    {
        init();
    }

    /// Set Evaluator and return old Evaluator.
    /**
     \param evaluator       The Evaluator for either blackbox or model evaluations-- \b IN.
     \return                The previous Evaluator.
     */
    std::shared_ptr<Evaluator> setEvaluator(std::shared_ptr<Evaluator> evaluator);
    const Evaluator* getEvaluator() { return _evaluator.get(); }
    std::shared_ptr<EvalParameters> getEvalParams() const;
    EvalType getEvalType() const;

    size_t getNbPointsInQueue() const { return _nbPointsInQueue; }
    void incNbPointsInQueue();
    void decNbPointsInQueue();
    void resetNbPointsInQueue() { _nbPointsInQueue = 0; }

    bool getDoneWithEval() const { return _doneWithEval; }
    void setDoneWithEval(const bool doneWithEval) { _doneWithEval = doneWithEval; }

    // Get and set parameters
    bool getOpportunisticEval() const;
    void setOpportunisticEval(const bool opportunisticEval);
    bool getUseCache() const;
    void setUseCache(const bool useCache);
    size_t getMaxBbEvalInSubproblem() const;
    void setMaxBbEvalInSubproblem(const size_t maxBbEval);
    bool getSurrogateOptimization() const;
    void setSurrogateOptimization(const bool surrogateOptimization);

    // Get and set counters
    // Resetting a counter also resets stop reason, if the stop reason was that the max was reached for this counter.
    void setLapMaxBbEval(const size_t maxBbEval) { _lapMaxBbEval = maxBbEval; }
    size_t getLapBbEval() const { return _lapBbEval; }
    size_t getLapMaxBbEval() const { return _lapMaxBbEval; }
    void incLapBbEval(const size_t countEval) { _lapBbEval += countEval; }
    void resetLapBbEval();
    size_t getModelEval() const { return _modelEval; }
    void resetModelEval();
    void incModelEval(const size_t countEval) { _modelEval += countEval; }
    size_t getBbEvalInSubproblem() const { return _subBbEval; }
    void incBbEvalInSubproblem(const size_t countEval) { _subBbEval += countEval; }
    void resetBbEvalInSubproblem();

    const std::shared_ptr<Barrier>& getBarrier() const { return _barrier; }
    void setBarrier(const std::shared_ptr<Barrier>& barrier) { _barrier = barrier; }

    const std::shared_ptr<EvalPoint>& getBestIncumbent() const;
    void setBestIncumbent(const std::shared_ptr<EvalPoint>& bestIncumbent);

    std::vector<EvalPoint> retrieveAllEvaluatedPoints();
    void addEvaluatedPoint(const EvalPoint& evaluatedPoint);
    bool remainsEvaluatedPoints() const { return !_evaluatedPoints.empty(); }
    void clearEvaluatedPoints() { _evaluatedPoints.clear(); }
    size_t getNbEvaluatedPoints() const { return _evaluatedPoints.size(); }

    const SuccessType& getSuccessType() const { return _success; }
    void setSuccessType(const SuccessType& success);

    size_t getCurrentlyRunning() const { return _currentlyRunning; }
    void incCurrentlyRunning();
    void decCurrentlyRunning();

    void setComputeType(const ComputeType &computeType) { _computeType = computeType; }
    const ComputeType& getComputeType() const { return _computeType; }

    void setLastSuccessfulFeasDir(const std::shared_ptr<Direction> &feasDir) { _lastSuccessfulFeasDir = feasDir; }
    const std::shared_ptr<Direction>& getLastSuccessfulFeasDir() const { return _lastSuccessfulFeasDir; }
    void setLastSuccessfulInfDir(const std::shared_ptr<Direction> &feasDir) { _lastSuccessfulInfDir = feasDir; }
    const std::shared_ptr<Direction>& getLastSuccessfulInfDir() const { return _lastSuccessfulInfDir; }

    void setStopReason(const EvalMainThreadStopType& s);
    const StopReason<EvalMainThreadStopType>& getStopReason() const { return _stopReason; }
    std::string getStopReasonAsString() const { return _stopReason.getStopReasonAsString(); }
    bool testIf(const EvalMainThreadStopType& s) const { return s == _stopReason.get(); }
    bool checkEvalTerminate() const { return _stopReason.checkTerminate(); }

private:
    /// Helper for constructor
    void init();
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_EVCMAINTHREADINFO__
