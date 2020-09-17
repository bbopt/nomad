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
