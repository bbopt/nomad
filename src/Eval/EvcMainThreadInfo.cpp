/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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

#include "../Eval/EvcMainThreadInfo.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Util/MicroSleep.hpp"

/*-------------------------*/
/* Class EvcMainThreadInfo */
/*-------------------------*/

void NOMAD::EvcMainThreadInfo::init()
{

    _evalSortType = _evalContParams->getAttributeValue<NOMAD::EvalSortType>("EVAL_QUEUE_SORT");
    
    _evalOpportunistic = _evalContParams->getAttributeValue<bool>("EVAL_OPPORTUNISTIC");
    _useCache = _evalContParams->getAttributeValue<bool>("EVAL_USE_CACHE");
    _subPbMaxBBEval = _evalContParams->getAttributeValue<size_t>("SUBPROBLEM_MAX_BB_EVAL");
    _evalSurrogateOptimization = _evalContParams->getAttributeValue<bool>("EVAL_SURROGATE_OPTIMIZATION");
}


bool NOMAD::EvcMainThreadInfo::hasEvaluator(NOMAD::EvalType evalType) const
{
    if (_evaluators.empty())
    {
        return false;
    }
    
    auto it = std::find_if(_evaluators.begin(),_evaluators.end(), [evalType](const NOMAD::EvaluatorPtr& e){ return e->getEvalType() == evalType; });
    
    if ( _evaluators.end() == it )
    {
        return false;
    }
    return true;
}

const NOMAD::Evaluator* NOMAD::EvcMainThreadInfo::getCurrentEvaluator() const
{
    if (_evaluators.empty())
    {
        std::string s = "Error in EvaluatorControl main thread management: no evaluator is registered.";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
    if (NOMAD::EvalType::UNDEFINED == _currentEvaluatorType)
    {
        std::string s = "Error in EvaluatorControl main thread management: current evaluator type is undefined.";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
    
    auto it = std::find_if(_evaluators.begin(),_evaluators.end(), [evalType=_currentEvaluatorType](const NOMAD::EvaluatorPtr& e){ return e->getEvalType() == evalType; });
    
    if ( _evaluators.end() == it )
    {
        std::string s = "Error in EvaluatorControl main thread management: evaluator with EvalType = " + NOMAD::evalTypeToString(_currentEvaluatorType);
        s += " is not available";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
    return (it->get());
}

void NOMAD::EvcMainThreadInfo::addEvaluator(const NOMAD::EvaluatorPtr& evaluator)
{
    if ( nullptr == evaluator )
    {
        std::string s = "Error in EvaluatorControl main thread management: cannot assign nullptr.";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
    
    auto evalType = evaluator->getEvalType();
    
    auto it = std::find_if(_evaluators.begin(),_evaluators.end(), [evalType](const NOMAD::EvaluatorPtr& e){ return e->getEvalType() == evalType; });
    
    if ( _evaluators.end() != it )
    {
        _evaluators.erase(it);
    }
    _evaluators.push_back(evaluator);
    _currentEvaluatorType = evaluator->getEvalType(); // The last added evaluator becomes the current evaluator
}

void NOMAD::EvcMainThreadInfo::incNbPointsInQueue()
{
    _nbPointsInQueue++;
}


void NOMAD::EvcMainThreadInfo::decNbPointsInQueue()
{
    if (0 == _nbPointsInQueue)
    {
        std::string s = "Error in EvaluatorControl main thread management: Trying to decrease number of points in queue which is already 0";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
    _nbPointsInQueue--;
}

void NOMAD::EvcMainThreadInfo::setEvalSortType(EvalSortType evalSortType)
{
    _evalSortType = evalSortType;
    _evalContParams->setAttributeValue("EVAL_QUEUE_SORT", evalSortType);
    _evalContParams->checkAndComply();
}


void NOMAD::EvcMainThreadInfo::setOpportunisticEval(const bool opportunisticEval)
{
    _evalOpportunistic = opportunisticEval;
    _evalContParams->setAttributeValue("EVAL_OPPORTUNISTIC", opportunisticEval);
    _evalContParams->checkAndComply();
}


void NOMAD::EvcMainThreadInfo::setUseCache(const bool useCache)
{
    _useCache = useCache;
    _evalContParams->setAttributeValue("EVAL_USE_CACHE", useCache);
    _evalContParams->checkAndComply();
}

void NOMAD::EvcMainThreadInfo::setMaxBbEvalInSubproblem(const size_t maxBbEval)
{
    _subPbMaxBBEval = maxBbEval;
    _evalContParams->setAttributeValue("SUBPROBLEM_MAX_BB_EVAL", maxBbEval);
    _evalContParams->checkAndComply();
}



void NOMAD::EvcMainThreadInfo::setSurrogateOptimization(const bool surrogateOptimization)
{
    _evalSurrogateOptimization = surrogateOptimization;
    _evalContParams->setAttributeValue("EVAL_SURROGATE_OPTIMIZATION", surrogateOptimization);
    _evalContParams->checkAndComply();
}


void NOMAD::EvcMainThreadInfo::resetLapBbEval()
{
    _lapBbEval = 0;
    if (NOMAD::EvalMainThreadStopType::LAP_MAX_BB_EVAL_REACHED == _stopReason.get())
    {
        _stopReason.set(NOMAD::EvalMainThreadStopType::STARTED);
    }
}


void NOMAD::EvcMainThreadInfo::resetModelEval()
{
    _modelEval = 0;
    if (NOMAD::EvalMainThreadStopType::MAX_MODEL_EVAL_REACHED == _stopReason.get())
    {
        _stopReason.set(NOMAD::EvalMainThreadStopType::STARTED);
    }
}


void NOMAD::EvcMainThreadInfo::resetBbEvalInSubproblem()
{
    _subBbEval = 0;
    if (NOMAD::EvalMainThreadStopType::SUBPROBLEM_MAX_BB_EVAL_REACHED == _stopReason.get())
    {
        _stopReason.set(NOMAD::EvalMainThreadStopType::STARTED);
    }
}


NOMAD::EvalPointPtr NOMAD::EvcMainThreadInfo::getBestIncumbent() const
{
    return _bestIncumbent;
}


void NOMAD::EvcMainThreadInfo::setBestIncumbent(const NOMAD::EvalPointPtr& bestIncumbent)
{
        _bestIncumbent = bestIncumbent;
}

std::vector<NOMAD::EvalPoint> NOMAD::EvcMainThreadInfo::retrieveAllEvaluatedPoints()
{
    std::vector<NOMAD::EvalPoint> allEvaluatedPoints;

    bool warningShown = false;
    while (_currentlyRunning > 0)
    {
        // Due to race conditions, retrieveAllEvaluatedPoints() might be called while
        // some points are still being evaluated. Wait for all points to be evaluated.
        OUTPUT_INFO_START
        if (!warningShown)
        {
            std::string s = "Warning: Calling retrieveAllEvaluatedPoints() while still ";
            s += NOMAD::itos(_currentlyRunning) + " currently running";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
            warningShown = true;
        }
        OUTPUT_INFO_END
        usleep(10);
    }

    allEvaluatedPoints.insert(allEvaluatedPoints.end(),
                              std::make_move_iterator(_evaluatedPoints.begin()),
                              std::make_move_iterator(_evaluatedPoints.end()));
    _evaluatedPoints.clear();

    return allEvaluatedPoints;
}


void NOMAD::EvcMainThreadInfo::addEvaluatedPoint(const NOMAD::EvalPoint& evaluatedPoint)
{
#ifdef _OPENMP
    #pragma omp critical(addEvaluatedPoint)
#endif // _OPENMP
    {
        _evaluatedPoints.push_back(evaluatedPoint);
    }
}


void NOMAD::EvcMainThreadInfo::setSuccessType(const NOMAD::SuccessType& success)
{
    // Note setSuccessType is already called inside a critical section, so do not
    // add another critical section here.
    _success = success;
}


void NOMAD::EvcMainThreadInfo::incCurrentlyRunning(size_t k)
{
    _currentlyRunning += k;
}


void NOMAD::EvcMainThreadInfo::decCurrentlyRunning(size_t k)
{
    if (0 == _currentlyRunning)
    {
        OUTPUT_DEBUG_START
        std::string s = "Trying to decrease number of currently running evaluations which is already 0. Evaluation results probably come from previous run cache file. This is not an error.";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
    }
    else
    {
        _currentlyRunning -=k;
    }
}


void NOMAD::EvcMainThreadInfo::setStopReason(const NOMAD::EvalMainThreadStopType& s)
{
    _stopReason.set(s);
}
