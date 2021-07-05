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

#include "../Eval/EvcMainThreadInfo.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Util/MicroSleep.hpp"

/*-------------------------*/
/* Class EvcMainThreadInfo */
/*-------------------------*/
void NOMAD::EvcMainThreadInfo::init()
{
}


std::shared_ptr<NOMAD::Evaluator> NOMAD::EvcMainThreadInfo::setEvaluator(std::shared_ptr<NOMAD::Evaluator> evaluator)
{
    auto previousEvaluator = _evaluator;
    _evaluator = evaluator;

    return previousEvaluator;
}


std::shared_ptr<NOMAD::EvalParameters> NOMAD::EvcMainThreadInfo::getEvalParams() const
{
    return (nullptr == _evaluator) ? nullptr : _evaluator->getEvalParams();
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


NOMAD::EvalType NOMAD::EvcMainThreadInfo::getEvalType() const
{
    return (nullptr == _evaluator) ? NOMAD::EvalType::UNDEFINED : _evaluator->getEvalType();
}


bool NOMAD::EvcMainThreadInfo::getOpportunisticEval() const
{
    while (true)
    {
        try
        {
            return _evalContParams->getAttributeValue<bool>("EVAL_OPPORTUNISTIC");
        }
        catch (NOMAD::ParameterToBeChecked&)
        {
            // Exception due to parameters being in process of checkAndComply().
            // While will loop - Retry
        }
    }
}


void NOMAD::EvcMainThreadInfo::setOpportunisticEval(const bool opportunisticEval)
{
    _evalContParams->setAttributeValue("EVAL_OPPORTUNISTIC", opportunisticEval);
    _evalContParams->checkAndComply();
}


bool NOMAD::EvcMainThreadInfo::getUseCache() const
{
    while (true)
    {
        try
        {
            return _evalContParams->getAttributeValue<bool>("EVAL_USE_CACHE");
        }
        catch (NOMAD::ParameterToBeChecked&)
        {
            // Exception due to parameters being in process of checkAndComply().
            // While will loop - Retry
        }
    }
}


void NOMAD::EvcMainThreadInfo::setUseCache(const bool useCache)
{
    _evalContParams->setAttributeValue("EVAL_USE_CACHE", useCache);
    _evalContParams->checkAndComply();
}


size_t NOMAD::EvcMainThreadInfo::getMaxBbEvalInSubproblem() const
{
    while (true)
    {
        try
        {
            return _evalContParams->getAttributeValue<size_t>("SUBPROBLEM_MAX_BB_EVAL");
        }
        catch (NOMAD::ParameterToBeChecked&)
        {
            // Exception due to parameters being in process of checkAndComply().
            // While will loop - Retry
        }
    }
}


void NOMAD::EvcMainThreadInfo::setMaxBbEvalInSubproblem(const size_t maxBbEval)
{
    _evalContParams->setAttributeValue("SUBPROBLEM_MAX_BB_EVAL", maxBbEval);
    _evalContParams->checkAndComply();
}


bool NOMAD::EvcMainThreadInfo::getSurrogateOptimization() const
{
    while (true)
    {
        try
        {
            return _evalContParams->getAttributeValue<bool>("EVAL_SURROGATE_OPTIMIZATION");
        }
        catch (NOMAD::ParameterToBeChecked&)
        {
            // Exception due to parameters being in process of checkAndComply().
            // While will loop - Retry
        }
    }
}


void NOMAD::EvcMainThreadInfo::setSurrogateOptimization(const bool surrogateOptimization)
{
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


const std::shared_ptr<NOMAD::EvalPoint>& NOMAD::EvcMainThreadInfo::getBestIncumbent() const
{
    return _bestIncumbent;
}


void NOMAD::EvcMainThreadInfo::setBestIncumbent(const std::shared_ptr<NOMAD::EvalPoint>& bestIncumbent)
{
    NOMAD::ComputeSuccessType computeSuccess(_evaluator->getEvalType(), _computeType);
    if (computeSuccess(bestIncumbent, _bestIncumbent) >= NOMAD::SuccessType::PARTIAL_SUCCESS)
    {
        _bestIncumbent = bestIncumbent;
    }
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


void NOMAD::EvcMainThreadInfo::incCurrentlyRunning()
{
    _currentlyRunning ++;
}


void NOMAD::EvcMainThreadInfo::decCurrentlyRunning()
{
    if (0 == _currentlyRunning)
    {
        // Note: we don't know the main thread number.
        std::string s = "Error in EvaluatorControl main thread management: Trying to decrease number of currently running evaluations which is already 0";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
    _currentlyRunning --;
}


void NOMAD::EvcMainThreadInfo::setStopReason(const NOMAD::EvalMainThreadStopType& s)
{
    _stopReason.set(s);
}
