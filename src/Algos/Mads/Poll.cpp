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

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/Mads/DoublePollMethod.hpp"
#include "../../Algos/Mads/NP1UniPollMethod.hpp"
#include "../../Algos/Mads/Ortho2NPollMethod.hpp"
#include "../../Algos/Mads/OrthoNPlus1PollMethod.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Algos/Mads/QR2NPollMethod.hpp"
#include "../../Algos/Mads/SinglePollMethod.hpp"
#include "../../Algos/Mads/UserPollMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"

#ifdef TIME_STATS
#include "../../Util/Clock.hpp"
// Initialize static variables
double NOMAD::Poll::_pollTime;
double NOMAD::Poll::_pollEvalTime;
#endif // TIME_STATS

/// <#Description#>
void NOMAD::Poll::init()
{
    setStepType(NOMAD::StepType::POLL);
    verifyParentNotNull();

    _trialPointMaxAddUp = 0;

    _hasSecondPass = false;
    _hasUserPollMethod = false;
    _hasUserFreePollMethod = false;
    if ( nullptr != _runParams)
    {
        // The direction types for primary and secondary poll centers
        _primaryDirectionTypes = _runParams->getAttributeValue<DirectionTypeList>("DIRECTION_TYPE");
        _secondaryDirectionTypes = _runParams->getAttributeValue<DirectionTypeList>("DIRECTION_TYPE_SECONDARY_POLL");

        // Ortho n+1 poll methods generate n trial points in a first pass and, if not successful, generate the n+1 th point (second pass)
        for (auto dirType : _primaryDirectionTypes)
        {
            if (NOMAD::DirectionType::ORTHO_NP1_NEG == dirType
                || NOMAD::DirectionType::ORTHO_NP1_QUAD == dirType)
            {
                _hasSecondPass = true;
            }
            if (NOMAD::DirectionType::USER_POLL == dirType)
            {
                _hasUserPollMethod = true;
            }
            if (NOMAD::DirectionType::USER_FREE_POLL == dirType)
            {
                _hasUserFreePollMethod = true;
            }
        }
        if (_hasUserPollMethod && _hasUserFreePollMethod)
        {
            throw NOMAD::Exception(__FILE__, __LINE__,"Cannot have both user regular poll and free poll methods.");
        }

        // Rho parameter of the progressive barrier. Used to choose if the primary frame center is the feasible or infeasible incumbent.
        _rho = _runParams->getAttributeValue<NOMAD::Double>("RHO");

        // Complete the primary and secondary poll directions to reach the given target number. Only for single pass direction type (ortho 2n). Managed by checkAndComply.
        _trialPointMaxAddUp = _runParams->getAttributeValue<size_t>("TRIAL_POINT_MAX_ADD_UP");

    }

    // Groups of variables.
    if ( nullptr != _pbParams)
    {
        _varGroups = _pbParams->getAttributeValue<NOMAD::ListOfVariableGroup>("VARIABLE_GROUP");
        if (!_varGroups.empty())
        {
            _mapDirTypeToVG = _runParams->getMapDirTypeToVG();
        }
    }

    // Unlike for search methods, we cannot generate poll methods during init because they depend on primary and secondary poll centers. The init is called once when instantiating Mads poll and poll centers change.


}


void NOMAD::Poll::startImp()
{
    // Sanity check.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    // Reset the current counters. The total counters are not reset (done only once when constructor is called).
    _trialPointStats.resetCurrentStats();

}


bool NOMAD::Poll::runImp()
{
    bool pollSuccessful = false;
    std::string s;

    // Sanity check. The runImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    // Go through all poll methods to generate points.
    OUTPUT_DEBUG_START
    s = "Generate points for " + getName();
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
#ifdef TIME_STATS
    double pollStartTime = NOMAD::Clock::getCPUTime();
    double pollEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS

    // 1- Create poll methods and generate points (first pass)
    generateTrialPoints();

    if ( ! _stopReasons->checkTerminate() )
    {
        // 1b- Count the points that would need eval (check cache and barrier)
        countTrialPointsThatNeedEval(this);

        // 1c- Generate extra trial points (using SINGLE direction) to reach a given number (if provided). Only for single pass direction type (ortho 2n).
        generateTrialPointsExtra();

        // 2- Complete trial points information (first pass)
        completeTrialPointsInformation();

        // 3- Evaluate points (first pass)
        evalTrialPoints(this);
        pollSuccessful = (_success >= NOMAD::SuccessType::FULL_SUCCESS);

    }

    // Second pass: only for Ortho N+1 and no success
    if (!_stopReasons->checkTerminate() && !pollSuccessful && _hasSecondPass)
    {
        // Generate points for the second pass. Case of the N+1 th point for Ortho N+1 methods.
        // Erase existing trial points to ensure that they
        // are not re-evaluated, then put them back for postProcessing().
        const auto& evaluatedTrialPoints = getTrialPoints();
        clearTrialPoints();

        // 1- Create poll methods and generate points (second pass)
        generateTrialPointsSecondPass();

        // 2- Complete trial points information (second pass)
        completeTrialPointsInformation();

        // Evaluate point.
        if (getTrialPointsCount() > 0)
        {
            evalTrialPoints(this);
            pollSuccessful = (_success >= NOMAD::SuccessType::FULL_SUCCESS);
        }

        // Add back trial points that are already evaluated.
        for (const auto& trialPoint: evaluatedTrialPoints)
        {
            insertTrialPoint(trialPoint);
        }
    }

    // Special pass for user poll method only
    if (!_stopReasons->checkTerminate() && _userPollMethodCallbackEnabled)
    {
        // Update at the end. Implemented only by user poll method.
        for (const auto& pollMethod : _pollMethods)
        {
            pollMethod->updateEndUserPoll();
        }
    }

#ifdef TIME_STATS
    _pollTime += NOMAD::Clock::getCPUTime() - pollStartTime;
    _pollEvalTime += NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime() - pollEvalStartTime;
#endif // TIME_STATS

    OUTPUT_INFO_START
    s = getName();
    s += (pollSuccessful) ? " is successful" : " is not successful";
    s += ". Stop reason: ";
    s += _stopReasons->getStopReasonAsString() ;
    AddOutputInfo(s);
    OUTPUT_INFO_END

    return pollSuccessful;
}


void NOMAD::Poll::endImp()
{
    // Sanity check. The endImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    // Compute hMax and update Barrier.
    postProcessing();

}


// Generate new points to evaluate
// Note: Used whether MEGA_SEARCH_POLL is true or false.
void NOMAD::Poll::generateTrialPointsImp()
{

    // Poll methods depend on poll centers
    createPollMethodsForPollCenters();

    // Use poll methods to create trial points
    for (const auto& pollMethod : _pollMethods)
    {
        if (_stopReasons->checkTerminate())
        {
            break;
        }

        // Generate trial points. Snap to bounds and project on mesh is also done.
        pollMethod->generateTrialPoints();

        // Reduce the number of trial points.
        // Currently, this is used only by Ortho N+1
        pollMethod->trialPointsReduction();

        // Pass the points from Poll method to Poll for later evaluation
        const auto& pollMethodPoints = pollMethod->getTrialPoints();
        for (const auto & point : pollMethodPoints)
        {
            insertTrialPoint(point);
        }
    }

    // Stopping criterion
    if (0 == getTrialPointsCount())
    {
        // Success type when no trial points are produced
        _success = NOMAD::SuccessType::NO_TRIALS;

        setMeshPrecisionStopType();
    }
}

void NOMAD::Poll::setMeshPrecisionStopType()
{
    auto madsStopReasons = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get(_stopReasons);
    madsStopReasons->set(NOMAD::MadsStopType::MESH_PREC_REACHED);
}

void NOMAD::Poll::generateTrialPointsSecondPass()
{
    for (const auto& pollMethod : _pollMethods)
    {
        if (_stopReasons->checkTerminate())
        {
            break;
        }
        if (!pollMethod->hasSecondPass())
        {
            continue;
        }

        // Second pass (optional) generation of trial points. Case for Ortho n+1 methods.
        pollMethod->generateTrialPointsSecondPass();

        // Add the additional poll method trial points to poll trial points.
        for (const auto& point : pollMethod->getTrialPoints())
        {
            insertTrialPoint(point);
        }
    }
}

void NOMAD::Poll::generateTrialPointsExtra()
{
    size_t tmpTrial = _nbEvalPointsThatNeedEval;
    const size_t maxLoops = 100;
    size_t loops =0;
    while (tmpTrial < _trialPointMaxAddUp && loops < maxLoops)
    {
        loops++;
        for (const auto & frameCenter : _frameCenters)
        {
            auto pollMethod = std::make_shared<NOMAD::SinglePollMethod>(this, frameCenter);

            pollMethod->generateTrialPoints();

            // Add the additional poll method trial points to poll trial points.
            for (const auto& point : pollMethod->getTrialPoints())
            {
                insertTrialPoint(point);
                tmpTrial ++;
            }
        }
        if (tmpTrial >= _trialPointMaxAddUp )
        {
            break;
        }
    }
}



void NOMAD::Poll::computePrimarySecondaryPollCenters(
                        std::vector<NOMAD::EvalPointPtr> &primaryCenters,
                        std::vector<NOMAD::EvalPointPtr> &secondaryCenters) const
{
    auto barrier = getMegaIterationBarrier();

    if (nullptr != barrier)
    {
        auto firstXIncFeas = barrier->getCurrentIncumbentFeas();
        auto firstXIncInf  = barrier->getCurrentIncumbentInf();
        bool primaryIsInf = false;

        // Negative rho means make no distinction between primary and secondary polls.
        const bool usePrimarySecondary = (_rho >= 0) && (nullptr != firstXIncFeas) && (nullptr != firstXIncInf);
        if (usePrimarySecondary)
        {
            const auto computeType = barrier->getFHComputeType();

            const NOMAD::ArrayOfDouble fFeas = firstXIncFeas->getFs(computeType);
            const NOMAD::ArrayOfDouble fInf  = firstXIncFeas->getFs(computeType);

            if (fFeas.isComplete() && fInf.isComplete())
            {
                // Single objective case
                if (fFeas.size() == 1 && (fFeas[0] - _rho) > fInf[0])
                {
                    // xFeas' f is too large, use xInf as primary poll instead.
                    primaryIsInf = true;
                }
                // Multiobjective case
                else if (fFeas.size() > 1)
                {
                    // 1- Compute dominance move value Psi (report to DMultiMADS-PB article).
                    NOMAD::Double domMove = NOMAD::INF;
                    for (const auto& xFeas: barrier->getAllXFeas())
                    {
                        // It means the current infeasible frame center is dominated by the
                        // previous element of the set of feasible points of the barrier.
                        // It will never be the primary frame center.
                        // No need to continue.
                        if (domMove == 0.0)
                            break;

                        NOMAD::Double currentDomMove = 0.0;
                        const NOMAD::ArrayOfDouble currentFFeas = xFeas->getFs(computeType);
                        for (size_t obj = 0; obj < fFeas.size(); ++obj)
                        {
                            currentDomMove += std::max(currentFFeas[obj].todouble() - fInf[obj].todouble(), 0.0);
                        }
                        domMove = std::min(domMove, currentDomMove);
                    }

                    // 2- Compute mu value (report to DMulti-MADS PB article).
                    if (domMove > 0)
                    {
                        NOMAD::Double mu = 0.0;
                        for (size_t obj = 0; obj < fFeas.size(); ++obj)
                        {
                            NOMAD::Double fMin = NOMAD::INF;
                            NOMAD::Double fMax = NOMAD::M_INF;
                            for (const auto& xFeas: barrier->getAllXFeas())
                            {
                                const NOMAD::Double xFeasFvalue = xFeas->getFs(computeType)[obj];
                                fMin = std::min(xFeasFvalue, fMin);
                                fMax = std::max(xFeasFvalue, fMin);
                            }
                            mu += (fMin == fMax) ? std::fabs(fMin.todouble())
                                                 : std::fabs(fMax.todouble() - fMin.todouble());
                        }

                        // xInf' f is ``far'' into the dominance zone of the set of feasible solutions,
                        // use it as primary poll instead.
                        if ((domMove - _rho * mu) > 0)
                            primaryIsInf = true;
                    }
                }
            }
        }

        if (usePrimarySecondary)
        {
            if (primaryIsInf)
            {
                primaryCenters.push_back( firstXIncInf );
                secondaryCenters.push_back(firstXIncFeas );
            }
            else
            {
                primaryCenters.push_back(firstXIncFeas);
                secondaryCenters.push_back(firstXIncInf);
            }
        }
        else
        {
            if (nullptr != firstXIncFeas)
            {
                primaryCenters.push_back(firstXIncFeas);
            }
            else
            {
                if (nullptr != firstXIncInf)
                {
                    primaryCenters.push_back(firstXIncInf);
                }
            }
        }
    }
}

void NOMAD::Poll::createPollMethods(const bool isPrimary, const EvalPointPtr& frameCenter)
{
    // Select the poll methods to be executed
    NOMAD::DirectionTypeList dirTypes = (isPrimary) ? _primaryDirectionTypes : _secondaryDirectionTypes;

    for (auto dirType : dirTypes)
    {
        std::shared_ptr<NOMAD::PollMethodBase> pollMethod;
        _frameCenters.push_back(frameCenter);
        switch (dirType)
        {
            case DirectionType::ORTHO_2N:
                pollMethod = std::make_shared<NOMAD::Ortho2NPollMethod>(this, frameCenter);
                break;
            case DirectionType::ORTHO_NP1_NEG:
                pollMethod =std::make_shared<NOMAD::OrthoNPlus1PollMethod>(this, frameCenter,false /* disable quad model opt. */);
                break;
            case DirectionType::ORTHO_NP1_QUAD:
                pollMethod = std::make_shared<NOMAD::OrthoNPlus1PollMethod>(this, frameCenter,true /* enable quad model opt. */);
                break;
            case DirectionType::QR_2N:
                pollMethod = std::make_shared<NOMAD::QR2NPollMethod>(this, frameCenter);
                break;
            case DirectionType::NP1_UNI:
                pollMethod = std::make_shared<NOMAD::NP1UniPollMethod>(this, frameCenter);
                break;
            case DirectionType::DOUBLE:
                pollMethod = std::make_shared<NOMAD::DoublePollMethod>(this, frameCenter);
                break;
            case DirectionType::SINGLE:
                pollMethod = std::make_shared<NOMAD::SinglePollMethod>(this, frameCenter);
                break;
            case DirectionType::USER_POLL:
                if ( !_userPollMethodCallbackEnabled )
                {
                    throw NOMAD::Exception(__FILE__, __LINE__,"Cannot add user poll method.");
                }
                pollMethod = std::make_shared<NOMAD::UserPollMethod>(this, frameCenter, false /* not user free poll*/);
                break;
            case DirectionType::USER_FREE_POLL:
                if ( !_userPollMethodCallbackEnabled )
                {
                    throw NOMAD::Exception(__FILE__, __LINE__,"Cannot add user free poll method.");
                }
                pollMethod = std::make_shared<NOMAD::UserPollMethod>(this, frameCenter, true /*a free user poll method*/);
                break;
            default:
                throw NOMAD::Exception(__FILE__, __LINE__,"Poll method " + directionTypeToString(dirType) + " is not available.");
                break;
        }

        // A variable group can be specifically managed by a direction type.
        std::map<NOMAD::DirectionType,NOMAD::ListOfVariableGroup>::const_iterator it;
        if (!_mapDirTypeToVG.empty() && (it = _mapDirTypeToVG.find(dirType)) != _mapDirTypeToVG.end())
        {
            pollMethod->setListVariableGroups(it->second);
        }
        // A default initialization with all variable groups is done for each direction type
        else
        {

            pollMethod->setListVariableGroups(_varGroups);
        }

        pollMethod->setIsPrimary(isPrimary);

        pollMethod->setIsPrimary(isPrimary);

        pollMethod->setIsPrimary(isPrimary);

        _pollMethods.push_back(pollMethod);
    }
    if (_userPollMethodCallbackEnabled && (_hasUserPollMethod || _hasUserFreePollMethod) && _hasSecondPass)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"User-Defined poll method cannot be combined with a poll method that works in two passes such as ORTHO N+1.");
    }
    if (_userPollMethodCallbackEnabled && _hasUserFreePollMethod && _pollMethods.size() <= 1)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"User-Defined free poll method must be combined with a regular poll method, except ORTHO N+1.");
    }
}


void NOMAD::Poll::createPollMethodsForPollCenters()
{
    // Compute primary and secondary poll centers
    std::vector<NOMAD::EvalPointPtr> primaryCenters, secondaryCenters;
    computePrimarySecondaryPollCenters(primaryCenters, secondaryCenters);

    // Add poll methods for primary polls
    _pollMethods.clear();
    _frameCenters.clear();
    for (const auto & pollCenter : primaryCenters)
    {
        createPollMethods(true, pollCenter);
    }
    // Add poll methods for secondary polls
    for (const auto & pollCenter : secondaryCenters)
    {
        createPollMethods(false, pollCenter);
    }
}
