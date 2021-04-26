/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0 has been created by                                        */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0 is owned by                               */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,            */
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
#include "../../Algos/Mads/DoublePollMethod.hpp"
#include "../../Algos/Mads/NP1UniPollMethod.hpp"
#include "../../Algos/Mads/Ortho2NPollMethod.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Algos/Mads/SinglePollMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"

#ifdef TIME_STATS
#include "../../Util/Clock.hpp"
// Initialize static variables
double NOMAD::Poll::_pollTime;
double NOMAD::Poll::_pollEvalTime;
#endif // TIME_STATS

void NOMAD::Poll::init()
{
    _name = "Poll";
    verifyParentNotNull();

    // Compute primary and secondary poll centers
    std::vector<NOMAD::EvalPoint> primaryCenters, secondaryCenters;
    computePrimarySecondaryPollCenters(primaryCenters, secondaryCenters);

    // Add poll methods for primary polls
    for (auto pollCenter : primaryCenters)
    {
        auto pollMethod = createPollMethod(true, pollCenter);
        _pollMethods.push_back(pollMethod);
    }
    // Add poll methods for secondary polls
    for (auto pollCenter : secondaryCenters)
    {
        auto pollMethod = createPollMethod(false, pollCenter);
        _pollMethods.push_back(pollMethod);
    }
}


void NOMAD::Poll::startImp()
{
    // Sanity check.
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);
}


bool NOMAD::Poll::runImp()
{
    bool pollSuccessful = false;
    std::string s;

    // Sanity check. The runImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    // Go through all poll methods to generate points.
    OUTPUT_DEBUG_START
    s = "Generate points for " + getName();
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
#ifdef TIME_STATS
    double pollStartTime = NOMAD::Clock::getCPUTime();
    double pollEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS

    // 1- Generate points for all poll methods
    generateTrialPoints();

    // 2- Evaluate points
    if ( ! _stopReasons->checkTerminate() )
    {
        evalTrialPoints(this);
        pollSuccessful = (getSuccessType() >= NOMAD::SuccessType::FULL_SUCCESS);
    }

    // 3- Future developments: Ortho N+1 would be here.

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
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    // Compute hMax and update Barrier.
    postProcessing();
}


// Generate new points to evaluate
// Note: Used whether MEGA_SEARCH_POLL is true or false.
void NOMAD::Poll::generateTrialPoints()
{
    for (auto pollMethod : _pollMethods)
    {
        if (_stopReasons->checkTerminate())
        {
            break;
        }

        OUTPUT_INFO_START
        AddOutputInfo("Generate points for " + pollMethod->getName(), true, false);
        OUTPUT_INFO_END

        pollMethod->generateTrialPoints();

        // Pass the points from Poll method to Poll for later evaluation
        auto pollMethodPoints = pollMethod->getTrialPoints();
        for (auto point : pollMethodPoints)
        {
            insertTrialPoint(point);
        }
        OUTPUT_INFO_START
        AddOutputInfo("Generated " + std::to_string(pollMethodPoints.size()) + " points");
        AddOutputInfo("Generate points for " + _name, false, true);
        OUTPUT_INFO_END

    }

    // Stopping criterion
    if (0 == getTrialPointsCount())
    {
        auto madsStopReasons = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get(_stopReasons);
        madsStopReasons->set(NOMAD::MadsStopType::MESH_PREC_REACHED);
    }
}


void NOMAD::Poll::computePrimarySecondaryPollCenters(
                        std::vector<NOMAD::EvalPoint> &primaryCenters,
                        std::vector<NOMAD::EvalPoint> &secondaryCenters) const
{
    auto barrier = getMegaIterationBarrier();
    if (nullptr != barrier)
    {
        auto firstXFeas = barrier->getFirstXFeas();
        auto firstXInf  = barrier->getFirstXInf();
        bool primaryIsInf = false;
        NOMAD::Double rho = _runParams->getAttributeValue<NOMAD::Double>("RHO");
        // Negative rho means make no distinction between primary and secondary polls.
        bool usePrimarySecondary = (rho >= 0) && (nullptr != firstXFeas) && (nullptr != firstXInf);
        if (usePrimarySecondary)
        {
            auto evc = NOMAD::EvcInterface::getEvaluatorControl();
            auto evalType = evc->getEvalType();
            auto computeType = evc->getComputeType();
            NOMAD::Double fFeas = firstXFeas->getF(evalType, computeType);
            NOMAD::Double fInf  = firstXInf->getF(evalType, computeType);
            if (fFeas.isDefined() && fInf.isDefined()
                && (fFeas - rho) > fInf)
            {
                // xFeas' f is too large, use xInf as primary poll instead.
                primaryIsInf = true;
            }
        }

        if (usePrimarySecondary)
        {
            if (primaryIsInf)
            {
                primaryCenters   = barrier->getAllXInf();
                secondaryCenters = barrier->getAllXFeas();
            }
            else
            {
                primaryCenters   = barrier->getAllXFeas();
                secondaryCenters = barrier->getAllXInf();
            }
        }
        else
        {
            // All points are primary centers.
            primaryCenters = barrier->getAllPoints();
        }
    }
}


std::shared_ptr<NOMAD::PollMethodBase> NOMAD::Poll::createPollMethod(const bool isPrimary, const NOMAD::EvalPoint& frameCenter) const
{
    std::shared_ptr<NOMAD::PollMethodBase> pollMethod;

    // Select the poll methods to be executed
    NOMAD::DirectionType dirType;
    if (isPrimary)
    {
        dirType = _runParams->getAttributeValue<DirectionType>("DIRECTION_TYPE");
    }
    else
    {
        dirType = _runParams->getAttributeValue<DirectionType>("DIRECTION_TYPE_SECONDARY_POLL");
    }

    switch (dirType)
    {
        case DirectionType::ORTHO_2N:
            pollMethod = std::make_shared<NOMAD::Ortho2NPollMethod>(this, frameCenter);
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
        default:
            throw NOMAD::Exception(__FILE__, __LINE__,"Poll method" + directionTypeToString(dirType) + " is not available.");
            break;
    }

    return pollMethod;
}
