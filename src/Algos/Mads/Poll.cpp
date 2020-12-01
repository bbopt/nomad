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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/DoublePollMethod.hpp"
#include "../../Algos/Mads/NP1UniPollMethod.hpp"
#include "../../Algos/Mads/Ortho2NPollMethod.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Algos/Mads/SinglePollMethod.hpp"
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

    // Select the single poll method to be executed
    auto dt = _pbParams->getAttributeValue<DirectionType>("DIRECTION_TYPE");
    switch (dt)
    {
        case DirectionType::ORTHO_2N:
            _pollMethod = std::make_shared<NOMAD::Ortho2NPollMethod>(this);
            break;
        case DirectionType::NP1_UNI:
            _pollMethod = std::make_shared<NOMAD::NP1UniPollMethod>(this);
            break;
        case DirectionType::DOUBLE:
            _pollMethod = std::make_shared<NOMAD::DoublePollMethod>(this);
            break;
        case DirectionType::SINGLE:
            _pollMethod = std::make_shared<NOMAD::SinglePollMethod>(this);
            break;
        default:
            throw NOMAD::Exception(__FILE__, __LINE__,"Poll method"+directionTypeToString(dt)+"not available.");
            break;
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

    NOMAD::SuccessType bestSuccessYet = NOMAD::SuccessType::NOT_EVALUATED;
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    // Go through all poll methods until we get a success.
    OUTPUT_DEBUG_START
    s = "Poll method " + _pollMethod->getName() + " is enabled";
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
#ifdef TIME_STATS
    double pollStartTime = NOMAD::Clock::getCPUTime();
    double pollEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS
    _pollMethod->start();
    _pollMethod->run();
    success = _pollMethod->getSuccessType();
    pollSuccessful = (success >= NOMAD::SuccessType::FULL_SUCCESS);
    if (success > bestSuccessYet)
    {
        bestSuccessYet = success;
    }
    _pollMethod->end();
#ifdef TIME_STATS
    _pollTime += NOMAD::Clock::getCPUTime() - pollStartTime;
    _pollEvalTime += NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime() - pollEvalStartTime;
#endif // TIME_STATS


    OUTPUT_INFO_START
    s = _pollMethod->getName();
    s += (pollSuccessful) ? " is successful" : " is not successful";
    s += ". Stop reason: ";
    s += _stopReasons->getStopReasonAsString() ;
    AddOutputInfo(s);
    OUTPUT_INFO_END

    setSuccessType(bestSuccessYet);

    return pollSuccessful;
}


void NOMAD::Poll::endImp()
{
    // Sanity check. The endImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    postProcessing(NOMAD::EvcInterface::getEvaluatorControl()->getEvalType());
}





// Generate new points to evaluate
void NOMAD::Poll::generateTrialPoints()
{
    // Sanity check. The generateTrialPoints function should be called only
    // when trial points are also generated for each search method (MegaSearchPoll).
    // After generation, all trials points evaluated at once.
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, true);

    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + _name, true, false);
    OUTPUT_INFO_END

    _pollMethod->generateTrialPoints();

    // Pass the points from Poll method to Poll for later evaluation
     auto pollMethodPoints = _pollMethod->getTrialPoints();
     for (auto point : pollMethodPoints)
     {
          insertTrialPoint(point);
    }

    // Sanity check
    verifyPointsAreOnMesh(getName());

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);
    OUTPUT_INFO_END
}
