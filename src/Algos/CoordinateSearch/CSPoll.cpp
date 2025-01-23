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
#include "../../Algos/CoordinateSearch/CSPoll.hpp"
#include "../../Algos/CoordinateSearch/CSPollMethod.hpp"

void NOMAD::CSPoll::init()
{
    setStepType(NOMAD::StepType::CS_POLL);
}

void NOMAD::CSPoll::startImp()
{
    // May need some refactoring. This code is similar to Poll::startImp().
    
    // Sanity check.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    // Reset the current counters. The total counters are not reset (done only once when constructor is called).
    _trialPointStats.resetCurrentStats();
    
    
    // Generate CS poll methods
    NOMAD::DirectionTypeList dirTypes = _runParams->getAttributeValue<DirectionTypeList>("DIRECTION_TYPE");
    
    if (dirTypes.size() != 1 || dirTypes[0] != DirectionType::CS)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"CS Poll method only support DirectionType::CS. " + directionTypeToString(dirTypes[0]) + " not supported for CS.");
    }
    
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

void NOMAD::CSPoll::createPollMethods(const bool isPrimary, const EvalPointPtr& frameCenter)
{
    _frameCenters.push_back(frameCenter);
    auto pollMethod = std::make_shared<NOMAD::CSPollMethod>(this, frameCenter);
    _pollMethods.push_back(pollMethod);
}


void NOMAD::CSPoll::setMeshPrecisionStopType()
{
    auto csStopReasons = NOMAD::AlgoStopReasons<NOMAD::CSStopType>::get(_stopReasons);
    csStopReasons->set(NOMAD::CSStopType::MESH_PREC_REACHED);
}
