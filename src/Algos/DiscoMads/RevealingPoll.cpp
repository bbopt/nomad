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
/**
 \file   RevealingPoll.cpp
 \brief  The DiscoMads algorithm poll step: implementation
 \author Solene Kojtych
 \see    RevealingPoll.hpp
 */

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/Mads/DoublePollMethod.hpp"
#include "../../Algos/Mads/NP1UniPollMethod.hpp"
#include "../../Algos/Mads/Ortho2NPollMethod.hpp"
#include "../../Algos/Mads/OrthoNPlus1PollMethod.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Algos/Mads/SinglePollMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"

#include "../../Algos/DiscoMads/RevealingPoll.hpp"

#ifdef TIME_STATS
#include "../../Util/Clock.hpp"
#endif // TIME_STATS

void NOMAD::RevealingPoll::init()
{
    setStepType(NOMAD::StepType::REVEALING_POLL);
    verifyParentNotNull();

    // Correct Poll attribute (from base class)
    _hasSecondPass = false;

    // Set revealing search parameters
    _nbPoints = _runParams->getAttributeValue<size_t>("DISCO_MADS_REVEALING_POLL_NB_POINTS"); 
    _searchRadius = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_REVEALING_POLL_RADIUS"); 
}




// Generate new points to evaluate
void NOMAD::RevealingPoll::generateTrialPointsImp()
{
    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"RevealingPoll must have an iteration ancestor");
    }

    auto barrier = getMegaIterationBarrier();

    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"RevealingPoll needs a barrier");
    }

    //Warning about groups of variables.
    auto varGroups = _pbParams->getAttributeValue<NOMAD::ListOfVariableGroup>("VARIABLE_GROUP");
    if (!varGroups.empty())
    {
        OUTPUT_INFO_START
        AddOutputInfo("VARIABLE_GROUP parameter is disabled during RevealingPoll to ensure density properties.");
        OUTPUT_INFO_END
    }

    // 1. Create random directions and manage group of variables 
    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    std::list<NOMAD::Direction> directionsFullSpace;
    generateDirections(directionsFullSpace,n);


    // 2. Choice of center for revealing search
    NOMAD::EvalPointPtr xFeas = barrier->getCurrentIncumbentFeas();
    NOMAD::EvalPointPtr xInf = barrier->getCurrentIncumbentInf();
    NOMAD::EvalPointPtr centerRevealingPoll;

    if (nullptr != xFeas)
        {
            centerRevealingPoll=xFeas;   
        }
    else if (nullptr != xInf)
        {
            centerRevealingPoll=xInf;
        }
    else
        {
            OUTPUT_INFO_START
            AddOutputInfo("DiscoMads: no feasible nor infeasible incumbent to use as center of revealing search."); 
            OUTPUT_INFO_END        
        }
    
    OUTPUT_DEBUG_START
    AddOutputDebug("Revealing poll center: " + centerRevealingPoll->display());
    OUTPUT_DEBUG_END

    // 3. Create trial points   
    for (auto& it : directionsFullSpace)
    {
        NOMAD::Point pt(n);
        //trial pt = revealing search center + direction
        for (size_t i = 0 ; i < n ; ++i )
        {
            pt[i] = (*centerRevealingPoll)[i] + it[i];
        }

        auto evalPoint = NOMAD::EvalPoint(pt);

        // Insert the point in eval list
        evalPoint.setPointFrom(std::make_shared<NOMAD::EvalPoint>(*centerRevealingPoll), NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        evalPoint.addGenStep(getStepType());
        insertTrialPoint(evalPoint);
    }


    // Snap the points to bounds and mesh
    const auto& searchMethodPoints = getTrialPoints();

    std::list<NOMAD::EvalPoint> snappedTrialPoints;

    for (auto evalPoint : searchMethodPoints)
    {
        ArrayOfDouble lb = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
        ArrayOfDouble ub = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
        if (snapPointToBoundsAndProjectOnMesh(evalPoint, lb, ub))
        {
            snappedTrialPoints.push_back(evalPoint);
            OUTPUT_INFO_START
            std::string s = "Snap point " + evalPoint.display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
    }

    // Re-insert snapped trial points
    clearTrialPoints();
    for (const auto& evalPoint : snappedTrialPoints)
    {
        insertTrialPoint(evalPoint);
    }
    
    OUTPUT_INFO_START
    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + getName());
    OUTPUT_INFO_END
}


// Generate revealing search directions
void NOMAD::RevealingPoll::generateDirections(std::list<NOMAD::Direction> &directions, const size_t n) const
{
    // Get nbPoints different directions in a sphere of radius _searchRadius
    while (directions.size()<_nbPoints)
    {
        directions.clear(); 
        for (size_t i=0; i<_nbPoints; i++)
        {
            NOMAD::Direction newDirection(n, 0.0);
            // random direction in unit sphere
            NOMAD::Direction::computeDirInUnitSphere(newDirection);

            // random direction in sphere of revealing search radius
            for (size_t j=0 ; j <n ; j++ )
                newDirection[j] *= _searchRadius;
              
            directions.push_back(newDirection);
        }

        // delete non unique directions
        directions.unique();
        if(directions.size()<_nbPoints)
        {
            OUTPUT_INFO_START 
            AddOutputInfo("Revealing search method: new computation of random directions to ensure unicity."); 
            OUTPUT_INFO_END   
        }
    }

    // Output
    OUTPUT_INFO_START
    AddOutputInfo("Direction(s) of revealing search method: ");
    for (auto& direction : directions)
    {
        AddOutputInfo(direction.display());
        
    }
    OUTPUT_INFO_END
}


void NOMAD::RevealingPoll::endImp()
{
    // Sanity check. The endImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    // Compute hMax and update Barrier. 
    // Only update hmax and incumbents if full success (because in this case Revealing Poll is the last step of iteration)
    _updateIncumbentsAndHMax = (_trialPointsSuccess >= NOMAD::SuccessType::FULL_SUCCESS);
    postProcessing();
}


