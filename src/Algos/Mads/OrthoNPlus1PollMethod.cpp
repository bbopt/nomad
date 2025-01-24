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
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/OrthoNPlus1PollMethod.hpp"
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#ifdef USE_SGTELIB
#include "../../Algos/QuadModel/QuadModelSinglePass.hpp"
#endif
#include "../../Algos/SubproblemManager.hpp"
#include "../../Math/Direction.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::OrthoNPlus1PollMethod::init()
{
    if (_flagUseQuadOpt)
    {
        setStepType(NOMAD::StepType::POLL_METHOD_ORTHO_NPLUS1_QUAD);
    }
    else
    {
        setStepType(NOMAD::StepType::POLL_METHOD_ORTHO_NPLUS1_NEG);
    }
    verifyParentNotNull();
    
    // The second pass directions do not need to be scaled on mesh.
    // Indeed, the n+1 the direction can be strictly inside the frame.
    _scaleAndProjectSecondPassDirectionOnMesh=false;


}

// Generate poll directions
void NOMAD::OrthoNPlus1PollMethod::generateUnitPollDirections(std::list<NOMAD::Direction> &directions, const size_t n) const
{
    directions.clear();

    generate2NDirections(directions, n);
}


// Ref in NOMAD 3: NOMAD::Mads::get_single_dynamic_direction.
// Get a single dynamic direction from incomplete poll by sum of negative (Ortho N+1 NEG) or by quad model optimization (Ortho N+1 QUAD).
void NOMAD::OrthoNPlus1PollMethod::generateSecondPassDirections(std::list<Direction> &directions) const
{
    
    // Gather directions for points which are already generated - which are currently in _trialPoints.
    NOMAD::Direction dirSec;
    
    if (!directions.empty())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "OrthoNPlus1PollMethod: directions is only for output.");
    }
    
    std::vector<NOMAD::Direction> allGeneratingDirs;
    
    // Test if we have enough directions
    if (_trialPoints.begin()->size() != _trialPoints.size())
    {
        OUTPUT_DEBUG_START
        const std::string s = "Insufficient number of trial points for second pass: " + std::to_string(_trialPoints.size());
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
        return;
    }
    for (const auto& trialPoint : _trialPoints)
    {
        
        auto dir = trialPoint.getDirection();
        allGeneratingDirs.push_back(*dir);
        if (nullptr != dir)
        {
            if (0 == dirSec.size())
            {
                dirSec = -(*dir);
            }
            else
            {
                dirSec -= (*dir);
            }
        }
    }
    
    NOMAD::Point prospect_point;
    if ( _flagUseQuadOpt )
        optimizeQuadModel( allGeneratingDirs, dirSec ); // dirSec is updated if optim is a success

    auto norm = dirSec.norm();
    if (dirSec.size() > 0 && norm > 0)
    {
        directions.push_back(dirSec);
    }
}

/// Reduce the number of trial points from 2n to n
void NOMAD::OrthoNPlus1PollMethod::trialPointsReduction()
{
    OUTPUT_DEBUG_START
    const std::string s = "Number of trial points before sort and reduction to form a basis: " + std::to_string(_trialPoints.size());
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
    
    if (_trialPoints.empty())
    {
        return;
    }
    
    // First step: sort the trial points. Put them into a vector.
    NOMAD::EvcInterface evcInterface(this);
    
    // Complete trial points information with MODEL or SURROGATE eval
    completeTrialPointsInformation();
    

    // If the mesh is finest we don't want to order the direction based on some informed rule. Some informed ordering can select always the same directions.
    // Forcing random will ensure that the Ortho n+1 direction will grow asymptotically dense
    // NOTE: With GMesh the mesh is refined less frequently then with XMesh (the Ortho N+1 paper was based on XMesh). So we more frequently force a random ordering. The quad model ordering is beneficial when opportunistic. Using the Frame size to decide is also not pertinent. When the mesh index is available we should use it instead of the mesh size.
    // bool forceRandom = meshIsFinest();
    bool forceRandom = false;
    std::vector<EvalPoint> sortedEvalPoints=evcInterface.getSortedTrialPoints(_trialPoints,forceRandom);
    
    size_t n = _frameCenter->size();
    std::vector<EvalPoint> spanningSortedEvalPoints;
    if (sortedEvalPoints.size() < n)
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("Not enough sorted trial points (because trimming). Add some points that have already been evaluated. Needed for second pass.");
        OUTPUT_DEBUG_END
        spanningSortedEvalPoints = sortedEvalPoints;
        for (const auto & trialPt: _trialPoints)
        {
            if ( find(spanningSortedEvalPoints.begin(),spanningSortedEvalPoints.end(),trialPt) == spanningSortedEvalPoints.end() )
            {
                spanningSortedEvalPoints.push_back(trialPt);
                
            }
            if ( spanningSortedEvalPoints.size() >= n)
                break;
        }
    }
    else
    {
        
        // Second step: loop on the points, keep points that increase the rank, up to n points to make a spanning set of directions.
        size_t currentRank=0, rank=0;
        for (const auto & ev: sortedEvalPoints)
        {
            spanningSortedEvalPoints.push_back(ev);
            rank = NOMAD::EvalPoint::getRank(spanningSortedEvalPoints);
            if ( rank > currentRank && rank <= n)
            {
                currentRank++;
            }
            else
            {
                spanningSortedEvalPoints.pop_back();
            }
            if (rank == n)
            {
                break;
            }
        }
    }
    
    // Final step: put back the Eval points in the set of trial points
    _trialPoints.clear();
    for (const auto & evalPoint : spanningSortedEvalPoints)
    {
        _trialPoints.insert(evalPoint);
    }
    
    OUTPUT_DEBUG_START
    const std::string s = "Number of trial points after reduction to form a basis: " + std::to_string(_trialPoints.size());
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
    NOMAD::OutputQueue::Flush();
    
}

bool NOMAD::OrthoNPlus1PollMethod::optimizeQuadModel(const std::vector<NOMAD::Direction> & allGeneratingDirs, NOMAD::Direction & dirSec) const
{
    
#ifdef USE_SGTELIB
    // Send trial EvalPoints to EvaluatorControl
    NOMAD::EvcInterface evcInterface(this);
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    // Reset the counter
    evc->resetModelEval();
    
    auto fullFixedVar = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this);
    OUTPUT_INFO_START
    std::string s = "Create QuadModelEvaluator with fixed variable = ";
    s += fullFixedVar.display();
    AddOutputInfo(s);
    OUTPUT_INFO_END
    
    // The trial points are generated for the frame center.
    
    auto madsIteration = getParentOfType<Iteration*>();
    auto barrier = madsIteration->getMegaIterationBarrier();
    const auto& computeType = barrier->getFHComputeType();
    
    if (nullptr != _frameCenter
        && _frameCenter->getF(computeType).isDefined()
        && _frameCenter->getF(computeType) < MODEL_MAX_OUTPUT)
    {
        NOMAD::QuadModelSinglePass singlePass(this, _frameCenter, madsIteration->getMesh(), allGeneratingDirs);
        
        // Generate the trial points (but we just need the points coordinates to output the direction)
        singlePass.generateTrialPoints();
        
        // Get a single n+1 th direction from the best feasible or best infeasible
        auto bestFeas = singlePass.getBestFeas();
        auto bestInf = singlePass.getBestInf();
        if(nullptr != bestFeas)
        {
            dirSec = *(bestFeas->getX()) - *(_frameCenter->getX());
            return true;
        }
        else if(nullptr != bestInf)
        {
            dirSec = *(bestInf->getX()) - *(_frameCenter->getX()) ;
            return true;
        }
    }

#endif
        
    return false;
}
