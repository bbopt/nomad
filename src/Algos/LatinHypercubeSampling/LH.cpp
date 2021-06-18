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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/LatinHypercubeSampling/LH.hpp"
#include "../../Algos/Mads/GMesh.hpp"
#include "../../Math/LHS.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::LH::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_LH);
    verifyParentNotNull();

}

NOMAD::ArrayOfPoint NOMAD::LH::suggest ()
{
    generateTrialPoints();
    
    NOMAD::ArrayOfPoint xs;
    for (auto trialPoint : _trialPoints)
    {
        xs.push_back(*trialPoint.getX());
    }
    return xs;
}


void NOMAD::LH::startImp()
{
    generateTrialPoints();
}


void NOMAD::LH::generateTrialPoints()
{
    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + getName(), true, false);
    OUTPUT_INFO_END

    auto lhEvals = _runParams->getAttributeValue<size_t>("LH_EVAL");
    if (NOMAD::INF_SIZE_T == lhEvals)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "The number of evaluations for LH cannot be infinite.");
    }

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");

    if (!lowerBound.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,getName() + " requires a complete lower bound vector");
    }

    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    if (!upperBound.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,getName() + " requires a complete upper bound vector");
    }

    // Apply Latin Hypercube algorithm
    NOMAD::LHS lhs(n, lhEvals, lowerBound, upperBound);
    auto pointVector = lhs.Sample();

    // Create a mesh and project points on this mesh. It will ensure that parameters
    // BB_INPUT_TYPE and GRANULARITY are satisfied.
    auto mesh = std::make_shared<NOMAD::GMesh>(_pbParams);
    mesh->setEnforceSanityChecks(false);
    // Modify mesh so it is the finest possible.
    // Note: GRANULARITY is already adjusted with regards to BB_INPUT_TYPE.
    NOMAD::ArrayOfDouble newMeshSize = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("GRANULARITY");
    auto eps = NOMAD::Double::getEpsilon();
    for (size_t i = 0; i < newMeshSize.size(); i++)
    {
        if (0 == newMeshSize[i])
        {
            // No granularity: Set mesh size to Epsilon.
            newMeshSize[i] = eps;
        }
    }

    mesh->setDeltas(newMeshSize, newMeshSize);
    auto center = NOMAD::Point(n, 0);

    for (auto point : pointVector)
    {
        // First, project on mesh.
        point = mesh->projectOnMesh(point, center);
        // Second, snap to bounds.
        point.snapToBounds(lowerBound, upperBound);

        NOMAD::EvalPoint evalPoint(point);
        // Test if the point is inserted correctly
        bool inserted = insertTrialPoint(evalPoint);
        OUTPUT_INFO_START
        std::string s = "Generated point";
        s += (inserted) ? ": " : " not inserted: ";
        s += evalPoint.display();
        AddOutputInfo(s);
        OUTPUT_INFO_END
    }

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + getName(), false, true);
    OUTPUT_INFO_END

}


bool NOMAD::LH::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }
    auto LHStopReasons = NOMAD::AlgoStopReasons<NOMAD::LHStopType>::get( _stopReasons );
    if (NOMAD::EvcInterface::getEvaluatorControl()->testIf(NOMAD::EvalMainThreadStopType::ALL_POINTS_EVALUATED))
    {
        LHStopReasons->set( NOMAD::LHStopType::ALL_POINTS_EVALUATED );
    }

    return foundBetter;
}


void NOMAD::LH::endImp()
{
}
