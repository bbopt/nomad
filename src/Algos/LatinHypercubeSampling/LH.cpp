/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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

#include "../../Algos/EvcInterface.hpp"

#include "../../Algos/LatinHypercubeSampling/LH.hpp"

#include "../../Math/LHS.hpp"
void NOMAD::LH::init()
{
    _name = "Latin Hypercube Sampling";
    verifyParentNotNull();

}

void NOMAD::LH::startImp()
{

    // Comment to appear at the end of stats lines
    NOMAD::MainStep::setAlgoComment("(LH)");

    generateTrialPoints();

}

void NOMAD::LH::generateTrialPoints()
{

    AddOutputInfo("Generate points for " + _name, true, false);

    auto lhEvals = _runParams->getAttributeValue<size_t>("LH_EVAL");
    if (NOMAD::INF_SIZE_T == lhEvals)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "The number of evaluations for LH cannot be infinite.");
    }

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");

    if (!lowerBound.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,_name + " requires a complete lower bound vector");
    }

    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    if (!upperBound.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,_name + " requires a complete upper bound vector");
    }

    int seed = _runParams->getAttributeValue<int>("SEED");

    // Apply Latin Hypercube algorithm
    NOMAD::LHS lhs(n, lhEvals, lowerBound, upperBound, seed);
    auto pointVector = lhs.Sample();

    for (auto point : pointVector)
    {
        // Make an EvalPoint from the Point.
        // We do not need the Eval part of EvalPoint right now,
        // but it will be used soon. Could be refactored, but
        // not high priority. Note that an EvalPointSet compares
        // the Point part of the EvalPoints only.

        // Projection without scale
        NOMAD::EvalPoint evalPoint(point);

        // Test if the point is inserted correctly
        bool inserted = insertTrialPoint(evalPoint);
        std::string s = "Generated point";
        s += (inserted) ? ": " : " not inserted: ";
        s += evalPoint.display();
        AddOutputInfo(s);
    }

    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);

}

bool NOMAD::LH::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }
    auto LHStopReasons = NOMAD::AlgoStopReasons<NOMAD::LHStopType>::get( _stopReasons );
    if (  _stopReasons->testIf( NOMAD::EvalStopType::ALL_POINTS_EVALUATED ) )
    {
        LHStopReasons->set( NOMAD::LHStopType::ALL_POINTS_EVALUATED );
    }

    return foundBetter;
}

void NOMAD::LH::endImp()
{
    // Remove any remaining points from eval queue.
    EvcInterface::getEvaluatorControl()->clearQueue();

    // reset to the previous stats comment
    NOMAD::MainStep::resetPreviousAlgoComment();

}
