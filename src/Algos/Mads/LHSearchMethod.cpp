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

#include "../../Algos/Mads/LHSearchMethod.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Math/LHS.hpp"
#include "../../Type/LHSearchType.hpp"

void NOMAD::LHSearchMethod::init()
{
    _name = "Latin Hypercube Search Method";
    setComment("(LH)");

    auto lhSearch = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
    setEnabled(lhSearch.isEnabled());
}


void NOMAD::LHSearchMethod::generateTrialPoints()
{
    AddOutputInfo("Generate points for " + _name, true, false);
    auto lhSearch = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
    const NOMAD::MadsIteration* iteration = dynamic_cast<const NOMAD::MadsIteration*>(getParentOfType<NOMAD::MadsIteration*>());
    auto mesh = iteration->getMesh();
    auto iterNumber = iteration->getK();

    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    size_t p = (0 == iterNumber) ? lhSearch.getNbInitial() : lhSearch.getNbIteration();
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    int seed = _runParams->getAttributeValue<int>("SEED");
    auto frameCenter = iteration->getFrameCenter();
    if (nullptr == frameCenter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"LHSearchMethod: frame center is NULL");
    }

    // Update undefined values of lower and upper bounds to use values based
    // on DeltaFrameSize.
    // Based on the code in NOMAD 3, but slightly different.
    // If we used INF values instead of these, we get huge values for the
    // generated points. It is not elegant.
    NOMAD::ArrayOfDouble deltaFrameSize = iteration->getMesh()->getDeltaFrameSize();
    NOMAD::Double scaleFactor = sqrt(-log(NOMAD::DEFAULT_EPSILON));

    for (size_t i = 0; i < n; i++)
    {
        if (!lowerBound[i].isDefined())
        {
            lowerBound[i] = (*frameCenter)[i] - 10.0 * deltaFrameSize[i] * scaleFactor;
        }
        if (!upperBound[i].isDefined())
        {
            upperBound[i] = (*frameCenter)[i] + 10.0 * deltaFrameSize[i] * scaleFactor;
        }
    }

    // Apply Latin Hypercube algorithm
    NOMAD::LHS lhs(n, p, lowerBound, upperBound, seed);
    auto pointVector = lhs.Sample();

    for (auto point : pointVector)
    {
        // Make an EvalPoint from the Point.
        // We do not need the Eval part of EvalPoint right now,
        // but it will be used soon. Could be refactored, but
        // not high priority. Note that an EvalPointSet compares
        // the Point part of the EvalPoints only.

        // Projection without scale 
        if (snapPointToBoundsAndProjectOnMesh(point, lowerBound, upperBound, frameCenter, mesh))
        {
            NOMAD::EvalPoint evalPoint(point);

            // Test if the point could be inserted correctly
            bool inserted = insertTrialPoint(evalPoint);
            std::string s = "Generated point";
            s += (inserted) ? ": " : " not inserted: ";
            s += evalPoint.display();
            AddOutputInfo(s);
        }
    }

    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);

}
