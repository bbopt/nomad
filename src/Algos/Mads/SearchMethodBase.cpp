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

#include "../../Algos/Mads/SearchMethodBase.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::SearchMethodBase::init()
{
    // A search method must have a parent
    verifyParentNotNull();

}


void NOMAD::SearchMethodBase::endImp()
{
    // Compute hMax and update Barrier.
    postProcessing();

    // Need to reimplement end() to set a stop reason for Mads based on the search method stop reason
}


void NOMAD::SearchMethodBase::generateTrialPoints()
{

    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + getName(), true, false);
    OUTPUT_INFO_END

    generateTrialPointsImp();

    // Snap the points to bounds and mesh
    auto searchMethodPoints = getTrialPoints();
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");

    std::list<NOMAD::EvalPoint> snappedTrialPoints;
    for (auto evalPoint : searchMethodPoints)
    {
        if (snapPointToBoundsAndProjectOnMesh(evalPoint, lowerBound, upperBound))
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
    for (auto evalPoint : snappedTrialPoints)
    {
        insertTrialPoint(evalPoint);
    }

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + getName(), false, true);
    OUTPUT_INFO_END
}
