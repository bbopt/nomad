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

#include "../../Algos/Mads/SearchMethod.hpp"

void NOMAD::SearchMethod::init()
{
    // A search method must have a parent
    verifyParentNotNull();
}


void NOMAD::SearchMethod::startImp()
{
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);


    if ( ! _stopReasons->checkTerminate() )
    {
        // Create EvalPoints
        generateTrialPoints();
        // NB verifyPointsAreOnMesh and updatePointsWithFrameCenter are
        // done here, so that if an user adds a new search method, they
        // will not have to think about adding these verifications.
        verifyPointsAreOnMesh(getName());
        updatePointsWithFrameCenter();
    }
}


bool NOMAD::SearchMethod::runImp()
{
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        // Show more information in the form of an AlgoComment.
        NOMAD::MainStep::setAlgoComment(getComment());

        foundBetter = evalTrialPoints(this);
        NOMAD::MainStep::resetPreviousAlgoComment();

    }
    return foundBetter;

}


void NOMAD::SearchMethod::endImp()
{
    // Compute hMax and update Barrier.
    postProcessing(getEvalType());

    // Need to reimplement end() to set a stop reason for Mads based on the search method stop reason
}
