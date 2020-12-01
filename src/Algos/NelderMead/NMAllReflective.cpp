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

#include "../../Algos/NelderMead/NMAllReflective.hpp"
#include "../../Algos/NelderMead/NMReflective.hpp"



void NOMAD::NMAllReflective::startImp()
{
    if ( ! _stopReasons->checkTerminate() )
    {
        // The iteration start function manages the simplex creation.
        NMIteration::startImp();

        // Generate REFLECT, EXPANSION, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION (no SHRINK)
        // All points are generated before evaluation
        verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, true);

        generateTrialPoints();
        verifyPointsAreOnMesh(getName());
        updatePointsWithFrameCenter();

    }
}

void NOMAD::NMAllReflective::generateTrialPoints ()
{
    NOMAD::NMReflective reflect( this );

    // Need to set the current step type before starting
    reflect.setCurrentNMStepType( NMStepType::REFLECT );

    // Create trial points but no evaluation
    reflect.start();
    reflect.end();
    auto trialPts = reflect.getTrialPoints();
    for ( const auto & pt : trialPts )
        insertTrialPoint( pt );

    // Expand simplex
    if ( ! _stopReasons->checkTerminate() )
    {
        reflect.setCurrentNMStepType( NMStepType::EXPAND );
        reflect.start();
        reflect.end();
        trialPts = reflect.getTrialPoints();
        for ( const auto & pt : trialPts )
            insertTrialPoint( pt );

    }

    // Inside contraction of simplex
    if ( ! _stopReasons->checkTerminate() )
    {
        reflect.setCurrentNMStepType( NMStepType::INSIDE_CONTRACTION );
        reflect.start();
        reflect.end();
        trialPts = reflect.getTrialPoints();
        for ( const auto & pt : trialPts )
            insertTrialPoint( pt );

    }

    // Outside contraction of simplex
    if ( ! _stopReasons->checkTerminate() )
    {
        reflect.setCurrentNMStepType( NMStepType::OUTSIDE_CONTRACTION );
        reflect.start();
        reflect.end();
        trialPts = reflect.getTrialPoints();
        for ( const auto & pt : trialPts )
            insertTrialPoint( pt );

    }

    // If everything is ok we terminate a single NM iteration completed anyway
    if ( ! _stopReasons->checkTerminate() )
    {
        auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( getAllStopReasons() );
        nmStopReason->set(NOMAD::NMStopType::NM_SINGLE_COMPLETED);
    }

}
