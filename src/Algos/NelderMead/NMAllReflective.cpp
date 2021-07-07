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

#include "../../Algos/NelderMead/NMAllReflective.hpp"
#include "../../Algos/NelderMead/NMReflective.hpp"
#include "../../Algos/SubproblemManager.hpp"


void NOMAD::NMAllReflective::startImp()
{
    if ( ! _stopReasons->checkTerminate() )
    {
        // The iteration start function manages the simplex creation.
        NMIteration::startImp();

        // Generate REFLECT, EXPANSION, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION (no SHRINK)
        // All points are generated before evaluation
        verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, true);

        generateTrialPoints();
        verifyPointsAreOnMesh(getName());
    }
}


void NOMAD::NMAllReflective::generateTrialPoints ()
{
    NOMAD::NMReflective reflect( this );

    // Need to set the current step type before starting
    reflect.setCurrentNMStepType( NOMAD::StepType::NM_REFLECT );

    // Create trial points but no evaluation
    reflect.start();
    reflect.end();
    auto trialPts = reflect.getTrialPoints();
    for (auto evalPoint : trialPts)
    {
        evalPoint.addGenStep(getStepType());
        insertTrialPoint(evalPoint);
    }

    // Expand simplex
    if ( ! _stopReasons->checkTerminate() )
    {
        reflect.setCurrentNMStepType( NOMAD::StepType::NM_EXPAND );
        reflect.start();
        reflect.end();
        trialPts = reflect.getTrialPoints();
        for (auto evalPoint : trialPts)
        {
            evalPoint.addGenStep(getStepType());
            insertTrialPoint(evalPoint);
        }

    }

    // Inside contraction of simplex
    if ( ! _stopReasons->checkTerminate() )
    {
        reflect.setCurrentNMStepType( NOMAD::StepType::NM_INSIDE_CONTRACTION );
        reflect.start();
        reflect.end();
        trialPts = reflect.getTrialPoints();
        for (auto evalPoint : trialPts)
        {
            evalPoint.addGenStep(getStepType());
            insertTrialPoint(evalPoint);
        }

    }

    // Outside contraction of simplex
    if ( ! _stopReasons->checkTerminate() )
    {
        reflect.setCurrentNMStepType( NOMAD::StepType::NM_OUTSIDE_CONTRACTION );
        reflect.start();
        reflect.end();
        trialPts = reflect.getTrialPoints();
        for (auto evalPoint : trialPts)
        {
            evalPoint.addGenStep(getStepType());
            insertTrialPoint(evalPoint);
        }

    }

    // If everything is ok we terminate a single NM iteration completed anyway
    if ( ! _stopReasons->checkTerminate() )
    {
        auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( getAllStopReasons() );
        nmStopReason->set(NOMAD::NMStopType::NM_SINGLE_COMPLETED);
    }
}
