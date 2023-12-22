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

#include <algorithm>    // For std::merge and std::unique

#include "../../nomad_platform.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/NelderMead/NMInitializeSimplex.hpp"
#include "../../Algos/NelderMead/NMIteration.hpp"
#include "../../Algos/NelderMead/NMUpdate.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"
#include "../../Algos/NelderMead/NMReflective.hpp"
#include "../../Algos/NelderMead/NMShrink.hpp"

void NOMAD::NMIteration::init()
{
    setStepType(NOMAD::StepType::NM_ITERATION);
    
    if (nullptr != _runParams)
    {
        _nmOpt = _runParams->getAttributeValue<bool>("NM_OPTIMIZATION");
        _nmSearchStopOnSuccess = _runParams->getAttributeValue<bool>("NM_SEARCH_STOP_ON_SUCCESS");
    }
    
}

void NOMAD::NMIteration::startImp()
{
    
    // Create the initial simplex Y if it is empty. Use a center pt and the cache
    NOMAD::NMInitializeSimplex initSimplex ( this );
    initSimplex.start();
    bool initSuccess = initSimplex.run();
    initSimplex.end();
    
    if ( ! initSuccess )
    {
        auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( getAllStopReasons() );
        nmStopReason->set( NOMAD::NMStopType::INITIAL_FAILED );
        
        _success = NOMAD::SuccessType::NO_TRIALS;
    }
}


bool NOMAD::NMIteration::runImp()
{
    // Sequential run of NM steps among INITIAL, ( REFLECT, EXPANSION, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION ), SHRINK

    // NMIteration cannot generate all points before evaluation
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    bool iterationSuccess = false;

    // Use a single NMReflect object to perform REFLECT, EXPAND, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION of this iteration
    NOMAD::NMReflective reflect( this );

    // Start with the REFLECT step
    NOMAD::StepType stepType = NOMAD::StepType::NM_REFLECT;



    // Running an NM iteration consists in performing
    // 1) A Reflect
    // 2) An Expansion or an Inside contraction or an Outside contraction or Continue to next iteration (no shrink).
    // 3) Possibly a Shrink.
    while (  ! _stopReasons->checkTerminate() && stepType != NOMAD::StepType::NM_CONTINUE && stepType != NOMAD::StepType::NM_SHRINK )
    {
        // Need to set the current step type before starting
        reflect.setCurrentNMStepType( stepType );

        // Create trial points and evaluate them
        reflect.start();
        reflect.run();
        reflect.end();

        // The NM step type for the next pass
        stepType = reflect.getNextNMStepType() ;

        // Update the type of success for passing to the MegaIteration
        NOMAD::SuccessType reflect_success = reflect.getSuccessType();

        if ( reflect_success > _success )
        {
            // NM Search can be stopped on success
            if ( reflect_success == NOMAD::SuccessType::FULL_SUCCESS && !_nmOpt && _nmSearchStopOnSuccess )
            {
                auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( _stopReasons );
                nmStopReason->set( NOMAD::NMStopType::NM_STOP_ON_SUCCESS );
            }
            iterationSuccess = true; // At least a partial success is obtained
            _success = reflect_success;
        }
    }

    // Perform SHRINK only for a standalone NM optimization
    // Shrink because Reflect step has failed (last chance) or when requested
    if ( (_stopReasons->checkTerminate() || stepType == NOMAD::StepType::NM_SHRINK)  &&
         _nmOpt )
    {
        NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get(_stopReasons)->setStarted();
        
        // Create shrink trial points and evaluate them
        NMShrink shrink ( this );
        shrink.start();
        shrink.run();
        shrink.end();

        // The success type is updated by Step::defaultEnd
        // Iteration success if at least a partial success is obtained
        if ( _success >= NOMAD::SuccessType::PARTIAL_SUCCESS )
        {
            iterationSuccess = true;
        }
    }
    // Perform SHRINK only for a standalone NM optimization ELSE stop NM
    if ( ! _stopReasons->checkTerminate() &&
         stepType == NOMAD::StepType::NM_SHRINK && !_nmOpt )
    {
        auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( _stopReasons );
        nmStopReason->set( NOMAD::NMStopType::NM_STOP_NO_SHRINK );

    }
    
    // End of the iteration: iterationSuccess is true if we have a success.
    return iterationSuccess;

}
