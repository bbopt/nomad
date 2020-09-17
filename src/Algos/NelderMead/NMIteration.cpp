
#include <algorithm>    // For std::merge and std::unique

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/NelderMead/NMIteration.hpp"
#include "../../Algos/NelderMead/NMUpdate.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"
#include "../../Algos/NelderMead/NMReflective.hpp"
#include "../../Algos/NelderMead/NMShrink.hpp"
#include "../../Algos/NelderMead/NMInitializeSimplex.hpp"

void NOMAD::NMIteration::init()
{
    _name = getAlgoName() + "Iteration";

    _bestSuccess = NOMAD::SuccessType::UNSUCCESSFUL;

}

void NOMAD::NMIteration::startImp()
{
    incK();

    // Update main barrier.
    NOMAD::NMUpdate update( this );
    update.start();
    update.run();
    update.end();

    // Create the initial simplex Y if it is empty. Use a center pt and the cache
    NOMAD::NMInitializeSimplex initSimplex ( this );
    initSimplex.start();
    bool initSuccess = initSimplex.run();
    initSimplex.end();

    if ( ! initSuccess )
    {
        auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( getAllStopReasons() );
        nmStopReason->set( NOMAD::NMStopType::INITIAL_FAILED );
    }
}


bool NOMAD::NMIteration::runImp()
{
    // Sequential run of NM steps among INITIAL, ( REFLECT, EXPANSION, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION ), SHRINK

    // NMIteration cannot generate all points before evaluation
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    bool iterationSuccess = false;

    // Use a single NMReflect object to perform REFLECT, EXPAND, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION of this iteration
    NOMAD::NMReflective reflect( this );

    // Start with the REFLECT step
    NMStepType stepType = NMStepType::REFLECT;

    bool nmOpt = _runParams->getAttributeValue<bool>("NM_OPTIMIZATION");
    bool nmSearchStopOnSuccess = _runParams->getAttributeValue<bool>("NM_SEARCH_STOP_ON_SUCCESS");

    // Running an NM iteration consists in performing
    // 1) A Reflect
    // 2) An Expansion or an Inside contraction or an Outside contraction or Continue to next iteration (no shrink).
    // 3) Possibly a Shrink.
    while (  ! _stopReasons->checkTerminate() && stepType != NMStepType::CONTINUE && stepType != NMStepType::SHRINK )
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
        NOMAD::SuccessType success = reflect.getSuccessType();

        if ( success > _bestSuccess )
        {
            // NM Search can be stopped on success
            if ( success == NOMAD::SuccessType::FULL_SUCCESS && !nmOpt && nmSearchStopOnSuccess )
            {
                auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( _stopReasons );
                nmStopReason->set( NOMAD::NMStopType::NM_STOP_ON_SUCCESS );
            }
            iterationSuccess = true; // At least a partial success is obtained
            _bestSuccess = success;
        }
    }

    // Perform SHRINK only for a standalone NM optimization
    if ( ! _stopReasons->checkTerminate() &&
         stepType == NMStepType::SHRINK  &&
         nmOpt )
    {
        // Create shrink trial points and evaluate them
        NMShrink shrink ( this );
        shrink.start();
        shrink.run();
        shrink.end();

        // Update the type of success for passing to the MegaIteration
        NOMAD::SuccessType success = shrink.getSuccessType();
        if ( success > _bestSuccess )
        {
            iterationSuccess = true; // At least a partial success is obtained
            _bestSuccess = success;
        }
    }

    if ( iterationSuccess )
    {
        // Update MegaIteration success type with best success found.
        getParentOfType<NOMAD::MegaIteration*>()->setSuccessType(_bestSuccess);
    }

    // End of the iteration: iterationSuccess is true if we have a success.
    return iterationSuccess;

}
