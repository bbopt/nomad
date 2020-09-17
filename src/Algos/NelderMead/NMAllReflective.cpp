
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
