
#include "../../Algos/QuadModel/QuadModelSinglePass.hpp"
#include "../../Algos/QuadModel/QuadModelOptimize.hpp"
#include "../../Algos/QuadModel/QuadModelUpdate.hpp"
#include "../../Cache/CacheBase.hpp"


void NOMAD::QuadModelSinglePass::generateTrialPoints ()
{
    // Select the sample points to construct the model. Use a center pt and the cache
    NOMAD::QuadModelUpdate update(this);
    update.start();
    bool updateSuccess = update.run();
    update.end();

    // Model Update is handled in start().
    if (!_stopReasons->checkTerminate() && updateSuccess && getModel()->is_ready() )
    {
        // Clear sgte value info from cache. For each pass we suppose we have a different quatric model and Sgte value must be re-evaluated.
        NOMAD::CacheBase::getInstance()->clearSgte();

        // Optimize to generate oracle points on this model
        // Initialize optimize member - model optimizer on sgte
        NOMAD::QuadModelOptimize optimize (this, _pbParams);

        optimize.start();
        // No run, the trial points are evaluated somewhere else.
        optimize.end();

        auto trialPts = optimize.getTrialPoints();
        for ( const auto & pt : trialPts )
            insertTrialPoint( pt );

    }

    // If everything is ok we set the stop reason.
    if (! _stopReasons->checkTerminate())
    {
        auto stopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );
        stopReason->set(NOMAD::ModelStopType::MODEL_SINGLE_PASS_COMPLETED);
    }

}
