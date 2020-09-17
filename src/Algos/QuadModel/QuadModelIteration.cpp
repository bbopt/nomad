
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelOptimize.hpp"
#include "../../Algos/QuadModel/QuadModelUpdate.hpp"
#include "../../../ext/sgtelib/src/Surrogate_Factory.hpp"

void NOMAD::QuadModelIteration::reset()
{
    if (nullptr != _model)
    {
        _model.reset();
    }

    if (nullptr != _trainingSet)
    {
        _trainingSet.reset();
    }
}

void NOMAD::QuadModelIteration::init()
{
    _name = getAlgoName() + NOMAD::Iteration::getName();

    // Count the number of constraints
    const auto bbot = NOMAD::QuadModelAlgo::getBBOutputType();
    size_t nbConstraints = NOMAD::getNbConstraints(bbot);

    // Init the TrainingSet
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    SGTELIB::Matrix empty_X("empty_X", 0, static_cast<int>(n));
    SGTELIB::Matrix empty_Z("empty_Z", 0, static_cast<int>(nbConstraints+1));
    _trainingSet = std::make_shared<SGTELIB::TrainingSet>(empty_X, empty_Z);

    // The quadratic model uses Sgtelib
    _model = std::shared_ptr<SGTELIB::Surrogate>(SGTELIB::Surrogate_Factory(*_trainingSet, "TYPE PRS"));

}


void NOMAD::QuadModelIteration::startImp()
{
    incK();

    // Select the sample points to construct the model. Use a center pt and the cache
    NOMAD::QuadModelUpdate update(this);
    update.start();
    bool updateSuccess = update.run();
    update.end();

    if ( ! updateSuccess )
    {
        auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );

        // The initial update is not a success. If the stop reason is not set to terminate we set a default stop reason for initialization.
        if ( !_stopReasons->checkTerminate() )
            qmsStopReason->set( NOMAD::ModelStopType::INITIAL_FAIL);
    }
}


bool NOMAD::QuadModelIteration::runImp()
{

    bool iterationSuccess = false;

    // Initialize optimize member - model optimizer on sgte
    NOMAD::QuadModelOptimize optimize (this, _pbParams);

    // Model Update is handled in start().
    if (!_stopReasons->checkTerminate() && _model->is_ready() )
    {
        // Optimize to find oracle points on this model
        optimize.start();
        iterationSuccess = optimize.run();
        optimize.end();
    }

    // Update MegaIteration success type (use deconstification!)
    NOMAD::SuccessType success = optimize.getSuccessType();
    auto megaIter = getParentOfType<NOMAD::MegaIteration*>();
    megaIter->setSuccessType(success);

    // End of the iteration: iterationSuccess is true if we have a success.
    return iterationSuccess;

}
