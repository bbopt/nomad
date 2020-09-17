
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelMegaIteration.hpp"
#include "../../Algos/QuadModel/QuadModelInitialization.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"

#include "../../../ext/sgtelib/src/Surrogate_Factory.hpp"
//

void NOMAD::QuadModelAlgo::init()
{
    setName("QuadModel");
    verifyParentNotNull();

    // Instanciate quad model initialization class
    _initialization = std::make_unique<NOMAD::QuadModelInitialization>(this);

}


/*-------------------------*/
/*       Destructor        */
/*-------------------------*/
NOMAD::QuadModelAlgo::~QuadModelAlgo()
{
}


// Start is executed when QuadModelAlgo is used as an algorithm on its own.
void NOMAD::QuadModelAlgo::startImp()
{
    // Default algorithm start. Manages initialization among other things.
    NOMAD::Algorithm::startImp();

    // Comment to appear at the end of stats lines
    setAlgoComment("(QuadModelAlgo)");

}


bool NOMAD::QuadModelAlgo::runImp()
{
    bool success = false;

    size_t k = 0;   // Iteration number

    if (!_termination->terminate(k))
    {
        // Barrier constructor automatically finds the best points in the cache.
        // Barrier is used for MegaIteration management.

        auto barrier = _initialization->getBarrier();
        if (nullptr == barrier)
        {
            auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
            barrier = std::make_shared<NOMAD::Barrier>(hMax, NOMAD::SubproblemManager::getSubFixedVariable(this), NOMAD::EvalType::BB);
        }

        NOMAD::SuccessType megaIterSuccessType = NOMAD::SuccessType::NOT_EVALUATED;

        // TODO fix this
        /*
        if (nullptr != _megaIteration)
        {
            // Case hot restart
            k       = _megaIteration->getK();
            barrier = _megaIteration->getBarrier();
            megaIterSuccessType = _megaIteration->getSuccessType();
        }
        */


        // A single megaiteration is done

        // Create an MegaIteration: manage multiple iterations around
        // different frame centers at the same time.
        NOMAD::QuadModelMegaIteration megaIteration(this, k, barrier, megaIterSuccessType);
        megaIteration.start();
        bool currentMegaIterSuccess = megaIteration.run();
        megaIteration.end();

        success = success || currentMegaIterSuccess;

        // Remember these values to construct the next MegaIteration.
        k       = megaIteration.getK();
        barrier = megaIteration.getBarrier();
        megaIterSuccessType = megaIteration.NOMAD::MegaIteration::getSuccessType();

        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }

        // member _megaIteration is used for hot restart (read and write)
        // Update it here.
        _megaIteration = std::make_shared<NOMAD::QuadModelMegaIteration>(this, k++, barrier, megaIterSuccessType);

    }

    _termination->start();
    _termination->run();
    _termination->end();

    NOMAD::OutputQueue::Flush();

    return success;
}


void NOMAD::QuadModelAlgo::endImp()
{
    // Remove any remaining points from eval queue.
    EvcInterface::getEvaluatorControl()->clearQueue();

    resetPreviousAlgoComment();
    NOMAD::Algorithm::endImp();
}
