
#include <algorithm>    // For std::merge and std::unique

#include "../Algos/Iteration.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Type/CallbackType.hpp"


NOMAD::Iteration::~Iteration()
{
    NOMAD::OutputQueue::Flush();
}


void NOMAD::Iteration::init()
{
    _name = getAlgoName() + "Iteration " + std::to_string(_k);
    verifyParentNotNull();
}



void NOMAD::Iteration::endImp()
{
    AddOutputInfo("Stop reason: " + _stopReasons->getStopReasonAsString() );

    if ( _runParams->getAttributeValue<bool>("USER_CALLS_ENABLED") )
    {
        bool stop = false;

        /// Callback user provided function to check if user requested a stop.
        runCallback(NOMAD::CallbackType::ITERATION_END, *this, stop);


        if (!_stopReasons->checkTerminate() && stop)
        {
            _stopReasons->set(NOMAD::BaseStopType::USER_STOPPED);
        }
    }
}




