
#include "../../Algos/Mads/SearchMethodSimple.hpp"

bool NOMAD::SearchMethodSimple::runImp()
{

    bool foundBetter = false;
    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }

    return foundBetter;
}


void NOMAD::SearchMethodSimple::startImp()
{
    if ( ! _stopReasons->checkTerminate() )
    {
        // Create EvalPoints and snap to bounds and snap on mesh
        generateTrialPoints();
    }
}
