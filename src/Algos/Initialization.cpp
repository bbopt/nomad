#include "../Algos/Initialization.hpp"
#include "../Output/OutputQueue.hpp"

NOMAD::Initialization::~Initialization()
{
    NOMAD::OutputQueue::Flush();
}


void NOMAD::Initialization::init()
{
    _name = getAlgoName() + "Initialization";
    verifyParentNotNull();
}
