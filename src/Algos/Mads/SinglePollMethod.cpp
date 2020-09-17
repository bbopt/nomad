
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/SinglePollMethod.hpp"

void NOMAD::SinglePollMethod::init()
{
    setName("Single Poll Method");
    verifyParentNotNull();
}

// Generate a poll direction
void NOMAD::SinglePollMethod::generateUnitPollDirections(std::list<NOMAD::Direction> &directions, size_t n) const
{
    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);
    directions.push_back(dirUnit);
}

