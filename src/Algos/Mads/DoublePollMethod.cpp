
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/DoublePollMethod.hpp"

void NOMAD::DoublePollMethod::init()
{
    setName("Double Poll Method");
    verifyParentNotNull();
}

// Generate poll directions
void NOMAD::DoublePollMethod::generateUnitPollDirections(std::list<NOMAD::Direction> &directions, size_t n) const
{
    directions.clear();

    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);
    directions.push_back(dirUnit);

    // insert the opposite direction
    dirUnit *=-1.0;
    directions.push_back(dirUnit);
}
