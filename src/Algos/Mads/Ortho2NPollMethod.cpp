
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/Ortho2NPollMethod.hpp"

void NOMAD::Ortho2NPollMethod::init()
{
    setName("Ortho 2N Poll Method");
    verifyParentNotNull();

}

// Generate poll directions
void NOMAD::Ortho2NPollMethod::generateUnitPollDirections(std::list<NOMAD::Direction> &directions, size_t n) const
{
    directions.clear();

    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);

    // Ortho MADS 2n
    // Householder Matrix
    // A vintage piece of code
    NOMAD::Direction** H = new NOMAD::Direction*[2*n];

    // Ordering D_k alternates Hk and -Hk instead of [H_k -H_k]
    for (size_t i = 0; i < n; ++i)
    {
        directions.push_back(NOMAD::Direction(n, 0.0));
        H[i]   = &(directions.back());
        directions.push_back(NOMAD::Direction(n, 0.0));
        H[i+n] = &(directions.back());
    }
    // Householder transformations on the 2n directions on a unit n-sphere
    NOMAD::Direction::householder(dirUnit, true, H);
    delete [] H;

}
