
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/NP1UniPollMethod.hpp"
#include "../../Algos/Mads/Ortho2NPollMethod.hpp"

void NOMAD::NP1UniPollMethod::init()
{
    setName("N+1 Uniform Poll Method");
    verifyParentNotNull();
}


void NOMAD::NP1UniPollMethod::generateUnitPollDirections(std::list<Direction> &directions, size_t n) const
{
    directions.clear();

    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);

    // Ortho MADS 2n
    // Householder Matrix
    // A vintage piece of code
    NOMAD::Direction** H = new NOMAD::Direction*[2*n];

    // Ordering D_k alternates Hk and -Hk instead of [H_k -H_k]
    std::list<Direction> vDirs;
    for (size_t i = 0; i < n; ++i)
    {
        vDirs.push_back(NOMAD::Direction(n, 0.0));
        H[i]   = &(vDirs.back());
        vDirs.push_back(NOMAD::Direction(n, 0.0));
        H[i+n] = &(vDirs.back());
    }
    // Householder transformations on the 2n directions on a unit n-sphere
    NOMAD::Direction::householder(dirUnit, true, H);

    // dir 0
    NOMAD::Direction dir0(*H[0]);
    for ( size_t i = 1 ; i < n ; ++i )
    {
        dir0=dir0+(*H[i]);
    }
    dir0*=-1.0/sqrt(double(n));
    directions.push_back(dir0);

    NOMAD::Double beta=(sqrt(double(n+1.0))-1.0)/sqrt(double(n));
    dir0*=beta;

    for ( size_t i = 0 ; i < n ; i++ )
    {
        NOMAD::Direction diri(*H[i]);
        diri*=sqrt(double(n+1));
        diri=diri+dir0;
        diri*=1.0/sqrt(double(n));

        directions.push_back(diri);
    }
    delete [] H;
}
 // end generateTrialPoints
