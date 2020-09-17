
#include "../Math/LHS.hpp"
#include "../Math/RNG.hpp"
#include "../Math/RandomPickup.hpp"
#include "../Util/Exception.hpp"

#include <algorithm>    // for shuffle


// Constructor
NOMAD::LHS::LHS(size_t n,
                size_t p,
                NOMAD::ArrayOfDouble lowerBound,
                NOMAD::ArrayOfDouble upperBound)
:   _n(n),
    _p(p),
    _lowerBound(lowerBound),
    _upperBound(upperBound)
{
    if (!_lowerBound.isComplete())
    {
        std::string s = "LHS Lower bound needs to be completely defined. Values given: ";
        s += lowerBound.display();
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
    if (!_upperBound.isComplete())
    {
        std::string s = "LHS Upper bound needs to be completely defined. Values given: ";
        s += upperBound.display();
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
}


// Do the sample
// Audet & Hare Algorithm 3.9 Latin Hypercube Sampling
std::vector<NOMAD::Point> NOMAD::LHS::Sample() const
{
    std::vector<NOMAD::Point> samplepoints;

    // 0 - Initialization
    // Let Pi be a n x p matrix in which each of its n rows
    // is a random permutation of the vector (1, 2, .., p).
    //
    std::vector<std::vector<size_t>> Pi;
    for (size_t i = 0; i < _n; i++)
    {
        std::vector<size_t> v = Permutation(_p);
        Pi.push_back(v);
    }

    // 1 - Sample construction
    for (size_t j = 0; j < _p; j++)
    {
        Point point(_n);
        for (size_t i = 0; i < _n; i++)
        {
            NOMAD::Double r_ij = RNG::rand(0,1);
            NOMAD::Double l_i = _lowerBound[i];
            NOMAD::Double Pi_ij( Pi[i][j] );
            NOMAD::Double pdouble( _p );
            NOMAD::Double u_i( _upperBound[i] );

            NOMAD::Double x_ij = l_i + (Pi_ij - r_ij) / pdouble * (u_i - l_i);
            point[i] = x_ij;

        }
        samplepoints.push_back(point);
    }

    return samplepoints;
}


// Input: p
// Output: Random permutation of the vector (1, 2, .., p)
std::vector<size_t> NOMAD::LHS::Permutation(const size_t p)
{
    NOMAD::RandomPickup rp(p);

    std::vector<size_t> v;
    for (size_t j = 0; j < p ; j++)
    {
        v.push_back(rp.pickup()+1);
    }

    return v;
}
