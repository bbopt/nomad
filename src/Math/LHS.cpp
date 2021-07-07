/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/

#include "../Math/LHS.hpp"
#include "../Math/RNG.hpp"
#include "../Math/RandomPickup.hpp"
#include "../Util/Exception.hpp"

#include <algorithm>    // for shuffle


// Constructor
NOMAD::LHS::LHS(const size_t n,
                const size_t p,
                const NOMAD::ArrayOfDouble& lowerBound,
                const NOMAD::ArrayOfDouble& upperBound,
                const NOMAD::Point& frameCenter,
                const NOMAD::ArrayOfDouble& deltaFrameSize,
                const NOMAD::Double& scaleFactor)
:   _n(n),
    _p(p),
    _lowerBound(lowerBound),
    _upperBound(upperBound)
{
    // Update undefined values of lower and upper bounds to use values based
    // on deltaFrameSize.
    // Based on the code in NOMAD 3, but slightly different.
    // Do not use INF values for bounds, that will generate points with huge
    // values. It is not elegant.
    if (frameCenter.isComplete() && deltaFrameSize.isComplete() && scaleFactor.isDefined())
    {
        for (size_t i = 0; i < n; i++)
        {
            if (!_lowerBound[i].isDefined())
            {
                _lowerBound[i] = frameCenter[i] - 10.0 * deltaFrameSize[i] * scaleFactor;
            }
            if (!_upperBound[i].isDefined())
            {
                _upperBound[i] = frameCenter[i] + 10.0 * deltaFrameSize[i] * scaleFactor;
            }
        }
    }

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
        NOMAD::Point point(_n);
        for (size_t i = 0; i < _n; i++)
        {
            NOMAD::Double r_ij = RNG::rand(0,1);
            NOMAD::Double l_i = _lowerBound[i];
            NOMAD::Double Pi_ij( (double)Pi[i][j] );
            NOMAD::Double pdouble( (double)_p );
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
