/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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
#include "../Util/Exception.hpp"

#include <random>       // for mt19937
#include <algorithm>    // for shuffle


// Constructor
NOMAD::LHS::LHS(size_t n,
                size_t p,
                NOMAD::ArrayOfDouble lowerBound,
                NOMAD::ArrayOfDouble upperBound,
                int seed)
:   _n(n),
    _p(p),
    _lowerBound(lowerBound),
    _upperBound(upperBound),
    _seed(seed)
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

    std::srand(_seed);
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
    std::vector<size_t> v;
    for (size_t j = 1; j <= p; j++)
    {
        v.push_back(j);
    }

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(v.begin(), v.end(), g);

    return v;
}
