/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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

/**
 \file   SimpleRNG.cpp
 \brief  Custom class for random number generator (implementation)
 \author Christophe Tribes and Sebastien Le Digabel
 \date   2011-09-28
 \see    SimpleRNG.hpp
 */

#include "SimpleRNG.hpp"
#include <cmath>

#ifdef _MSC_VER
#include <io.h>
#include <process.h>
#else
#include <unistd.h>
#endif


void NOMAD::SimpleRNG::setSeed(int s)
{
    if (s <= 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "NOMAD::SimpleRNG::setSeed: invalid seed. Seed should be in (0,INT_MAX]");
    }
    _x_def = s;
    
    // Reset the random number generator state according to the seed
    reset();
}

uint32_t NOMAD::SimpleRNG::rand()
{
    // http://madrabbit.org/~ray/code/xorshf96.c //period 2^96-1
    uint32_t t;
    _x ^= _x << 16;
    _x ^= _x >> 5;
    _x ^= _x << 1;
    
    t = _x;
    _x = _y;
    _y = _z;
    _z = t ^ _x ^ _y;
    
    return _z;
}

/*----------------------------------------*/
/*          normal random generators       */
/*----------------------------------------*/
double NOMAD::SimpleRNG::normalRand(double mean, double var)
{
    // Box-Muller transformation~\cite{BoMu58}

    double x1, w;

    do
    {
        x1 = rand(-1.0, 1.0);
        double x2 = rand(-1.0, 1.0);
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);

    return pow(var, 0.5) * x1 * w + mean;
}

double NOMAD::SimpleRNG::normalRandMean0(double Var, int Nsample)
{
    double sum = 0.0;
    double a = pow(3.0 * Var, 0.5);
    for (int i = 0; i < Nsample; i++)
    {
        sum += rand(-a, a);
    }
    return sum / pow(Nsample, 0.5);
}
