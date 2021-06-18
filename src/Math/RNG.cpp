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
/**
 \file   RNG.cpp
 \brief  Custom class for random number generator (implementation)
 \author Christophe Tribes and Sebastien Le Digabel
 \date   2011-09-28
 \see    rng.hpp
 */

#include "../Math/RNG.hpp"
#include <math.h>

// Default values for the provided number seed
int NOMAD::RNG::_s = 0;

uint32_t NOMAD::RNG::x_def = 123456789;
uint32_t NOMAD::RNG::y_def = 362436069;
uint32_t NOMAD::RNG::z_def = 521288629;
uint32_t NOMAD::RNG::_x = x_def;
uint32_t NOMAD::RNG::_y = y_def;
uint32_t NOMAD::RNG::_z = z_def;


void NOMAD::RNG::setSeed(int s)
{
    if (s <= INT_MAX && s >= 0)
    {
        _s = s;
    }
    else
    {
        throw NOMAD::Exception (__FILE__, __LINE__,
                                "NOMAD::RNG::setSeed(): invalid seed. Seed should be in [0,INT_MAX]");
    }

    // Always reset the private seed to default values
    resetPrivateSeedToDefault();

    for (int i = 0; i < _s; i++)
    {
        NOMAD::RNG::rand();
    }
}


uint32_t NOMAD::RNG::rand()
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
double NOMAD::RNG::normalRand(double mean, double var)
{
    // Box-Muller transformation~\cite{BoMu58}

    double x1, w;

    do
    {
        x1 = NOMAD::RNG::rand(-1.0,1.0);
        double x2 = NOMAD::RNG::rand(-1.0,1.0);
        w  = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );

    return pow(var,0.5) * x1 * w + mean;
}


double NOMAD::RNG::normalRandMean0(double Var, int Nsample)
{
    double sum = 0.0;
    double a = pow(3.0 * Var, 0.5);
    for (int i = 0 ; i < Nsample ; i++)
    {
        sum += NOMAD::RNG::rand(-a,a);
    }
    return sum / pow(Nsample, 0.5);
}
