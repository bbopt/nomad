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
/**
 \file   RNG.hpp
 \brief  Custom class for random number generator
 \author Christophe Tribes and Sebastien Le Digabel 
 \date   2011-09-28
 \see    RNG.cpp
 */

#ifndef __NOMAD400_RNG__
#define __NOMAD400_RNG__

#include "../Util/defines.hpp"
#include "../Util/Exception.hpp"

using namespace std;

#include "../nomad_nsbegin.hpp"


/// Class for random number generator
/**
This class is used to set a seed for the random number generator and get a random integer or a random double between two values. \n
 http://madrabbit.org/~ray/code/xorshf96.c with period 2^96-1
 */
class RNG {
    
public:
     /// Get current seed
    /**
     \return An integer in [0,UINT32_MAX].
     */
    static int getSeed()
    {
        return _s;
    }
    
    /// Set seed
    /**
     \param s   The seed -- \b IN.
     */
    static void setSeed(int s);

    
    /// Get a random integer
    /**
     \return    An integer in the interval [0,UINT32_MAX].
     */
    static uint32_t rand();
    
    
    /// Get a random number having a normal distribution as double
    /**
     \param a   Lower bound  -- \b IN.
     \param b   Upper bound  -- \b IN.
     \return    A double in the interval [a,b].
     */
    static double rand(double a, double b)
    {
        return a+((b-a)*RNG::rand())/UINT32_MAX;
    }
    
    /// Get a random number using a normal distribution centered on 0
    /**
     * Get a random number approaching a normal distribution (N(0,Var)) as double
     *
     *
     \param Var     Variance of the target normal distribution    -- \b IN.
     \param Nsample Number of samples for averaging                -- \b IN.
     \return        A double in the interval [-sqrt(3*Var);+sqrt(3*Var)].
     */
    static double normalRandMean0(double Var = 1, int Nsample = 12);
    
    
    /// Get a random number approaching a normal distribution N(Mean,Var) as double.
    /**
     A series of Nsample random numbers Xi in the interval [-sqrt(3*Var);+sqrt(3*Var)] is used -> E[Xi] = 0, Var(Xi) = var. \n
     See http://en.wikipedia.org/wiki/Central_limit_theorem
     
     \param Mean    Mean of the target normal distribution        -- \b IN.
     \param Var     Variance of the target normal distribution    -- \b IN.
     \return        A random number.
     */
    static double normalRand(double Mean = 0, double Var = 1);
    
    /// Reset seed to its default value
    static void resetPrivateSeedToDefault()
    {
        _x = x_def;
        _y = y_def;
        _z = z_def;
    }

    /// Get private values
    static void getPrivateSeed(uint32_t &x, uint32_t &y, uint32_t &z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    /// Reset seed to given values.
    /**
     Used to set back RNG in a known state.
     */
    static void setPrivateSeed(uint32_t x, uint32_t y, uint32_t z)
    {
        _x = x;
        _y = y;
        _z = z;
    }
    
private:
    
    static uint32_t x_def, y_def, z_def;    ///< Initial values for the random number generator
    static uint32_t _x, _y, _z;             ///< Current values for the random number generator
    
    static int _s;
};

#include "../nomad_nsend.hpp"


#endif // __NOMAD400_RNG__
