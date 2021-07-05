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
 \file   RNG.hpp
 \brief  Custom class for random number generator
 \author Christophe Tribes and Sebastien Le Digabel
 \date   2011-09-28
 \see    RNG.cpp
 */

#ifndef __NOMAD_4_0_RNG__
#define __NOMAD_4_0_RNG__

#include "../nomad_platform.hpp"
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
    typedef uint32_t result_type;

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return UINT32_MAX; }

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
     * The set seed works like a reset. The private seed used by RNG is always reset.
     \param s   The seed -- \b IN.
     */
    static void setSeed(int s);

    /// Get a random integer
    /**
     \return    An integer in the interval [0,UINT32_MAX].
     */
    static result_type rand();

    /// Functor to get a random integer
    /**
     \return    An integer in the interval [0,UINT32_MAX].
     */
    result_type operator()() { return rand(); }

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
    
    /// Get private values. Used by PyNomad to obtain RNG state.
    static string getPrivateSeedAsString()
    {
        return std::to_string(_x)+ " " + std::to_string(_y)+ " " + std::to_string(_z);
    }

    /// Reset seed to given values.
    /**
     Used by PyNomad to set back RNG in a known state.
     */
    static void setPrivateSeed(uint32_t x, uint32_t y, uint32_t z)
    {
        _x = x;
        _y = y;
        _z = z;
    }
    
    /// Reset seed to given values provided as a single string.
    /**
     Used to set back RNG in a known state.
     */
    static void setPrivateSeedAsString(const std::string & private_seeds)
    {
        istringstream ss(private_seeds);
        uint32_t ps;
        ss >> ps;
        if (ps <= UINT32_MAX )
            _x = ps;
        ss >> ps;
        if (ps <= UINT32_MAX )
            _y = ps;
        ss >> ps;
        if (ps <= UINT32_MAX )
            _z = ps;

    }
    

private:

    DLL_UTIL_API static uint32_t x_def, y_def, z_def;    ///< Initial values for the random number generator
    DLL_UTIL_API static uint32_t _x, _y, _z;             ///< Current values for the random number generator

    DLL_UTIL_API static int _s;


};

#include "../nomad_nsend.hpp"


#endif // __NOMAD_4_0_RNG__
