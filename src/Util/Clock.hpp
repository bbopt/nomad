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
 \file   Clock.hpp
 \brief  Clock class (headers)
 \author Sebastien Le Digabel
 \date   2010-04-02
 \see    Clock.cpp
 */
#ifndef __NOMAD_4_0_CLOCK__
#define __NOMAD_4_0_CLOCK__

#include <ctime>

#include "../nomad_platform.hpp"
#include "../nomad_nsbegin.hpp"

/// Clock class.
/**
 Time measurement.\n\n
 \b Example:
 \code
 std::cout << "elapsed real time = " << Clock::getRealTime() << std::endl;
 std::cout << "elapsed CPU time  = " << Clock::getCPUTime()  << std::endl;
 \endcode
 */
class DLL_UTIL_API Clock {

private:

    static time_t       _real_t0;           ///< Wall clock time measurement.
    static clock_t      _CPU_t0;            ///< CPU time measurement.
    static const double _D_CLOCKS_PER_SEC;  ///< System constant for CPU time measurement.

public:
    // No need for constructor. All is static.

    /// Reset the clock.
    static void reset();

    /// Get wall clock time.
    /**
     \return The time elapsed since _real_t0
     */
    static size_t getRealTime();

    /// Get the CPU time.
    /**
     \return The CPU time elapsed since _CPU_t0
     */
    static double getCPUTime()
    {
        return ( clock() - _CPU_t0 ) / _D_CLOCKS_PER_SEC;
    }

    /// Get time since start or reset.
    /**
     \return The time elapsed since _real_t0
     */
    static size_t getTimeSinceStart() { return getRealTime(); }
};

#include "../nomad_nsend.hpp"


#endif // __NOMAD_4_0_CLOCK__
