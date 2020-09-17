/**
 \file   Clock.hpp
 \brief  Clock class (headers)
 \author Sebastien Le Digabel
 \date   2010-04-02
 \see    Clock.cpp
 */
#ifndef __NOMAD400_CLOCK__
#define __NOMAD400_CLOCK__

#include <ctime>

#include "../nomad_nsbegin.hpp"

#ifdef _MSC_VER
#pragma warning(disable:4275)
#pragma warning(disable:4251)
#ifdef DLL_EXPORTS
#define DLL_API __declspec(dllexport)
#else
#define DLL_API __declspec(dllimport)
#endif
#else
#define DLL_API
#endif


/// Clock class.
/**
 Time measurement.\n\n
 \b Example:
 \code
 std::cout << "elapsed real time = " << Clock::getRealTime() << std::endl;
 std::cout << "elapsed CPU time  = " << Clock::getCPUTime()  << std::endl;
 \endcode
 */
class DLL_API Clock {

private:

    static time_t       _real_t0;          ///< Wall clock time measurement.
    static clock_t      _CPU_t0;           ///< CPU time measurement.
    static const double _D_CLOCKS_PER_SEC; ///< System constant for CPU time measurement.

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


#endif // __NOMAD400_CLOCK__
