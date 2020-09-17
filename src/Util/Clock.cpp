/**
 \file   Clock.cpp
 \brief  Clock class (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-02
 \see    Clock.hpp
 */
#include "Clock.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
time_t NOMAD::Clock::_real_t0;
clock_t NOMAD::Clock::_CPU_t0 = clock();
const double NOMAD::Clock::_D_CLOCKS_PER_SEC = static_cast<double>(CLOCKS_PER_SEC);

/*---------------------*/
/*   Reset the clock   */
/*---------------------*/
void NOMAD::Clock::reset()
{
    time(&_real_t0);
    _CPU_t0 = clock();
}


/*---------------------------------------------------------*/
/*  compute the wall-clock time (real time) elapsed since  */
/*  the construction of the Clock object                   */
/*---------------------------------------------------------*/
size_t NOMAD::Clock::getRealTime()
{
    time_t t2;
    time(&t2);
    return static_cast<size_t>(difftime(t2, _real_t0));
}
