/**
 \file   RandomPickup.cpp
 \brief  Class for randomly picking up integers
 \author Sebastien Le Digabel
 \date   2010-04-07
 \see    RandomPickup.hpp
 */

#include "../Math/RandomPickup.hpp"
#include "../Math/RNG.hpp"

/*---------------------------------------------------------*/
/*                         constructor                     */
/*---------------------------------------------------------*/
NOMAD::RandomPickup::RandomPickup(const size_t n)
  : _n0(n),
    _n(n),
    _elems(new size_t[n])
{
    for (size_t i = 0; i < n; ++i)
    {
        _elems[i] = i;
    }
}


/*---------------------------------------------------------*/
/*                           reset                         */
/*---------------------------------------------------------*/
void NOMAD::RandomPickup::reset()
{
    _n = _n0;
    for (size_t i = 0 ; i < _n ; ++i)
    {
        _elems[i] = i;
    }
}


/*---------------------------------------------------------*/
/*                randomly pick up an element              */
/*---------------------------------------------------------*/
size_t NOMAD::RandomPickup::pickup()
{
    if (_n == 0)
    {
        return 0;
    }
    size_t ind = NOMAD::RNG::rand() % _n;
    size_t tmp = _elems[ind];
    if (ind < _n-1)
    {
        _elems[ind ] = _elems[_n-1];
        _elems[_n-1] = tmp;
    }
    --_n;

    return tmp;
}


