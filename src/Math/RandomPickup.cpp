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


