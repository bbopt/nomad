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
 \file   RandomPickup.hpp
 \brief  Class for randomly picking up integers
 \author Sebastien Le Digabel
 \date   2010-04-07
 \see    RandomPickup.cpp
 */

#ifndef __NOMAD_4_0_RANDOM_PICKUP__
#define __NOMAD_4_0_RANDOM_PICKUP__

#include <cstdlib>
#include "../Util/Uncopyable.hpp"

using namespace std;

#include "../nomad_nsbegin.hpp"


/// Class for randomly picking up integers.
/**
   - The integers are chosen in [0;n-1] and are distinct.
   - Example displaying 5 different integers in [0;4]:
   \code
   NOMAD::RandomPickup rp(5);
   for (size_t i = 0 ; i < 5 ; ++i)
   {
       std::cout << rp.pickup() << std::endl;
   }
   \endcode
*/
class RandomPickup : private Uncopyable {

private:
    size_t _n0;     ///< Initial value of \c n.
    size_t _n;      ///< Current value of \c n.
    size_t *_elems; ///< Elements that have not been chosen yet.

public:
    /// Constructor.
    /**
       \param n -- The unsigned integer \c n defining the range
                   of values that can be picked up -- \b IN.
    */
    explicit RandomPickup(const size_t n);

    /// Destructor.
    virtual ~RandomPickup() { delete [] _elems; }

    /// Get number of remaining values
    size_t getN() const { return _n; }

    /// Reset.
    void reset();

    /// Randomly pick up an element in [0;n-1].
    /**
       \return The element.
    */
    size_t pickup();

};

#include "../nomad_nsend.hpp"


#endif // __NOMAD_4_0_RANDOM_PICKUP__
