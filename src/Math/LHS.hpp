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

#ifndef __NOMAD_4_0_LHS__
#define __NOMAD_4_0_LHS__

#include <vector>
#include "../Math/Point.hpp"
using namespace std;

#include "../nomad_nsbegin.hpp"


/// \brief Latin Hypercube Sampling class.
/**
 * Input:
 *  n dimension
 *  p number of desired samples
 *  lowerBound and upperBound of type ArrayOfDouble indicating lower and upper bounds. Must be completely defined.
 *  seed (optional)
 *
 * Output:
 *  p points of dimension n, distributed in the l x u hyper-rectangle of R^n
 */
class LHS
{
private:
    size_t _n;  ///< dimension
    size_t _p;  ///< number of samples
    ArrayOfDouble    _lowerBound; ///< lower bounds
    ArrayOfDouble    _upperBound; ///< upper bounds

public:
    /// Constructor
    /**
     \param n               Dimension -- \b IN.
     \param p               Number of samples -- \b IN.
     \param lowerBound      Lower bounds -- \b IN.
     \param upperBound      Upper bounds -- \b IN.
     */
    explicit LHS(const size_t n,
                 const size_t p,
                 const ArrayOfDouble& lowerBound,
                 const ArrayOfDouble& upperBound,
                 const NOMAD::Point& frameCenter = NOMAD::Point(),
                 const NOMAD::ArrayOfDouble& deltaFrameSize = NOMAD::ArrayOfDouble(),
                 const NOMAD::Double& scaleFactor = NOMAD::Double());

    /// Get lower bound
    /**
     \return Lower bound as \c ArrayOfDouble.
     */
    ArrayOfDouble getLowerBound() const                  { return _lowerBound; }

    /// Set lower bound
    /**
     \param lowerBound  An \c ArrayOfDouble for lower bound -- \b IN.
     */
    void setLowerBound(const ArrayOfDouble lowerBound)   { _lowerBound = lowerBound; }

    /// Get upper bound
    /**
     \return Upper bound as \c ArrayOfDouble.
     */
    ArrayOfDouble    getUpperBound(void) const           { return _upperBound; }

    /// Set upper bound
    /**
     \param upperBound  An \c ArrayOfDouble for upper bound -- \b IN.
     */
    void setUpperBound(const ArrayOfDouble upperBound)   { _upperBound = upperBound; }

    /// Do the sampling
    /**
     \return A vector \c Points.
     */
    std::vector<Point> Sample() const;

    /// Random permutation of the vector (1, 2, .., p)
    /**
     \param p   Number of positive integers elements in series -- \b IN.
     \return    Vector of positive integers
     */
    static std::vector<size_t> Permutation(const size_t p);
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_LHS__
