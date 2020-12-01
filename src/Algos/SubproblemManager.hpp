/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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

#ifndef __NOMAD400_SUBPROBLEMMANAGER__
#define __NOMAD400_SUBPROBLEMMANAGER__

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include "../Algos/Algorithm.hpp"
#include "../Algos/Subproblem.hpp"

#include "../nomad_nsbegin.hpp"


/// Class to associate Algorithms with Subproblems
/**
 Ease the passage between sub-dimension and full dimension. Algorithm works in
 a sub dimentsion and does not know the full dimension.
 */
class SubproblemManager
{
private:
    static std::map<const Algorithm*, const Subproblem> _map;
#ifdef _OPENMP
    static omp_lock_t _mapLock;
#endif // _OPENMP

public:
    /// Constructor
    explicit SubproblemManager()
    {
        init();
    }

    /// Destructor
    virtual ~SubproblemManager()
    {
        destroy();
    }

    static const Subproblem& getSubproblem(const Step* step);

    static const Point& getSubFixedVariable(const Step* step);

    static void addSubproblem(const Algorithm* algo, const Subproblem& subproblem);

    static void removeSubproblem(const Algorithm* algo);

    static void reset();

private:
    /// Helper for constructor
    void init();

    /// Helper for destructor
    void destroy();
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_SUBPROBLEMMANAGER__
