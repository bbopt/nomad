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
#ifndef __NOMAD_4_0_NMSEARCHMETHOD__
#define __NOMAD_4_0_NMSEARCHMETHOD__

#include "../../Algos/Mads/SearchMethodAlgo.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to perform a Search method using Nelder Mead simplex algorithm.
/**
 If Nelder Mead search is enabled (check is done in NMSearchMethod::init), the NMSearchMethod::run function manages the execution (start, run, end) of the NM algorithm. \n
 The new trial points can be generated during a single pass of all Nelder Mead reflective steps (REFLECT, EXPAND, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION) ( generateTrialPoint ) or as a NM optimization.
 */
class NMSearchMethod final : public SearchMethodAlgo
{
public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit NMSearchMethod(const Step* parentStep )
      : SearchMethodAlgo(parentStep )
    {
        init();
    }


    /**
     Execute (start, run, end) of the NM algorithm. Returns a \c true flag if the algorithm found better point.
     */
    virtual bool runImp() override ;


private:

    /// Helper for constructor.
    /**
     Test if the NM search is enabled or not. Set the maximum number of trial points.
     */
    void init();

    ///Generate new points (no evaluation)
    /**
     \copydoc SearchMethodAlgo::generateTrialPointsImp \n
     Perform one iteration of all reflective steps (Reflect, Expansion, Inside and Outside Contraction). This is just portion of the NM algorithm without iteration.
     */
    virtual void generateTrialPointsImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_NMSEARCHMETHOD__

