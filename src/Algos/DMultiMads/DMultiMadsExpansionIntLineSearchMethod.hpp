/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
#ifndef __NOMAD_4_5_DMULTIMADSEXPANSIONINTLINESEARCHMETHOD__
#define __NOMAD_4_5_DMULTIMADSEXPANSIONINTLINESEARCHMETHOD__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/SearchMethodSimple.hpp"
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"

#include "../../nomad_nsbegin.hpp"

/// Expansion integer linesearch method for DMultiMads
/**
 The expansion integer linesearch method is a backtracking linesearch method.
 It tries to expand the current set of solutions by exploring along a direction
 for integer variables only.
 */
class DMultiMadsExpansionIntLineSearchMethod final : public SearchMethodSimple
{
private:
    std::shared_ptr<NOMAD::BarrierBase> _ref_dmads_barrier;

    NOMAD::ArrayOfDouble _lb;
    NOMAD::ArrayOfDouble _ub;
    std::vector<NOMAD::BBInputType> _bbInputTypes;

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit DMultiMadsExpansionIntLineSearchMethod(const NOMAD::Step* parentStep)
      : SearchMethodSimple(parentStep)
    {
        init();
    }

private:

    /// Helper for constructor.
    /**
     Test if the DmultiMads expansion integer linesearch is enabled or not.
     */
    void init();

    void preRunValidations();

    NOMAD::Direction computePrimitiveDirection(const NOMAD::Point& frameCenter,
                                               const NOMAD::Point& pointFrom,
                                               int& initStepSize) const;

    int computeMaxStepSize(const NOMAD::Point& frameCenter,
                           const NOMAD::Direction& dir,
                           const NOMAD::ArrayOfDouble& lb,
                           const NOMAD::ArrayOfDouble& ub) const;

    bool isInBarrier(const NOMAD::Point& x) const;

    /// Generate new points (no evaluation)
    /**
     \copydoc SearchMethodAlgo::generateTrialPointsFinal 
     The expansion integer linesearch method generates points along a direction
     for integer variables only.
     */
     void generateTrialPointsFinal() override;
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_DMULTIMADSEXPANSIONINTLINESEARCHMETHOD__
