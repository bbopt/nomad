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
#ifndef __NOMAD_4_3_ORTHO_N_PLUS_1_POLLMETHOD__
#define __NOMAD_4_3_ORTHO_N_PLUS_1_POLLMETHOD__

#include "../../Algos/Mads/PollMethodBase.hpp"
#include "../../nomad_nsbegin.hpp"

/// Class to perform Ortho N+1 NEG/MOD Poll.
// Ortho MADS N+1 NEG/MOD:
// 1- Generate 2N points
// 2- Sort points and keep only the first N points not already evaluated and that form a basis.
// 3- Evaluate N points
// 4- If no success found, compute N+1 th direction using NEG or MOD
// 5- Evaluate 1 point
class OrthoNPlus1PollMethod final : public PollMethodBase
{
private:
    bool _flagUseQuadOpt;  ///< Flag to use Quad model optim. for N+1 th direction. Otherwise use sum of negative directions.
    
    const EvalPointPtr _frameCenter;

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit OrthoNPlus1PollMethod(const Step* parentStep,
                                     const EvalPointPtr frameCenter,
                                   bool flagUseQuadOpt)
      : PollMethodBase(parentStep, frameCenter, true), // true: hasSecondPass
        _frameCenter(frameCenter),
        _flagUseQuadOpt(flagUseQuadOpt)
    {
        init();
    }
    
private:

    /// Helper for constructor.
    /**
     Test if the OrthoNPlus1 Poll is enabled or not and set step type.
     */
    void init();

    ///Generate 2N polls direction on a unit N-Sphere (no evaluation)
    /**
    The reduction to N directions is done by function trialPointsReduction.
     \param directions  The 2N directions obtained for this poll -- \b OUT.
     \param n           The dimension of the variable space  -- \b IN.
      */
    void generateUnitPollDirections(std::list<Direction> &directions, const size_t n) const override final;

    
    /// Compute N+1th direction and add it to the vector of directions.
    /**
     \param directions  The N+1 th direction obtained for this poll -- \b OUT.
      */
    void generateSecondPassDirections(std::list<Direction> &directions) const override final;
    
    /// Reduce the number of trial points
    /*
     This is used to obtain the N most promising points from 2N points.
     If the mesh is not finest the points are sorted first.
     */
    void trialPointsReduction() override ;
    
    /// Helper for generate second pass directions (only for ortho n+1 quad)
    /*
    \param allGenerationDirections  The first N directions obtained for this poll -- \b IN.
    \param dirSec  The N+1 th direction: IN, negative sum, OUT, result of quad model opt. or negative sum (if not SUCCESS) -- \b IN / \b OUT.
     /return flag for success of optimization
     */
    bool optimizeQuadModel(const std::vector<Direction> & allGeneratingDirs, Direction & dirSec) const;
    

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_3_ORTHO_N_PLUS_1_POLLMETHOD__
