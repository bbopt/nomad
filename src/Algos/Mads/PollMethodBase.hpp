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
#ifndef __NOMAD_4_0_POLLMETHODBASE__
#define __NOMAD_4_0_POLLMETHODBASE__

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Step.hpp"
#include "../../Math/Direction.hpp"

#include "../../nomad_nsbegin.hpp"


/// Class for generic poll method of MADS. Run by Poll.
/**
 */
class PollMethodBase: public Step, public IterationUtils
{
private:
    const EvalPoint _frameCenter;
    const bool _hasNPlus1;      // Second pass for Ortho N+1 methods

public:
    /// Constructor
    /**
     /param parentStep      The parent of this poll step -- \b IN.
     */
    explicit PollMethodBase(const Step* parentStep,
                            const EvalPoint& frameCenter,
                            const bool hasNPlus1 = false)
      : Step(parentStep),
        IterationUtils(parentStep),
        _frameCenter(frameCenter),
        _hasNPlus1(hasNPlus1)
    {
        init();
    }

    bool hasNPlus1() const { return _hasNPlus1; }
    
    
    /// Implementation of startImp.
    /**
      Do nothing.
      Point generation is done in Poll.
     */
    void startImp() override {}

    /// Implementation of endImp.
    /**
      Do nothing.
      Evaluation is done in Poll.
     */
    bool runImp() override { return true; }

    /// Implementation of endImp
    /**
      Do nothing.
      postProcessing is done in Poll.
    */
    void endImp() override {}

    /// Intermediate function used by generateTrialPoints
    std::list<NOMAD::Direction> generateFullSpaceScaledDirections(bool isNPlus1, std::shared_ptr<NOMAD::MeshBase> mesh = nullptr);
    
    /// Intermediate function (not yet the implementation that generates the trial points)
    /**
     - Display before and after generation comments.
     - Launch the implementation of the poll method to generate the trial points (::generateTrialPointsImp).
     - Snap the points to bounds and mesh.
     */
    void generateTrialPoints() override;
    
    /// Intermediate function to compute N+1th point, to be used when _hasNPlus1 is true.
    /*
       \param inputTrialPoints  Trial points generated from the "N" part, evaluated, and
                                not successful.
    */
    virtual void generateTrialPointsNPlus1(const NOMAD::EvalPointSet& inputTrialPoints);

    /// Generate poll directions on a unitary frame. See derived classes (Ortho2nPollMethod, Np1UniPollMethod,...) for implementations.
    virtual void generateUnitPollDirections(std::list<Direction> &directions, const size_t dim) const = 0;

protected:
    void init();
    
    /// Compute 2n directions (from which n directions will be chosen).
    /// Used in Ortho 2N and in Ortho N+1.
    /**
     \param directions  The 2n directions obtained for this poll -- \b OUT.
     \param n           The dimension of the variable space  -- \b IN.
      */
    void generate2NDirections(std::list<NOMAD::Direction> &directions, size_t n) const;
    
    /// Generate n+1th direction. Optionally reimplemented (in Ortho N+1).
    virtual void generateNPlus1Direction(std::list<Direction> &directions) const {};

private:
    
    /// Private method to handle general case and also N+1 case
    void generateTrialPointsInternal(const bool isNPlus1 = false);
    
    /// Scale and project on mesh poll directions.
    /**
     /param dirs      The unit directions to be scaled and projected on mesh -- \b IN/OUT.
     */
    void scaleAndProjectOnMesh(std::list<Direction> & dirs, std::shared_ptr<NOMAD::MeshBase> mesh = nullptr );


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_POLLMETHODBASE__

