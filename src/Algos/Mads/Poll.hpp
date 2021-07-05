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
#ifndef __NOMAD_4_0_POLL__
#define __NOMAD_4_0_POLL__

#include <set>

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Mads/PollMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the poll (step 3) of MADS algorithm.
/**
 Generate the trial points (Poll::startImp), launch evaluation (Poll::runImp) and postprocecssing (Poll::endImp).
 */
class Poll: public Step, public IterationUtils
{
private:
    std::vector<std::shared_ptr<PollMethodBase>> _pollMethods;  ///< Unlike for Search, Poll methods generate all their points and only then they are evaluated.
#ifdef TIME_STATS
    DLL_ALGO_API static double  _pollTime;      ///< Total time spent running the poll
    DLL_ALGO_API static double  _pollEvalTime;  ///< Total time spent evaluating poll points
#endif // TIME_STATS


public:
    /// Constructor
    /**
     \param parentStep The parent of this poll step
     */
    explicit Poll(const Step* parentStep)
      : Step(parentStep),
        IterationUtils(parentStep),
        _pollMethods()
    {
        init();
    }
    virtual ~Poll() {}

    /// Generate new points to evaluate
    /**
     The trial points are obtained by:
        - adding poll directions (Poll::setPollDirections) to the poll center (frame center).
        - snaping points (and directions) to bounds.
        - projecting points on mesh.
     */
    void generateTrialPoints() override;

    /// Generate N+1th point for ORTHO N+1 methods
    /**
      The trial point is obtained from the evaluations of the first N points.
      */
    void generateTrialPointsNPlus1(const NOMAD::EvalPointSet& inputTrialPoints);

#ifdef TIME_STATS
    /// Time stats
    static double getPollTime()       { return _pollTime; }
    static double getPollEvalTime()   { return _pollEvalTime; }
#endif // TIME_STATS



private:
    /// Helper for constructor
    void init();

    /// Implementation for start tasks for MADS poll.
    /**
     Call to generate trial points and test for mesh precision
     */
    virtual void    startImp() override ;

    /// Implementation for run tasks for MADS poll.
    /**
     Start trial points evaluation.
     \return Flag \c true if found better solution \c false otherwise.
     */
    virtual bool    runImp() override;

    /// Implementation for end tasks for MADS poll.
    /**
     Call the IterationUtils::postProcessing of the points.
     */
    virtual void    endImp() override ;

    /// Helper for start: get lists of Primary and Secondary Polls
    void computePrimarySecondaryPollCenters(std::vector<EvalPoint> &primaryCenters, std::vector<EvalPoint> &secondaryCenters) const;

    /// Helper for start: create poll methods
    std::vector<std::shared_ptr<PollMethodBase>> createPollMethods(const bool isPrimary, const EvalPoint& frameCenter) const;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_POLL__
