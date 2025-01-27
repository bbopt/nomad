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
/**
 \file   RevealingPoll.hpp
 \brief  The DiscoMads algorithm poll step
 \author Solene Kojtych
 \see    RevealingPoll.cpp
 */
#ifndef __NOMAD_4_3_REVEALING_POLLFROMPOLL__
#define __NOMAD_4_3_REVEALING_POLLFROMPOLL__

#include <set>

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Mads/SearchMethodBase.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../nomad_nsbegin.hpp"

/// Class for RevealingPoll of DiscoMads algorithm
/**
 * The "revealing poll" in DiscoMads theory is implemented as inherited from Poll as it is required for the convergence analysis
 Generate the trial points (RevealingPoll::startImp), launch evaluation (RevealingPoll::runImp) and postprocessing (RevealingPoll::endImp).
 */
class RevealingPoll: public Poll
{
private:
#ifdef TIME_STATS
    DLL_ALGO_API static double  _pollTime;      ///< Total time spent running the poll
    DLL_ALGO_API static double  _pollEvalTime;  ///< Total time spent evaluating poll points
#endif // TIME_STATS

    size_t _nbPoints;                 ///< nb of points to generate during revealing poll
    NOMAD::Double _searchRadius;      ///< radius of the revealing poll

public:
    /// Constructor
    /**
     \param parentStep The parent of this poll step
     */
    explicit RevealingPoll(const Step* parentStep)
      : Poll(parentStep)
    {
        init();
    }
    virtual ~RevealingPoll() {}

#ifdef TIME_STATS
    /// Time stats
    static double getPollTime()       { return _pollTime; }
    static double getPollEvalTime()   { return _pollEvalTime; }
#endif // TIME_STATS


private:

    /// Helper for constructor
    void init() ;

    /// Implementation for end tasks for revealing poll.
    /**
     Call the IterationUtils::postProcessing of the points with flag to not update hmax and incumbents, except if full success
     */
    virtual void  endImp() override ;
    
    /// Generate new points to evaluate
    /**
     Implementation called by IterationUtils::generateTrialPoints.
     The trial points are obtained by:
        - adding poll directions (Poll::setPollDirections) to the poll center (frame center).
        - snapping points (and directions) to bounds.
        - projecting points on mesh.
     */
    void generateTrialPointsImp() override;


    ///Helper for generateTrialPointsImp
    // Generate random poll directions for the revealing search
    /**
     \param directions  The direction computed -- \b OUT.
     \param n           The dimension of the variable space -- \b IN.
      */
     void generateDirections(std::list<Direction> &directions, const size_t n) const;
    

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_3_REVEALING_POLLFROMPOLL__
