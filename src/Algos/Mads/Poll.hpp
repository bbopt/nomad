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
#ifndef __NOMAD_4_4_POLL__
#define __NOMAD_4_4_POLL__

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
#ifdef TIME_STATS
    DLL_ALGO_API static double  _pollTime;      ///< Total time spent running the poll
    DLL_ALGO_API static double  _pollEvalTime;  ///< Total time spent evaluating poll points
#endif // TIME_STATS
    
    DirectionTypeList _primaryDirectionTypes, _secondaryDirectionTypes;  ///< The poll methods implement different direction types for primary and secondary poll centers.

    NOMAD::Double _rho; ///< Rho parameter of the progressive barrier. Used to choose if the primary frame center is the feasible or infeasible  incumbent.
    
    size_t _trialPointMaxAddUp; ///< Add new trial points to the ones produced by the selected direction type up to reached a given number.
    
protected:
    bool _hasSecondPass;   ///<  Ortho n+1 poll methods generate n trial points in a first pass and, if not successful, generate the n+1 th point (second pass)
    
    std::vector<std::shared_ptr<PollMethodBase>> _pollMethods;  ///< Unlike for Search, Poll methods generate all their points and only then they are evaluated.

    std::vector<EvalPointPtr> _frameCenters;  ///< The frame centers (primary and secondary) of the poll methods. See createPollMethods. 
    
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



    /// Second pass of point generation after first pass is not successfull.
    /**
      For Ortho N+1 methods the N+1th point is produced after evaluating N sorted points without success.
      In this second pass, the trial points are determined using the evaluations of the first N points.
      */
    void generateTrialPointsSecondPass();
    
    
    /// Extra trial point generation after first pass.
    /**
      Add trial points to reached a prescribed trial point number (TRIAL_POINT_MAX_ADD_UP).
     Only for single pass direction type (ortho 2n).
      */
    void generateTrialPointsExtra();

#ifdef TIME_STATS
    /// Time stats
    static double getPollTime()       { return _pollTime; }
    static double getPollEvalTime()   { return _pollEvalTime; }
#endif // TIME_STATS

protected:
    /// Helper for start: get lists of Primary and Secondary Polls
    void computePrimarySecondaryPollCenters(std::vector<EvalPointPtr> &primaryCenters, std::vector<EvalPointPtr> &secondaryCenters) const;
    
    /// Helper for start: create poll methods
    virtual void createPollMethods(const bool isPrimary, const EvalPointPtr frameCenter);
    
    /// Helper for generateTrialPoints
    ///  Set the stop type for the Algorithm (can be reimplemented, for example CS)
    virtual void setMeshPrecisionStopType();
    
    /// Helper to create poll methods for current poll centers
    virtual void createPollMethodsForPollCenters();
    

private:
    /// Helper for constructor
    void init();


    /// Implementation for start tasks for MADS poll.
    /**
     Call to generate the poll methods
     */
    virtual void    startImp() override ;

    /// Implementation for run tasks for MADS poll.
    /**
     Call poll methods and perform trial points evaluation.
     \return Flag \c true if found better solution \c false otherwise.
     */
    virtual bool    runImp() override;

    /// Implementation for end tasks for MADS poll.
    /**
     Call the IterationUtils::postProcessing of the points.
     */
    virtual void    endImp() override ;

    /// Generate new points to evaluate
    /**
     Implementation called by IterationUtils::generateTrialPoints.
     The trial points are obtained by:
        - adding poll directions (Poll::setPollDirections) to the poll center (frame center).
        - snaping points (and directions) to bounds.
        - projecting points on mesh.
     */
    void generateTrialPointsImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_POLL__
