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

#ifndef __NOMAD_4_3_ITERATIONUTILS__
#define __NOMAD_4_3_ITERATIONUTILS__

#include <stdexcept>


#include "../Algos/Iteration.hpp"
#include "../Algos/TrialPointStats.hpp"
#include "../Algos/MegaIteration.hpp"
#include "../Algos/Step.hpp"
#include "../Eval/SuccessStats.hpp"

#ifdef USE_IBEX
#include "__build__/ibex.h"
#endif

#include "../nomad_nsbegin.hpp"


/// Class of utils (attributes and helper functions) for some phases of an algorithm that involve Iteration.
/**
    The class is in charge of the trial points produced during an algorithm iteration or an algorithm without interation.
    The derived classes for different algorithms provide more additional utils (NMIterationUtils).
    The trial points should be inserted, snapped (to bounds and mesh) and then evaluated. An exception is triggered if points are added after evaluation.
 */
class IterationUtils
{
protected:

    EvalPointSet _trialPoints; ///< The points generated during the start(). Used for run() and postProcessing().

    size_t _nbEvalPointsThatNeedEval;

    const Step * _parentStep;

    /**
     Iteration ancestor (may be _parentStep).
     If _fromAlgo is true, this is NULL.
     */
    Iteration* _iterAncestor;

    /**
     Mega Iteration ancestor (may be _parentStep).
     If _fromAlgo is true, this is NULL.
     */
    MegaIteration* _megaIterAncestor;

    /**
         Stats for trial points generated and evaluated.
     */
    TrialPointStats _trialPointStats;
    
    #ifdef USE_IBEX
    NOMAD::Point projectWithIbex(NOMAD::Point point);
    #endif
    

    
    
private:
    /**
     Flag: True if the direct parent step is an Algorithm. False otherwise. \n
     Used when evaluating trial points without mesh and frame center.
     */
    bool _fromAlgo;

    SuccessType _trialPointsSuccess; ///< Success type of trial points evaluation.
    
public:
    /// Constructor
    /**
     \param parentStep      The calling iteration Step -- \b IN.
     */
    explicit IterationUtils(const Step* parentStep)
      : _trialPoints(),
        _nbEvalPointsThatNeedEval(0),
        _parentStep(parentStep),
        _trialPointsSuccess(SuccessType::UNDEFINED),
        _iterAncestor(nullptr),
        _trialPointStats( parentStep ),
        _fromAlgo(false)
    {
        init();
    }

    /// Destructor
    /**
     Clear trial points
     */
    virtual ~IterationUtils()
    {
        _trialPoints.clear();
    }

    /*---------*/
    /* Get/Set */
    /*---------*/
    const SuccessType& getTrialPointsSuccessType() const       { return _trialPointsSuccess; }

    size_t getTrialPointsCount() const              { return _trialPoints.size(); }
    const EvalPointSet& getTrialPoints() const      { return _trialPoints; }
    
    
    

    /*---------------*/
    /* Other methods */
    /*---------------*/
    /// Insert a trial point
    bool insertTrialPoint(const EvalPoint &evalPoint);

    /// Clear trial points
    void clearTrialPoints( void )
    {
        _trialPoints.clear();
    }

    /// Helper for end task
    /**
     Post-processing of the points after evaluation.
     For instance, computation of a new hMax and update of the MegaIteration's Barrier.
     \return \c true if some value changed (ex. Barrier, hMax), \c false if nothing happened.
     */
    virtual bool postProcessing();

    /// Helper for start task
    /**
     Verify that all points in trialPoints are on the current mesh.
     If a point is not on the mesh the function returns false.
     All points are snapped to the bounds and the mesh. Sometimes, this verification fails anyway. Depending on the algorithm calling we may trigger an exception or display a warning.
     */
    bool verifyPointsAreOnMesh(const std::string& name) const;

    /// Snap a given trial point to the bounds and project on mesh
    /**
     * Used by classes that generate points: SearchMethods, Poll, etc,
     * to make the point satisfy the bounds and mesh requisites,
     * before sending it to evaluation.
     \param evalPoint    The point to process
     \param lowerBound   The lower bounds.
     \param upperBound   The upper bounds
     \return             \c true if the function worked, the point is now on mesh and inside bounds
     */
    bool snapPointToBoundsAndProjectOnMesh(EvalPoint& evalPoint,
                                           const ArrayOfDouble& lowerBound,
                                           const ArrayOfDouble& upperBound);

    /// Start evaluation of the trial points
    /**
     * Called by run.
     \param step    Current step.
     \param keepN   Number of points to keep (by default keep all)
     \param removeStepType Remove trial points of this type after evaluation (default UNDEFINED, do not remove).
     \return true if a success was found, false otherwise.
     */
    bool evalTrialPoints(const Step* step,
                         const size_t keepN = INF_SIZE_T,
                         StepType removeStepType = StepType::UNDEFINED);

    /// Get the number of evaluation points in the queue for evaluation
    size_t getNbEvalPointsThatNeededEval() const { return _nbEvalPointsThatNeedEval; }

    /// Generate the trial points of an algorithm iteration before evaluation.
    /** Algorithm iteration steps must implement generateTrialPointImp
     */
    void generateTrialPoints();
    
    /// Add evaluation information from static surrogate or model to trial points
    /*
     This is used for sorting. Sorting can be used to reduce the number of trial points (Poll methods) and also to prioritize the evaluations (Poll and Search methods).
     */
    virtual void completeTrialPointsInformation();

    /// Generate the second pass trial points of an algorithm iteration before evaluation.
    void generateTrialPointsSecondPass();
    
    /// Implementation for first pass to generate the trial points (pure virtual because Search Methods and Poll Methods must implement it)
    virtual void generateTrialPointsImp() = 0;
    
    /// Implementation for second pass to generate the trial points (currently, only Poll Methods must implement it, the default does nothing)
    virtual void generateTrialPointsSecondPassImp() {} ;
    
    void updateStats(TrialPointStats & trialPointStats);
    
protected:
    bool meshIsFinest() const;
    
private:
       
    /// Helper for constructor
    void init();
    
    /// Helper for evalTrialPoints
    void updateStepSuccessStats(const Step* step);
    
    
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_3_ITERATIONUTILS__
