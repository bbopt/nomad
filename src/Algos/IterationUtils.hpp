/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#ifndef __NOMAD400_ITERATIONUTILS__
#define __NOMAD400_ITERATIONUTILS__

#include <stdexcept>

#include "../Algos/Iteration.hpp"
#include "../Algos/Algorithm.hpp"
#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

/// Class of utils (attributes and helper functions) for some phases of an algorithm that involve Iteration start(), run() and end().
/**
    The class is in charge of the trial points.
    The derived classes for different algorithms (MadsIterationUtils, NMIterationUtils) provide more utils.
 */
class IterationUtils
{
private:

    EvalPointSet _trialPoints; ///< The points generated during the start(). Used for run() and postProcessing().

    size_t _nbEvalPointsThatNeedEval;

protected:

    const Step * _parentStep;

    SuccessType _success; ///< Success type of this step.

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


private:
    /**
     Flag: True if the direct parent step is an Algorithm. False otherwise. \n
     Used when evaluating trial points without mesh and frame center.
     */
    bool _fromAlgo;


public:
    /// Constructor
    /**
     \param parentStep      The calling iteration Step -- \b IN.
     */
    explicit IterationUtils(const Step* parentStep)
      : _trialPoints(),
        _nbEvalPointsThatNeedEval(0),
        _parentStep(parentStep),
        _success(SuccessType::NOT_EVALUATED),
        _iterAncestor(nullptr),
        _fromAlgo(false)
    {
        init();
    }

    /// Destructor
    /**
     Clear trial points
     */
    virtual ~IterationUtils() {
        _trialPoints.clear();

    }

    /*---------*/
    /* Get/Set */
    /*---------*/
    const SuccessType& getSuccessType() const       { return _success; }
    void setSuccessType(const SuccessType& success) { _success = success; }

    size_t getTrialPointsCount() const              { return _trialPoints.size(); }
    const EvalPointSet & getTrialPoints() const     { return _trialPoints; }

    /*---------------*/
    /* Other methods */
    /*---------------*/
    /// Insert a trial point
    bool insertTrialPoint(const EvalPoint &evalPoint);

    /// Clear trial points
    void clearTrialPoints( void ) { _trialPoints.clear() ; }

    /// Helper for end task
    /**
     Post-processing of the points after evaluation.
     For instance, computation of a new hMax and update of the MegaIteration's Barrier.
     Use evaluations of the type given by input parameter evalType.
     */
    virtual void postProcessing(const EvalType& evalType);

    /// Helper for start task
    /**
     Verify that all points in trialPoints are on the current mesh.
     If a point is not on the mesh, issue a warning and remove it from the set.
     */
    virtual void verifyPointsAreOnMesh(const std::string& name);

    /// Snap a given trial point to the bounds and project on mesh
    /**
     * Used by classes that generate points: SearchMethods, Poll, etc,
     * to make the point satisfy the bounds and mesh requisites,
     * before sending it to evaluation.
     \param point        The point to process
     \param lowerBound   The lower bounds.
     \param upperBound   The upper bounds
     \param frameCenter  The frame center (can be null)
     \param mesh         The mesh (can be null)
     \return             \c true if the function worked, the point is now on mesh and inside bounds
     */
    static bool snapPointToBoundsAndProjectOnMesh(Point& point,
                                                  const ArrayOfDouble& lowerBound,
                                                  const ArrayOfDouble& upperBound,
                                                  const std::shared_ptr<Point> frameCenter,
                                                  const std::shared_ptr<MeshBase> mesh);

    /// Start evaluation of the trial points
    /**
     Called by run.
     \note Complete the documentation
     */
    bool evalTrialPoints(Step * step);

    /// Get the number of evaluation points in the queue for evaluation
    size_t getNbEvalPointsThatNeededEval() const { return _nbEvalPointsThatNeedEval ; }

    /// Generate the trial points of an algorithm iteration before evaluation.
    /** Virtual function that algorithm iteration steps must implement
     */
    virtual void generateTrialPoints() = 0;

    /// Add current frame center as originator of each point in trialPoints
    virtual void updatePointsWithFrameCenter() ;

private:

    /// Helper for constructor
    void init();


};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_ITERATIONUTILS__
