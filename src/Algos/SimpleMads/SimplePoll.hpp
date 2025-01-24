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
#ifndef __NOMAD_4_5_SIMPLEPOLL__
#define __NOMAD_4_5_SIMPLEPOLL__

#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/Mads/PollMethodBase.hpp"
#include "../../Algos/QuadModel/QuadModelEvaluator.hpp"
#include "../../Algos/SimpleMads/SimpleProgressiveBarrier.hpp"

#include "../../../ext/sgtelib/src/Surrogate.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the poll (step 3) of MADS algorithm.
/**
 Generate the trial points (Poll::startImp), launch evaluation (Poll::runImp) and postprocessing (Poll::endImp).
 */
class SimplePoll: public Iteration
{
private:

    std::vector<NOMAD::SimpleEvalPoint> _trialPoints; ///< The points generated during the start(). Used for run() and postProcessing().


    NOMAD::BBOutputTypeList _bbot;

    std::unique_ptr<SimpleProgressiveBarrier> _barrier;
    std::shared_ptr<GMesh> _mesh;

    std::shared_ptr<SGTELIB::Surrogate> _model;

    const singleOutputComputeFType & _singleObjCompute;

    std::vector<SimpleEvalPoint> _frameCenters;  ///< The frame centers (primary and secondary) of the poll methods.

    DirectionType _primaryDirectionType, _secondaryDirectionType;  ///< The poll methods implement different direction types for primary and secondary poll centers.

    NOMAD::Double _rho; ///< Rho parameter of the progressive barrier. Used to choose if the primary frame center is the feasible or infeasible  incumbent.

    size_t _nbEval;

    std::vector<std::shared_ptr<PollMethodBase>> _pollMethods;  ///< Unlike for Search, Poll methods generate all their points and only then they are evaluated.

    size_t _n, _nSimple, _nbOutputs;  ///< Pb dimension

    Point _fixedVariable;

    bool _phaseOneSearch;

    std::function<bool(std::vector<NOMAD::SimpleEvalPoint>&)> _eval_x; ///< Function for outputs evaluation

public:
    /// Constructor #1 for model optim
    /**
     \param parentStep The parent of this poll step
     */
    explicit SimplePoll(const Step* parentStep, const std::shared_ptr<SGTELIB::Surrogate> & model, const NOMAD::BBOutputTypeList & bbot, const singleOutputComputeFType & singleObjCompute)
      : Iteration(parentStep, 0),
        _barrier(nullptr),
        _mesh(nullptr),
        _model(model),
        _bbot(bbot),
        _nbEval(0),
        _phaseOneSearch(false),
        _singleObjCompute(singleObjCompute)
    {
        init();
    }
    /// Constructor #2, for user function optim
    /**
     \param parentStep The parent of this poll step
     */
    explicit SimplePoll(const Step* parentStep, const NOMAD::BBOutputTypeList & bbot, std::function<bool(std::vector<NOMAD::SimpleEvalPoint> &)> eval_x)
      : Iteration(parentStep, 0),
        _barrier(nullptr),
        _mesh(nullptr),
        _bbot(bbot),
        _nbEval(0),
        _phaseOneSearch(false),
        _eval_x(eval_x),
        _singleObjCompute(NOMAD::defaultEmptySingleOutputCompute)
    {
        init();
    }

    virtual ~SimplePoll() {}

    const std::unique_ptr<SimpleProgressiveBarrier> & getBarrier() const { return _barrier; }

    const MeshBasePtr getMesh() const override { return _mesh; }

    size_t getNbEval () const { return _nbEval;}

    bool getPhaseOneSearch () const { return _phaseOneSearch; }

protected:
    /// Helper for start: get lists of Primary and Secondary Polls
    void computePrimarySecondaryPollCenters(SimpleEvalPoint & primaryCenter, SimpleEvalPoint & secondaryCenter) const;

    /// Helper for start: create a poll method
    // virtual void createPollMethod(const bool isPrimary, const EvalPointPtr frameCenter);
    virtual void createPollMethod(const bool isPrimary, const SimpleEvalPoint & frameCenter);

    /// Helper to create poll methods for current poll centers
    virtual void createPollMethodsForPollCenters();


private:
    /// Helper for constructor
    void init();

    /// Helpers
    void generateTrialPoints();
    void evalTrialPoints();

    NOMAD::Double getF(const NOMAD::ArrayOfDouble & out) const;
    NOMAD::Double getH(const NOMAD::ArrayOfDouble & out) const;

    /// Implementation for start tasks for MADS poll.
    /**
     Call to generate the poll methods
     */
    virtual void    startImp()  override;

    /// Implementation for run tasks for MADS poll.
    /**
     Call poll methods and perform trial points evaluation.
     \return Flag \c true if found better solution \c false otherwise.
     */
    virtual bool    runImp() override;

    /// Implementation for end tasks for MADS simple poll.
    /**
     Call the IterationUtils::postProcessing of the points.
     */
    virtual void    endImp() override ;



};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_SIMPLEPOLL__
