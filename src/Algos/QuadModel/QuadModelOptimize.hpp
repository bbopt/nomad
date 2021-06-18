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
#ifndef __NOMAD_4_0_QUAD_MODEL_OPTIMIZE__
#define __NOMAD_4_0_QUAD_MODEL_OPTIMIZE__

#include "../../Algos/Step.hpp"
#include "../../Algos/QuadModel/QuadModelIterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to create trial points by performing quadratic model optimization using Mads
/**
 - Start, run and end tasks are performed.
 - Start: the quadratic model optimization problem is setup and solved by calling startImp. Call ::generateTrialPoints.
 - Run: trial (oracle) points are evaluated with EvalType::BB. Set the stop reason.
 - End: Remove from cache EvalType::MODEL only cache points.
 */
class QuadModelOptimize : public Step, public QuadModelIterationUtils
{
private:

    OutputLevel                         _displayLevel;

    ArrayOfDouble _modelLowerBound; ///> Lower bound: min of trainingSet points
    ArrayOfDouble _modelUpperBound; ///> Upper bound: max of trainingSet points
    Point         _modelFixedVar;   ///> Fixed variables: fixed variables detected from trainingSet

    Point _modelCenter;

    const std::shared_ptr<PbParameters> _refPbParams; ///< Reference to the original problem parameters.

    std::shared_ptr<RunParameters>      _optRunParams; ///< run parameters for model optimization
    std::shared_ptr<PbParameters>       _optPbParams; ///< pb parameters for model optimization


public:
    /// Constructor
    /* Parent must explicitely be a (pointer to a) QuadModelAlgo.
     * Run parameters will be recomputed for model optimization.
     */
    explicit QuadModelOptimize(const Step* parentStep,
                               const std::shared_ptr<PbParameters>               refPbParams )
      : Step(parentStep),
      QuadModelIterationUtils (parentStep),
        _displayLevel(OutputLevel::LEVEL_INFO),
        _modelLowerBound(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelUpperBound(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelFixedVar(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelCenter(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _refPbParams(refPbParams),
        _optRunParams(nullptr),
        _optPbParams(nullptr)
    {
        init();
    }

    /// Generate new points to evaluate
    /**
     - Setup the evaluator control parameters.
     - Manage display of sub-optimization.
     - Setup evaluator (EvalType::MODEL) and success type identification function.
     - Setup the bounds and fixed variables from the trainingSet of the quadratic model.
     - Setup run and pb parameters for Mads
     - Perform start, run and end tasks on Mads.
     - best feasible and best infeasible (if available) are inserted as trial points.
     */
    void generateTrialPoints() override;


private:
    void init();

    virtual void startImp() override; ///< The quadratic model optimization problem is setup and solved by calling startImp. Calls ::generateTrialPoints.
    virtual bool runImp() override; ///< Trial (oracle) points are evaluated with EvalType::BB. Set the stop reason.
    virtual void endImp() override; ///< Remove from cache EvalType::MODEL only cache points.

    // Helpers
    void setupRunParameters();
    void setupPbParameters();
    void setModelBoundsAndFixedVar();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_QUAD_MODEL_OPTIMIZE__
