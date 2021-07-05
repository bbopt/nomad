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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelEvaluator.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelOptimize.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"
#include "../../Type/LHSearchType.hpp"

void NOMAD::SgtelibModelOptimize::init()
{
    setStepType(NOMAD::StepType::OPTIMIZE);
    verifyParentNotNull();

    // Set and verify run parameter values
    setupRunParameters();
    // Pb Parameters are set and verified later, they need arguments
}


void NOMAD::SgtelibModelOptimize::startImp()
{
    auto modelDisplay = _runParams->getAttributeValue<std::string>("SGTELIB_MODEL_DISPLAY");
    _displayLevel = (std::string::npos != modelDisplay.find("O"))
                        ? NOMAD::OutputLevel::LEVEL_INFO
                        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

    OUTPUT_INFO_START
    std::string s;
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();
    s = "SGTELIB_MODEL_MAX_EVAL: " + std::to_string(evcParams->getAttributeValue<size_t>("SGTELIB_MODEL_MAX_EVAL"));
    AddOutputInfo(s, _displayLevel);
    s = "BBOT: " + NOMAD::BBOutputTypeListToString(NOMAD::SgtelibModel::getBBOutputType());
    AddOutputInfo(s, _displayLevel);
    s = "Formulation: " + NOMAD::SgtelibModelFormulationTypeToString(_runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION"));
    AddOutputInfo(s, _displayLevel);

    std::ostringstream oss;
    oss << "Run Parameters for SgtelibModelOptimize:" << std::endl;
    _optRunParams->display(oss, false);
    AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
    OUTPUT_INFO_END
}


bool NOMAD::SgtelibModelOptimize::runImp()
{
    bool optimizeOk = false;
    std::string s;

    auto modelFormulation = _runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION");
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    if (NOMAD::SgtelibModelFormulationType::EXTERN == modelFormulation)
    {
        std::string modelDefinition = _runParams->getAttributeValue<std::string>("MODEL_DEFINITION");
        evc->getEvalParams()->setAttributeValue("BB_EXE", modelDefinition);
    }
    else
    {
        // Enforce no opportunism and use no cache.
        auto previousOpportunism = evc->getOpportunisticEval();
        auto previousUseCache = evc->getUseCache();
        evc->setOpportunisticEval(false);
        evc->setUseCache(false);

        auto modelDisplay = _runParams->getAttributeValue<std::string>("SGTELIB_MODEL_DISPLAY");
        auto diversification = _runParams->getAttributeValue<NOMAD::Double>("SGTELIB_MODEL_DIVERSIFICATION");
        auto modelFeasibility = _runParams->getAttributeValue<NOMAD::SgtelibModelFeasibilityType>("SGTELIB_MODEL_FEASIBILITY");
        double tc = _runParams->getAttributeValue<NOMAD::Double>("SGTELIB_MODEL_SEARCH_EXCLUSION_AREA").todouble();

        if (nullptr == _modelAlgo)
        {
            s = "Error: In SgtelibModelOptimize, need a SgtelibModel parent.";
            throw NOMAD::Exception(__FILE__, __LINE__, s);
        }

        auto ev = std::make_shared<NOMAD::SgtelibModelEvaluator>(
                                                evc->getEvalParams(), _modelAlgo, modelDisplay,
                                                diversification, modelFeasibility, tc,
                                                NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));

        // Replace the EvaluatorControl's evaluator with this one
        // we just created
        auto previousEvaluator = evc->setEvaluator(ev);
        if (nullptr == previousEvaluator)
        {
            std::cerr << "Warning: QuadModelOptimize: Could not set MODEL Evaluator" << std::endl;
            return false;
        }

        auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

        // Create a Mads step
        // Parameters for mads (_optRunParams and _optPbParams) are already updated.
        _mads = std::make_shared<NOMAD::Mads>(this, madsStopReasons, _optRunParams, _optPbParams);
        //_mads->setName(_mads->getName() + " (SgtelibModelOptimize)");
        _mads->setEndDisplay(false);
        evc->resetModelEval();
        _mads->start();
        optimizeOk = _mads->run();
        _mads->end();
        evc->resetModelEval();
        evc->setEvaluator(previousEvaluator);

        // Note: No need to check the Mads stop reason: It is not a stop reason
        // for SgtelibModel.

        // Get the solutions
        updateOraclePoints();

        // Reset opportunism to previous values.
        evc->setOpportunisticEval(previousOpportunism);
        evc->setUseCache(previousUseCache);
    }

    if (!optimizeOk)
    {
        auto sgetlibModelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        sgetlibModelStopReasons->set(NOMAD::ModelStopType::MODEL_OPTIMIZATION_FAIL);
    }

    return optimizeOk;
}


void NOMAD::SgtelibModelOptimize::endImp()
{
}


void NOMAD::SgtelibModelOptimize::setupRunParameters()
{
    _optRunParams = std::make_shared<NOMAD::RunParameters>(*_refRunParams);

    // Ensure there is no model used in model optimization.
    _optRunParams->setAttributeValue("SGTELIB_MODEL_SEARCH", false);
    _optRunParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    _optRunParams->setAttributeValue("VNS_MADS_SEARCH", false);

    // Set direction type to Ortho 2n
    _optRunParams->setAttributeValue("DIRECTION_TYPE",NOMAD::DirectionType::ORTHO_2N);

    // Use isotropic mesh
    _optRunParams->setAttributeValue("ANISOTROPIC_MESH", false);

    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();
    std::string lhstr = std::to_string(int(evcParams->getAttributeValue<size_t>("SGTELIB_MODEL_MAX_EVAL") * 0.3));
    lhstr += " 0";
    NOMAD::LHSearchType lhSearch(lhstr);
    _optRunParams->setAttributeValue("LH_SEARCH", lhSearch);

    // No hMax in the context of SgtelibModel
    _optRunParams->setAttributeValue("H_MAX_0", NOMAD::Double(NOMAD::INF));

    // Disable user calls
    _optRunParams->setAttributeValue("USER_CALLS_ENABLED", false);

    _optRunParams->checkAndComply(evcParams, _pbParams);
}


void NOMAD::SgtelibModelOptimize::setupPbParameters(const NOMAD::ArrayOfDouble& lowerBound,
                                                    const NOMAD::ArrayOfDouble& upperBound)
{
    _optPbParams = std::make_shared<NOMAD::PbParameters>(*_refPbParams);

    _optPbParams->setAttributeValue("LOWER_BOUND", lowerBound);
    _optPbParams->setAttributeValue("UPPER_BOUND", upperBound);

    // Reset initial mesh and frame sizes
    // The initial mesh and frame sizes will be calculated from bounds and X0
    _optPbParams->resetToDefaultValue("INITIAL_MESH_SIZE");
    _optPbParams->resetToDefaultValue("INITIAL_FRAME_SIZE");
    // Use default min mesh and frame sizes
    _optPbParams->resetToDefaultValue("MIN_MESH_SIZE");
    _optPbParams->resetToDefaultValue("MIN_FRAME_SIZE");

    // Granularity is set to 0 and bb_input_type is set to all continuous variables. Candidate points are projected on the mesh before evaluation.
    _optPbParams->resetToDefaultValue("GRANULARITY");
    _optPbParams->resetToDefaultValue("BB_INPUT_TYPE");

    // No variable groups are considered for suboptimization
    _optPbParams->resetToDefaultValue("VARIABLE_GROUP");

    // Find best points (MODEL evals) and use them as X0s to optimize models.
    NOMAD::CacheInterface cacheInterface(this);
    std::vector<NOMAD::EvalPoint> evalPointFeasList;
    std::vector<NOMAD::EvalPoint> evalPointInfList;
    NOMAD::Double hMax = _modelAlgo->getHMax();

    // Only looking into model evaluations here
    cacheInterface.findBestFeas(evalPointFeasList,
                                NOMAD::EvalType::MODEL,
                                NOMAD::ComputeType::STANDARD,
                                nullptr);
    cacheInterface.findBestInf(evalPointInfList,
                               hMax,
                               NOMAD::EvalType::MODEL,
                               NOMAD::ComputeType::STANDARD,
                               nullptr);

    NOMAD::ArrayOfPoint x0s;
    for (auto evalPointX0 : evalPointFeasList)
    {
        if (evalPointX0.inBounds(lowerBound, upperBound))
        {
            x0s.push_back(*(evalPointX0.getX()));
        }
    }
    for (auto evalPointX0 : evalPointInfList)
    {
        if (evalPointX0.inBounds(lowerBound, upperBound))
        {
            x0s.push_back(*(evalPointX0.getX()));
        }
    }
    // Fallback: No MODEL points found. Use points from an upper Mads barrier
    // (_barrierForX0s).
    if (0 == x0s.size())
    {
        // Get best points from upper Mads
        for (auto evalPointX0 : _modelAlgo->getX0s())
        {
            if (evalPointX0.inBounds(lowerBound, upperBound))
            {
                x0s.push_back(*(evalPointX0.getX()));
            }
        }
    }
    _optPbParams->setAttributeValue("X0", x0s);

    // We do not want certain warnings appearing in sub-optimization.
    _optPbParams->doNotShowWarnings();

    _optPbParams->checkAndComply();

}


// Oracle points are the best feasible and infeasible points found by
// running Mads member. Get them using the barrier.
void NOMAD::SgtelibModelOptimize::updateOraclePoints()
{
    _oraclePoints.clear();

    std::shared_ptr<NOMAD::Barrier> barrier;
    if (nullptr != _mads && nullptr != _mads->getMegaIteration())
    {
        barrier = _mads->getMegaIteration()->getBarrier();
    }

    if (barrier)
    {
        auto allBestPoints = barrier->getAllPoints();

        for (auto evalPoint : allBestPoints)
        {
            _oraclePoints.insert(evalPoint);
        }
    }

}


