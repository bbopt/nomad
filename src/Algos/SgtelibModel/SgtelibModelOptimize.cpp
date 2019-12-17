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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelEvaluator.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelOptimize.hpp"
#include "../../Type/LHSearchType.hpp"

void NOMAD::SgtelibModelOptimize::init()
{
    _name = getAlgoName() + "Optimize";
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

    std::string s;
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlParams();
    s = "SGTELIB_MODEL_EVAL_NB: " + std::to_string(evcParams->getAttributeValue<size_t>("SGTELIB_MODEL_EVAL_NB"));
    AddOutputInfo(s, _displayLevel);
    s = "BBOT: " + NOMAD::BBOutputTypeListToString(NOMAD::SgtelibModel::getBBOutputType());
    AddOutputInfo(s, _displayLevel);
    s = "Formulation: " + NOMAD::SgtelibModelFormulationTypeToString(_runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION"));
    AddOutputInfo(s, _displayLevel);

    std::ostringstream oss;
    oss << "Run Parameters for SgtelibModelOptimize:" << std::endl;
    _optRunParams->display(oss, false);
    AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
}


bool NOMAD::SgtelibModelOptimize::runImp()
{
    bool optimizeOk = false;
    std::string s;

    auto modelFormulation = _runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION");
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlParams();


    if (NOMAD::SgtelibModelFormulationType::EXTERN == modelFormulation)
    {
        std::string modelDefinition = _runParams->getAttributeValue<std::string>("MODEL_DEFINITION");
        evcParams->setAttributeValue("BB_EXE", modelDefinition);
    }
    else
    {
        auto evalParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvalParams();

        auto modelDisplay = _runParams->getAttributeValue<std::string>("SGTELIB_MODEL_DISPLAY");
        auto diversification = _runParams->getAttributeValue<NOMAD::Double>("SGTELIB_MODEL_DIVERSIFICATION");
        auto modelFeasibility = _runParams->getAttributeValue<NOMAD::SgtelibModelFeasibilityType>("SGTELIB_MODEL_FEASIBILITY");
        double tc = _runParams->getAttributeValue<NOMAD::Double>("SGTELIB_MODEL_EXCLUSION_AREA").todouble();

        if (nullptr == _modelAlgo)
        {
            s = "Error: In SgtelibModelOptimize, need a SgtelibModel parent.";
            throw NOMAD::Exception(__FILE__, __LINE__, s);
        }

        auto ev = std::make_unique<NOMAD::SgtelibModelEvaluator>(                                                                                  evalParams, _modelAlgo, modelDisplay,
                                                diversification, modelFeasibility, tc,
                                                getSubFixedVariable());

        // Replace the EvaluatorControl's evaluator with this one
        // we just created
        auto mainEvaluator = EvcInterface::getEvaluatorControl()->getEvaluatorUPtr();
        EvcInterface::getEvaluatorControl()->setEvaluator(std::move(ev));

        auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

        // Create a Mads step
        // Parameters for mads (_optRunParams and _optPbParams) are already updated.
        _mads = std::make_shared<NOMAD::Mads>(this, madsStopReasons, _optRunParams, _optPbParams);
        _mads->setName(_mads->getName() + " (SgtelibModelOptimize)");
        EvcInterface::getEvaluatorControl()->resetSgteEval();
        NOMAD::MainStep::setAlgoComment("(SgtelibModelOptimize)");
        _mads->start();
        optimizeOk = _mads->run();
        _mads->end();
        NOMAD::MainStep::resetPreviousAlgoComment();

        // Note: No need to check the Mads stop reason: It is not a stop reason
        // for SgtelibModel.
        
        // Get the solutions
        updateOraclePoints();

        // When we are done, restore mainEvaluator
        EvcInterface::getEvaluatorControl()->setEvaluator(std::move(mainEvaluator));
    }

    if (!optimizeOk)
    {
        auto SgtelibModelStopReasons = NOMAD::AlgoStopReasons<NOMAD::SgtelibModelStopType>::get(_stopReasons);
        SgtelibModelStopReasons->set(NOMAD::SgtelibModelStopType::MODEL_OPTIMIZER_FAIL);
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
    _optRunParams->setAttributeValue("SGTELIB_SEARCH", false);
    _optRunParams->setAttributeValue("DISABLE", std::string("MODELS"));

    // Use isotropic mesh
    _optRunParams->setAttributeValue("ANISOTROPIC_MESH", false);

    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlParams();
    std::string lhstr = std::to_string(int(evcParams->getAttributeValue<size_t>("SGTELIB_MODEL_EVAL_NB") * 0.3));
    lhstr += " 0";
    NOMAD::LHSearchType lhSearch(lhstr);
    _optRunParams->setAttributeValue("LH_SEARCH", lhSearch);

    // Disable user calls
    _optRunParams->setAttributeValue("USER_CALLS_ENABLED", false);

    _optRunParams->checkAndComply(evcParams, _pbParams);
}


void NOMAD::SgtelibModelOptimize::setupPbParameters(const NOMAD::ArrayOfDouble& lowerBound,
                                                    const NOMAD::ArrayOfDouble& upperBound,
                                                    const NOMAD::ArrayOfDouble& initialMeshSize,
                                                    const NOMAD::ArrayOfDouble& initialFrameSize)
{
    _optPbParams = std::make_shared<NOMAD::PbParameters>(*_refPbParams);

    _optPbParams->setAttributeValue("LOWER_BOUND", lowerBound);
    _optPbParams->setAttributeValue("UPPER_BOUND", upperBound);

    // Find best points (SGTE evals) and use them as X0s to optimize models.
    NOMAD::CacheInterface cacheInterface(this);
    std::vector<NOMAD::EvalPoint> evalPointFeasList;
    std::vector<NOMAD::EvalPoint> evalPointInfList;
    NOMAD::Double hMax = _modelAlgo->getHMax();

    // Only looking into sgte evaluations here
    cacheInterface.findBestFeas(evalPointFeasList,
                                getSubFixedVariable(),
                                NOMAD::EvalType::SGTE);
    cacheInterface.findBestInf(evalPointInfList,
                               hMax,
                               getSubFixedVariable(),
                               NOMAD::EvalType::SGTE);

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
    // Fallback: No SGTE points found. Use points from an upper Mads barrier
    // (_barrierForX0s).
    if (0 == x0s.size())
    {
        // Get best points from upper Mads 
        for (auto evalPointX0 : _modelAlgo->getX0s())
        {
            if (evalPointX0->inBounds(lowerBound, upperBound))
            {
                x0s.push_back(*(evalPointX0->getX()));
            }
        }
    }
    _optPbParams->setAttributeValue("X0", x0s);

    if (initialMeshSize.isDefined() && initialMeshSize.isComplete())
    {
        _optPbParams->setAttributeValue("INITIAL_MESH_SIZE", initialMeshSize);
    }
    if (initialFrameSize.isDefined() && initialFrameSize.isComplete())
    {
        _optPbParams->setAttributeValue("INITIAL_FRAME_SIZE", initialFrameSize);
    }

    // Ensure that optimization will be done using SGTE evals. All sub
    // steps will use SGTE.
    _optPbParams->setAttributeValue("EVAL_TYPE", NOMAD::EvalType::SGTE);

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
        auto allBestFeas = barrier->getAllXFeas();
        auto allBestInf  = barrier->getAllXInf();

        for (auto bestFeas : allBestFeas)
        {
            _oraclePoints.insert(*bestFeas);
        }
        for (auto bestInf : allBestInf)
        {
            _oraclePoints.insert(*bestInf);
        }
    }

}


