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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/PSDMads/PSDMadsMegaIteration.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"
#include "../../Type/LHSearchType.hpp"


void NOMAD::PSDMadsMegaIteration::destroy()
{
    _madsOnSubPb.reset();
    setStepType(NOMAD::StepType::MEGA_ITERATION);
}


void NOMAD::PSDMadsMegaIteration::startImp()
{
    auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
    bool isPollster = (0 == NOMAD::getThreadNum());

    // Set parameters for subproblems
    auto subProblemPbParams = std::make_shared<NOMAD::PbParameters>(*_pbParams);
    auto subProblemRunParams = std::make_shared<NOMAD::RunParameters>(*_runParams);
    setupSubproblemParams(subProblemPbParams, subProblemRunParams, isPollster);

    // Create Mads for this subproblem
    _madsOnSubPb = std::make_shared<NOMAD::Mads>(this, madsStopReasons, subProblemRunParams, subProblemPbParams);
    /*
    std::string madsName = "Mads ";
    if (isPollster)
    {
        madsName += "pollster";
    }
    else
    {
        if (_fixedVariable.size() <= 10)
        {
            madsName += "with fixed variable ";
            madsName += _fixedVariable.display();
        }
        else
        {
            madsName += "with ";
            madsName += NOMAD::itos(_fixedVariable.size() - _fixedVariable.nbDefined());
            madsName += " fixed variables";
        }
    }
    _madsOnSubPb->setName(madsName);
    */
    _madsOnSubPb->setStepType(NOMAD::StepType::ALGORITHM_PSD_MADS_SUBPROBLEM);
}



bool NOMAD::PSDMadsMegaIteration::runImp()
{
    // Run Mads
    // Note: Pollster is always run whenever thread 0 is available.
    // However, contrary to the NOMAD 3 version, mesh is not updated at each pollster run.

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    int mainThreadNum = NOMAD::getThreadNum();

    std::string s = "Running " + _madsOnSubPb->getName();
    s += " on thread " + NOMAD::itos(mainThreadNum);
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_NORMAL);
    _madsOnSubPb->start();
    bool madsSuccessful = _madsOnSubPb->run();   // If this run is successful, barrier will be updated.
    _madsOnSubPb->end();

    s = "Done running " + _madsOnSubPb->getName();
    s += " on thread " + NOMAD::itos(mainThreadNum) + ". ";
    s += "Number of evaluations: " + NOMAD::itos(evc->getBbEvalInSubproblem()) + ". ";
    s += "Found a new success: " + NOMAD::boolToString(madsSuccessful) + ".";
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_NORMAL);
    evc->resetBbEvalInSubproblem();

    return madsSuccessful;
}


void NOMAD::PSDMadsMegaIteration::setupSubproblemParams(std::shared_ptr<NOMAD::PbParameters> &subProblemPbParams,
                                           std::shared_ptr<NOMAD::RunParameters> &subProblemRunParams,
                                           const bool isPollster)
{
    auto mainFrameSize = _mainMesh->getDeltaFrameSize();
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    // Note: If n >= 50, models are disabled. They could be re-enabled on
    // subproblems with lesser dimension. See issue #370.

    subProblemPbParams->doNotShowWarnings();
    if (isPollster)
    {
        subProblemRunParams->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::SINGLE);
        subProblemPbParams->setAttributeValue("INITIAL_FRAME_SIZE", mainFrameSize);

        // Disable all searches
        subProblemRunParams->setAttributeValue("LH_SEARCH", NOMAD::LHSearchType("0 0"));
        subProblemRunParams->setAttributeValue("NM_SEARCH", false);
        subProblemRunParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
        subProblemRunParams->setAttributeValue("SGTELIB_MODEL_SEARCH", false);
        subProblemRunParams->setAttributeValue("SPECULATIVE_SEARCH", false);
        subProblemRunParams->setAttributeValue("VNS_MADS_SEARCH", false);  // VNS has static member. Problematic with threads. See issue # 604
        
    }
    else
    {
        auto initialFrameSize = _mainMesh->getDeltaFrameSizeCoarser();
        subProblemPbParams->setAttributeValue("INITIAL_FRAME_SIZE", initialFrameSize);

        // The main frame size is used as minFrameSize for the subproblem.
        // Initial and min must be compatible -> adjust.
        for (size_t i = 0; i < initialFrameSize.size(); i++)
        {
            if (initialFrameSize[i] < mainFrameSize[i])
            {
                OUTPUT_INFO_START
                AddOutputInfo("Set initial frame size to main frame size.");
                OUTPUT_INFO_END
                subProblemPbParams->setAttributeValue("INITIAL_FRAME_SIZE", mainFrameSize);
                break;
            }
        }

        subProblemPbParams->setAttributeValue("FIXED_VARIABLE", _fixedVariable);
        subProblemPbParams->setAttributeValue("X0", _x0);
        subProblemPbParams->setAttributeValue("MIN_FRAME_SIZE", mainFrameSize);
        subProblemRunParams->setAttributeValue("PSD_MADS_NB_VAR_IN_SUBPROBLEM", _fixedVariable.size() - _fixedVariable.nbDefined());
    }

    // Set max number of bb evals per subproblem.
    // Strategy to be refined.
    if (isPollster)
    {
        evc->setMaxBbEvalInSubproblem(1);
    }
    else
    {
        size_t totalBudget = evc->getEvaluatorControlGlobalParams()->getAttributeValue<size_t>("MAX_BB_EVAL");
        size_t nbSubproblem = _runParams->getAttributeValue<size_t>("PSD_MADS_NB_SUBPROBLEM");
        size_t budgetPerSubproblem = (totalBudget-1) / nbSubproblem;
        size_t maxBbEvalInSubproblem = evc->getMaxBbEvalInSubproblem();
        if (budgetPerSubproblem < maxBbEvalInSubproblem)
        {
            maxBbEvalInSubproblem = budgetPerSubproblem;
        }
        evc->setMaxBbEvalInSubproblem(maxBbEvalInSubproblem);
    }

    subProblemPbParams->checkAndComply();
    auto evcParams = evc->getEvaluatorControlGlobalParams();
    subProblemRunParams->checkAndComply(evcParams, subProblemPbParams);
}
