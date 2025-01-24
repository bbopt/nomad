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

#include "../../Algos/Mads/QPSolverAlgoSearchMethod.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverAlgoSinglePass.hpp"
#include "../../Algos/SubproblemManager.hpp"

void NOMAD::QPSolverAlgoSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_QUAD_MODEL);
    verifyParentNotNull();

    const auto parentSearch = getParentStep()->getParentOfType<NOMAD::QPSolverAlgoSearchMethod*>(false);

    // For some testing, it is possible that _runParams is null or evaluator control is null
    setEnabled((nullptr == parentSearch)
               && (nullptr !=_runParams)
               && _runParams->getAttributeValue<bool>("QP_SEARCH")
               &&  (nullptr != EvcInterface::getEvaluatorControl()));
#ifndef USE_SGTELIB
    if (isEnabled())
    {
        OUTPUT_INFO_START
        AddOutputInfo(getName() + " cannot be performed because NOMAD is compiled without sgtelib library");
        OUTPUT_INFO_END
        setEnabled(false);
    }
#endif

#ifdef USE_SGTELIB
    // Check that there is exactly one objective
    if (isEnabled())
    {
        auto nbObj = NOMAD::Algorithm::getNbObj();
        if (0 == nbObj)
        {
            OUTPUT_INFO_START
            AddOutputInfo(getName() + " not performed when there is no objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }
        else if (nbObj > 1)
        {
            OUTPUT_INFO_START
            AddOutputInfo(getName() + " not performed on multi-objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }

        const auto& modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");
        _displayLevel = modelDisplay.empty()
                            ? NOMAD::OutputLevel::LEVEL_DEBUGDEBUG
                            : NOMAD::OutputLevel::LEVEL_INFO;
    }
#endif
}


void NOMAD::QPSolverAlgoSearchMethod::generateTrialPointsFinal()
{
#ifdef USE_SGTELIB
    // The trial points are generated for a feasible frame center and an infeasible one.
    if ( ! _stopReasons->checkTerminate() )
    {
        auto madsIteration = getParentOfType<MadsIteration*>();

        // MegaIteration's barrier member is already in sub dimension.
        auto bestXFeas = madsIteration->getMegaIterationBarrier()->getCurrentIncumbentFeas();
        auto bestXInf  = madsIteration->getMegaIterationBarrier()->getCurrentIncumbentInf();

        const auto& computeType = madsIteration->getMegaIterationBarrier()->getFHComputeType();
        if (nullptr != bestXFeas
            && bestXFeas->getF(computeType).isDefined()
            && bestXFeas->getF(computeType) < MODEL_MAX_OUTPUT)
        {
            NOMAD::QPSolverAlgoSinglePass singlePassFeas(this, bestXFeas, madsIteration->getMesh(),{} /* no scaling direction */);

            // Generate the trial points
            singlePassFeas.generateTrialPoints();

            // Pass the generated trial pts to this
            const auto& trialPtsSinglePassFeas = singlePassFeas.getTrialPoints();
            for (auto evalPoint : trialPtsSinglePassFeas)
            {
                evalPoint.setPointFrom(bestXFeas, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
                insertTrialPoint(evalPoint);
            }
        }
        if (nullptr != bestXInf
            && bestXInf->getF(computeType).isDefined()
            && bestXInf->getF(computeType) < MODEL_MAX_OUTPUT
            && bestXInf->getH(computeType).isDefined()
            && bestXInf->getH(computeType) < MODEL_MAX_OUTPUT)
        {
            NOMAD::QPSolverAlgoSinglePass singlePassInf(this, bestXInf, madsIteration->getMesh(),{} /* no scaling direction */);

            // Generate the trial points
            singlePassInf.generateTrialPoints();

            // Pass the generated trial pts to this
            const auto& trialPtsSinglePassInf = singlePassInf.getTrialPoints();
            for (auto evalPoint : trialPtsSinglePassInf)
            {
                evalPoint.setPointFrom(bestXInf, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
                insertTrialPoint(evalPoint);
            }
        }
    }
#endif
}

