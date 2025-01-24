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

#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/QuadSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"

#ifdef USE_SGTELIB
#include "../../Algos/QuadModel/QuadModelSinglePass.hpp"
#endif

//
// Reference: File Sgtelib_Model_Search.cpp in NOMAD 3.9.1
// Author: Bastien Talgorn

void NOMAD::QuadSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_QUAD_MODEL);
    verifyParentNotNull();

    const auto parentSearch = getParentStep()->getParentOfType<NOMAD::QuadSearchMethod*>(false);

    // For some testing, it is possible that _runParams is null or evaluator control is null
    setEnabled((nullptr == parentSearch)
               && (nullptr !=_runParams)
               && (_runParams->getAttributeValue<bool>("QUAD_MODEL_SEARCH") || _runParams->getAttributeValue<bool>("QUAD_MODEL_SEARCH_SIMPLE_MADS"))
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


void NOMAD::QuadSearchMethod::generateTrialPointsFinal()
{
#ifdef USE_SGTELIB
    // The trial points are generated for a feasible frame center and an infeasible one.
    if ( ! _stopReasons->checkTerminate() )
    {
        auto madsIteration = getParentOfType<MadsIteration*>();

        // MegaIteration's barrier member is already in sub dimension.
        auto bestXIncFeas = madsIteration->getMegaIterationBarrier()->getCurrentIncumbentFeas();
        auto bestXIncInf  = madsIteration->getMegaIterationBarrier()->getCurrentIncumbentInf();
        const auto& computeType = madsIteration->getMegaIterationBarrier()->getFHComputeType();
        if (nullptr != bestXIncFeas
            && bestXIncFeas->getF(computeType).isDefined()
            && bestXIncFeas->getF(computeType) < MODEL_MAX_OUTPUT)
        {
            NOMAD::QuadModelSinglePass singlePassFeas(this, bestXIncFeas, madsIteration->getMesh(),{} /* no scaling direction */);

            // Generate the trial points
            singlePassFeas.generateTrialPoints();

            // Pass the generated trial pts to this
            const auto& trialPtsSinglePassFeas = singlePassFeas.getTrialPoints();
            for (auto evalPoint : trialPtsSinglePassFeas)
            {
                evalPoint.setPointFrom(bestXIncFeas, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
                evalPoint.addGenStep(getStepType(), false /*do not inherit -> just qms step*/);
                insertTrialPoint(evalPoint);
            }
        }
        if (nullptr != bestXIncInf
            && bestXIncInf->getF(computeType).isDefined()
            && bestXIncInf->getF(computeType) < MODEL_MAX_OUTPUT
            && bestXIncInf->getH(computeType).isDefined()
            && bestXIncInf->getH(computeType) < MODEL_MAX_OUTPUT)
        {
            NOMAD::QuadModelSinglePass singlePassInf(this, bestXIncInf, madsIteration->getMesh(),{} /* no scaling direction */);

            // Generate the trial points
            singlePassInf.generateTrialPoints();

            // Pass the generated trial pts to this
            const auto& trialPtsSinglePassInf = singlePassInf.getTrialPoints();
            for (auto evalPoint : trialPtsSinglePassInf)
            {
                evalPoint.setPointFrom(bestXIncInf, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
                evalPoint.addGenStep(getStepType(), false /*do not inherit-> just qms step */);
                insertTrialPoint(evalPoint);
            }
        }
    }
#endif
}
 // end generateTrialPoints
