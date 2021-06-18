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

#include "../Algos/EvcInterface.hpp"
#include "../Algos/SurrogateEvaluation.hpp"
#include "../Algos/SurrogateEvaluator.hpp"
#include "../Output/OutputQueue.hpp"

void NOMAD::SurrogateEvaluation::init()
{
    setStepType(NOMAD::StepType::SURROGATE_EVALUATION);
    verifyParentNotNull();
}


void NOMAD::SurrogateEvaluation::startImp()
{
}


bool NOMAD::SurrogateEvaluation::runImp()
{
    // Evaluation using static surrogate. The evaluation will be used for sorting afterwards.
    // Setup evaluation for SURROGATE:
    //  - Set opportunistic evaluation to false
    //  - Set the Evaluator to SURROGATE
    //  - The point's evalType will be set to SURROGATE in evalTrialPoints().
    // Evaluate the points using the surrogate
    // Reset for BB:
    //  - Reset opportunism
    //  - Reset Evaluator
    // And proceed - the sort using surrogate will be done afterwards.

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    auto previousOpportunism = evc->getOpportunisticEval();
    evc->setOpportunisticEval(false);

    auto surrogateEvaluator = std::make_shared<NOMAD::SurrogateEvaluator>(evc->getEvalParams());

    // Replace the EvaluatorControl's evaluator with this one
    // we just created
    auto previousEvaluator = evc->setEvaluator(surrogateEvaluator);
    if (nullptr == previousEvaluator)
    {
        std::cerr << "Warning: Could not set SURROGATE Evaluator" << std::endl;
        return false;
    }

    // get parent as IterationUtils, so that evalTrialPoints() can be used.
    auto stepAsIterationUtilsConst = dynamic_cast<const NOMAD::IterationUtils*>(getParentStep());
    auto stepAsIterationUtils = const_cast<NOMAD::IterationUtils*>(stepAsIterationUtilsConst);
    stepAsIterationUtils->evalTrialPoints(this);

    evc->setEvaluator(previousEvaluator);
    evc->setOpportunisticEval(previousOpportunism);

    // Points are now all evaluated using SURROGATE.
    // Points are still in the parent step's trial points. There is no update to do here.
    // They are ready to be sorted using their surrogate values, and then
    // evaluated using the BB.

    return true;
}


void NOMAD::SurrogateEvaluation::endImp()
{
}
