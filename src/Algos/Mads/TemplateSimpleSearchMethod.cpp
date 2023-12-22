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
/**
 \file   TemplateSimpleSearchMethod.cpp
 \brief  A template for simple random search without iteration (implementation)
 \author Christophe Tribes
 \date   2022-05-25
 */

#include "../../Algos/Mads/TemplateSimpleSearchMethod.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoSinglePass.hpp"
#include "../../Output/OutputQueue.hpp"

/*-------------------------------------------------------------*/
/*    Template for a simple (no iterations) search method      */
/*-------------------------------------------------------------*/
/*  Can be called (RANDOM_SIMPLE_SEARCH yes)                   */
/*  Generate random points around best incumbent.              */
/*  TEMPLATE use for a new search method: copy and rename the  */
/*  file and the class name. Adapt the code to your needs.     */
/*-------------------------------------------------------------*/

void NOMAD::TemplateSimpleSearchMethod::init()
{
    // TEMPLATE use for a new search method: define a specific step type in ../Type/StepType.hpp.
    setStepType(NOMAD::StepType::SEARCH_METHOD_SIMPLE_RANDOM);

    bool enabled = false;
    // For some testifng, it is possible that _runParams is null
    if (nullptr != _runParams)
    {
        // TEMPLATE use for a new search method: a new parameter must be defined to enable or not the search method (see ../Attributes/runAttributesDefinition.txt)
        enabled = _runParams->getAttributeValue<bool>("RANDOM_SIMPLE_SEARCH");
    }

    setEnabled(enabled);
}


void NOMAD::TemplateSimpleSearchMethod::generateTrialPointsFinal()
{
    // The trial points of one iteration of this template algorithm (random) are generated (not evaluated).

    // Note: Use first point of barrier as center.
    NOMAD::TemplateAlgoSinglePass randomAlgo(this,
                                             getMegaIterationBarrier()->getFirstPoint());
    randomAlgo.start();
    randomAlgo.end();

    // Pass the generated trial pts to this
    auto trialPtsNM = randomAlgo.getTrialPoints();
    for (auto point : trialPtsNM)
    {
        insertTrialPoint(point);
    }
    
}
