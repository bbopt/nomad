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

#include "../../Algos/Mads/TemplateAlgoSearchMethod.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgo.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoSinglePass.hpp"

void NOMAD::TemplateAlgoSearchMethod::init()
{
    // For some testing, it is possible that _runParams is null or evaluator control is null
    bool randomAlgoSearch = false;
    if ( nullptr != _runParams && nullptr != NOMAD::EvcInterface::getEvaluatorControl() )
    {
        if ( _runParams->getAttributeValue<bool>("MEGA_SEARCH_POLL") )
        {
            setStepType(NOMAD::StepType::SEARCH_METHOD_ALGO_RANDOM);
        }
        else
        {
            setStepType(NOMAD::StepType::ALGORITHM_RANDOM);
        }
        // TEMPLATE use for a new search method: a new parameter must be defined to enable or not the search method (see ../Attributes/runAttributesDefinition.txt)
        randomAlgoSearch = _runParams->getAttributeValue<bool>("RANDOM_ALGO_SEARCH");
    }
    setEnabled(randomAlgoSearch);
    
    
    if (randomAlgoSearch)
    {
        // TEMPLATE for a new search method: parameters can be defined to control the search method (see ../Attributes/runAttributesDefinition.txt)
        auto dummyFactor = _runParams->getAttributeValue<size_t>("RANDOM_ALGO_DUMMY_FACTOR");
        auto dim = _pbParams->getAttributeValue<size_t>("DIMENSION");
        if (dummyFactor < NOMAD::INF_SIZE_T)
        {
            NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( dim*dummyFactor ); // In this example, the single pass (lap) max bb eval is set.
        }
        
        // The algorithm has its own stop reasons.
        // TEMPLATE for a new search method: adapt for the new Algo.
        _randomAlgoStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::RandomAlgoStopType>>();
        
        // Create the algorithm with its own stop reason
        // TEMPLATE use for a new search method: adapt for the new Algo.
        _randomAlgo = std::make_unique<NOMAD::TemplateAlgo>(this,
                                              _randomAlgoStopReasons ,
                                              _runParams,
                                              _pbParams);
        
    }
}


bool NOMAD::TemplateAlgoSearchMethod::runImp()
{
   
    _randomAlgo->setEndDisplay(false);

    _randomAlgo->start();
    bool foundBetter = _randomAlgo->run();
    _randomAlgo->end();

    // Maybe use _randomAlgoStopReason to update parent algorithm
    
    return foundBetter;
}


void NOMAD::TemplateAlgoSearchMethod::generateTrialPointsFinal()
{
    // The trial points of one iteration of this template algorithm (random) are generated (not evaluated).


    // Use first point of barrier as simplex center.
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
