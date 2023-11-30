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

#include "../../Algos/Mads/SimpleLineSearchMethod.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/SimpleLineSearch/SimpleLineSearch.hpp"


void NOMAD::SimpleLineSearchMethod::init()
{
    // For some testing, it is possible that _runParams is null or evaluator control is null

    if ( nullptr != _runParams && nullptr != NOMAD::EvcInterface::getEvaluatorControl() )
    {
        setStepType(NOMAD::StepType::SEARCH_METHOD_SIMPLE_LINE_SEARCH);

        bool enabled = _runParams->getAttributeValue<bool>("SIMPLE_LINE_SEARCH");
        
        if (enabled && _runParams->getAttributeValue<bool>("SPECULATIVE_SEARCH"))
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"SimpleLineSearchMethod: cannot work with speculative search.");
        }
        
        setEnabled(enabled);

        
        _simpleLineSearchStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::SimpleLineSearchStopType>>();
        
        _simpleLineSearch = std::make_unique<NOMAD::SimpleLineSearch>(this,
                                            _simpleLineSearchStopReasons ,
                                            _runParams,
                                            _pbParams);
        
    }
}


bool NOMAD::SimpleLineSearchMethod::runImp()
{
   
    _simpleLineSearch->setEndDisplay(false);

    _simpleLineSearch->start();
    bool foundBetter = _simpleLineSearch->run();
    _simpleLineSearch->end();

    // Maybe use _simpleLineSearchStopReason to update parent algorithm
    
    return foundBetter;
}


void NOMAD::SimpleLineSearchMethod::generateTrialPointsFinal()
{
    throw NOMAD::Exception(__FILE__,__LINE__,"SimpleLineSearchMethod: cannot work with MegaSearchPoll.");
}
