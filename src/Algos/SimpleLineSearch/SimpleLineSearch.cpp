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

#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Util/fileutils.hpp"

#include "../../Algos/Mads/SimpleLineSearchMethod.hpp"
#include "../../Algos/SimpleLineSearch/SimpleLineSearch.hpp"
#include "../../Algos/SimpleLineSearch/SimpleLineSearchMegaIteration.hpp"

void NOMAD::SimpleLineSearch::init()
{

    setStepType(NOMAD::StepType::SEARCH_METHOD_SIMPLE_LINE_SEARCH);
    verifyParentNotNull();

    const auto parentSearch = getParentOfType<NOMAD::SimpleLineSearchMethod*>(false);

    if (nullptr == parentSearch)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"SimpleLineSearch: cannot be a standalone method. Must be part of a mads search.");
    }
    
    
}

bool NOMAD::SimpleLineSearch::runImp()
{
    _algoSuccessful = false;
    
    
    if ( ! _stopReasons->checkTerminate() )
    {
        auto barrier = getParentOfType<NOMAD::MegaIteration*>()->getBarrier();
        
        if (nullptr == barrier)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"SimpleLineSearch needs a barrier from a Mega Iteration");
        }
        
        // Create a single MegaIteration: manage multiple iterations.
        NOMAD::SimpleLineSearchMegaIteration megaIteration(this, 0, barrier,NOMAD::SuccessType::UNDEFINED);
        
        megaIteration.start();
        bool currentMegaIterSuccess = megaIteration.run();
        megaIteration.end();
        
        _algoSuccessful = _algoSuccessful || currentMegaIterSuccess;
        
        auto algoStopReason = NOMAD::AlgoStopReasons<NOMAD::SimpleLineSearchStopType>::get ( _stopReasons );
        algoStopReason->set( NOMAD::SimpleLineSearchStopType::ALL_POINTS_EVALUATED); // This will stop iterations.
    
        // _refMegaIteration is used to keep values used in Mads::end(). Update it here.
        _refMegaIteration = std::make_shared<NOMAD::SimpleLineSearchMegaIteration>(this, 0, barrier, _success);
        
        _termination->start();
        _termination->run();
        _termination->end();
    }
    
    return _algoSuccessful;
}


void NOMAD::SimpleLineSearch::readInformationForHotRestart()
{
    if (_runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES"))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"SimpleLineSearch: cannot be used with hot restart.");
    }
}

