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

#include <algorithm>    // For std::merge and std::unique

#include "../../Algos/CoordinateSearch/CSIteration.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::CSIteration::init()
{
    setStepType(NOMAD::StepType::ITERATION);
    
    if (nullptr != _runParams && _runParams->getAttributeValue<bool>("MEGA_SEARCH_POLL"))
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"CS does not support Mega Search Poll. ");
    }
    else
    {
        _csPoll = std::make_unique<NOMAD::CSPoll>(this);
    }
    
}


void NOMAD::CSIteration::startImp()
{
}


bool NOMAD::CSIteration::runImp()
{
    bool iterationSuccess = false;
    
    OUTPUT_INFO_START
    std::string s = "No search method called by Coordinate Search algorithm.";
    AddOutputDebug(s);
    OUTPUT_INFO_END
    
    if ( ! _stopReasons->checkTerminate() )
    {
        // 2. CS Poll
        _csPoll->start();
        // Iteration is a success if either a better xFeas or
        // a better xInf (partial success or dominating) xInf was found.
        // See Algorithm 12.2 from DFBO.
        iterationSuccess = _csPoll->run();
        _csPoll->end();
        
    }
    
    // End of the iteration: iterationSuccess is true iff we have a full success.
    return iterationSuccess;
}



