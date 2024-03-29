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

#include <algorithm>    // For std::merge and std::unique

#include "../Algos/Iteration.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Type/CallbackType.hpp"


NOMAD::Iteration::~Iteration()
{
    NOMAD::OutputQueue::Flush();
}


void NOMAD::Iteration::init()
{
    setStepType(NOMAD::StepType::ITERATION);
    verifyParentNotNull();
    
    _userCallbackEnabled = false;
    if (nullptr != _runParams)
    {
        _userCallbackEnabled = _runParams->getAttributeValue<bool>("USER_CALLS_ENABLED");
    }
}


std::string NOMAD::Iteration::getName() const
{
    return getAlgoName() + NOMAD::stepTypeToString(_stepType) + " " + std::to_string(_k);
}


void NOMAD::Iteration::endImp()
{
    OUTPUT_INFO_START
    AddOutputInfo("Stop reason: " + _stopReasons->getStopReasonAsString() );
    OUTPUT_INFO_END
    if ( _userCallbackEnabled )
    {
        bool stop = false;

        // Callback user provided function to check if user requested a stop.
        runCallback(NOMAD::CallbackType::ITERATION_END, *this, stop);
        if (!_stopReasons->checkTerminate() && stop)
        {
            _stopReasons->set(NOMAD::BaseStopType::USER_GLOBAL_STOP);
        }
        
        // Reset user iteration stop reason
        if (_stopReasons->testIf(NOMAD::IterStopType::USER_ITER_STOP))
        {
            _stopReasons->set(NOMAD::IterStopType::STARTED);
        }
    }
}




