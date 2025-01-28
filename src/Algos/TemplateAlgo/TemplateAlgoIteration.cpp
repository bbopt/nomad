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

#include "../../nomad_platform.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoIteration.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoMegaIteration.hpp"

void NOMAD::TemplateAlgoIteration::init()
{
    setStepType(NOMAD::StepType::ITERATION);
    
    _templateAlgoRandom = std::make_unique<NOMAD::TemplateAlgoRandom>(this);
    _templateAlgoUpdate = std::make_unique<NOMAD::TemplateAlgoUpdate> (this);

}

void NOMAD::TemplateAlgoIteration::startImp()
{
    // For illustration purpose the Update is used to update the center point
    // (the best feasible or best infeasible) around which the trial points are generated.
    _templateAlgoUpdate->start();
    bool updateSuccess = _templateAlgoUpdate->run();
    _templateAlgoUpdate->end();
    
    if ( ! updateSuccess )
    {
        auto stopReason = NOMAD::AlgoStopReasons<NOMAD::RandomAlgoStopType>::get ( getAllStopReasons() );

        // The update is not a success. If the global stop reason is not set to terminate we set a default stop reason for initialization.
        if ( !_stopReasons->checkTerminate() )
            stopReason->set( NOMAD::RandomAlgoStopType::UPDATE_FAILED);
    }
}


bool NOMAD::TemplateAlgoIteration::runImp()
{
    // Iteration cannot generate all points before evaluation
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    bool iterationSuccess = false;
    
    
    _templateAlgoRandom->start();
    
    // Iteration is a success if either a better xFeas or
    // a better xInf (partial success or dominating) xInf was found.
    iterationSuccess = _templateAlgoRandom->run();

    _templateAlgoRandom->end();
    
    if ( iterationSuccess )
    {
        // Update the MegaIteration best success type with success found.
        getParentOfType<NOMAD::MegaIteration*>()->setSuccessType(_templateAlgoRandom->getSuccessType());
    }

    // End of the iteration: iterationSuccess is true if we have a success.
    return iterationSuccess;

}
