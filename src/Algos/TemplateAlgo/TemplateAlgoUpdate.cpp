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

#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoMegaIteration.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoUpdate.hpp"

void NOMAD::TemplateAlgoUpdate::init()
{
    setStepType(NOMAD::StepType::UPDATE);
    verifyParentNotNull();
}


std::string NOMAD::TemplateAlgoUpdate::getName() const
{
    return getAlgoName() + NOMAD::stepTypeToString(_stepType);
}

bool NOMAD::TemplateAlgoUpdate::runImp()
{
    
    // Update the Algo Iteration best point from MegaIteration barrier
    // THIS IS TO ILLUSTRATE AN UPDATE. It could have been done more direct.
    
    bool updateSuccess = false;
    auto barrier = getMegaIterationBarrier();
    auto iter = getParentOfType<NOMAD::TemplateAlgoIteration*>();
    
    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Update must have a barrier in the MegaIteration among its ancestors.");
    }
    if (nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Update must have an Iteration among its ancestors.");
    }
    auto bestXFeas = barrier->getFirstXFeas();
    auto bestXInf  = barrier->getFirstXInf();
    
    // This is the bestXFeas or bestXInf from MegaIteration barrier.
    // Update them (done at the start of Iteration)
    if (nullptr != bestXFeas)
    {
        iter->setFrameCenter(bestXFeas);
        updateSuccess = true;
    }
    else if (nullptr != bestXInf)
    {
        iter->setFrameCenter(bestXInf);
        updateSuccess = true;
    }
    OUTPUT_DEBUG_START
    auto frameCenter = iter->getFrameCenter();
    AddOutputDebug("Current frame center: " + (frameCenter ? frameCenter->display() : "NULL"));
    OUTPUT_DEBUG_END

    return updateSuccess;
}
