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


#include "../../Algos/CatMads/CatSpeculativeSearchMethod.hpp"
#include "../../Algos/Mads/SpeculativeSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"

void NOMAD::CatSpeculativeSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_CAT_SPECULATIVE);
    
    bool enabled = false;
    // For some testing, it is possible that _runParams is null
    if (nullptr != _runParams)
    {
        enabled = _runParams->getAttributeValue<bool>("SPECULATIVE_SEARCH");
    }
    setEnabled(enabled);

}


void NOMAD::CatSpeculativeSearchMethod::generateTrialPointsFinal()
{
    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatSpeculativeSearchMethod: must have an iteration ancestor");
    }
    
     // The frame center is only used to compute bounds, if they are not defined.
     // Use the first available point.
     auto barrier = getMegaIterationBarrier();
     if (nullptr == barrier)
     {
         throw NOMAD::Exception(__FILE__,__LINE__,"CatSpeculativeSearchMethod: must have a (M)adsMegaIteration ancestor with a barrier");
     }
    
     std::vector<std::shared_ptr<NOMAD::EvalPoint>> frameCenters;
     auto firstXIncFeas = barrier->getCurrentIncumbentFeas();
     auto firstXIncInf  = barrier->getCurrentIncumbentInf();
     if (firstXIncFeas)
     {
         frameCenters.push_back(firstXIncFeas);
     }
     if (firstXIncInf)
     {
         frameCenters.push_back(firstXIncInf);
     }
    
    
     for (const auto & frameCenter : frameCenters)
     {

         // Speculative search only on quantitative success
         auto genStep = frameCenter->getGenStep();
        
         if (genStep == NOMAD::StepType::POLL_METHOD_CAT || genStep == NOMAD::StepType::SEARCH_METHOD_CATALGO_CAT)
         {
             OUTPUT_INFO_START
             AddOutputInfo("This frame center is from a cat success. Skip speculative search.");
             OUTPUT_INFO_END
             continue;
         }
         else
         {
             OUTPUT_INFO_START
             std::string s = "Speculative search with frame center " + frameCenter->display(NOMAD::defaultFHComputeTypeS) + " obtained by " + NOMAD::stepTypeToString(genStep);
             AddOutputInfo(s);
             OUTPUT_INFO_END
             generateTrialPointsFromFC(frameCenter);
         }

     }
}
