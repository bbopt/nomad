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


#include "../../Algos/CatMads/CatMads.hpp"
#include "../../Algos/CatAds/CatAds.hpp"
#include "../../Algos/CatMads/CatSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"

void NOMAD::CatSearchMethod::init()
{
    // Only run for CatMads or CatAds, not Mads
    
    setStepType(NOMAD::StepType::SEARCH_METHOD_CATALGO_CAT);
    setEnabled(true);
}


void NOMAD::CatSearchMethod::generateTrialPointsFinal()
{
    
    // Use the feasible and incumbent feasible to obtain new trial points
    auto barrier = getMegaIterationBarrier();
    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatSearchMethod: must have a MadsMegaIteration ancestor with a barrier");
    }
    auto feasFrameCenter = barrier->getCurrentIncumbentFeas();
    auto infeasFrameCenter = barrier->getCurrentIncumbentInf();
    
    if (feasFrameCenter == nullptr && infeasFrameCenter == nullptr)
    {
        return;
    }
    
    // CatMads or CatAds algorithm is holding the wrapper function to access the python implementation
    auto catAlgo = dynamic_cast<const NOMAD::CatAlgoUtils*>(getRootAlgorithm());
    if (!catAlgo)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatSearchMethod: must have a CatMads or a CatAds ancestor");
    }
    std::vector<NOMAD::EvalPoint> trialPoints = catAlgo->runCatCategoricalSearch(feasFrameCenter, infeasFrameCenter);
    
    OUTPUT_INFO_START
    AddOutputInfo("CatSearchMethod: generate trial points (done): "+ std::to_string(trialPoints.size()));
    OUTPUT_INFO_END
        
    // -- insert points for processing and evaluation-- //
    for (auto& tp : trialPoints)
    {
        if (nullptr != feasFrameCenter)
        {
            tp.setPointFrom(feasFrameCenter, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        }
        else
        {
            tp.setPointFrom(infeasFrameCenter, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        }
        tp.addGenStep(getStepType());
        OUTPUT_INFO_START
        std::string s = tp.display();
        AddOutputInfo("CatSearchMethod: trial point = " + s);
        OUTPUT_INFO_END
        insertTrialPoint(tp);
        
    }

}
