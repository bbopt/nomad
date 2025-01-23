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
#include "../../Algos/Mads/UserSearchMethod.hpp"

#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/SubproblemManager.hpp"

void NOMAD::UserSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_USER);

    // Query the enabling parameter here
    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"UserSearchMethod (" + std::to_string(_id) +"): must have an iteration ancestor");
    }

    if (_id > 2 || _id <= 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"UserSearchMethod (" + std::to_string(_id) +"): only id =1 and id =2 are supported");
    }
    auto mads = dynamic_cast<const NOMAD::Mads*>(_iterAncestor->getRootAlgorithm());
    if ( nullptr != mads)
    {
        setEnabled(mads->hasUserSearchMethod());
    }
    else
    {
        setEnabled(false);
    }
}


void NOMAD::UserSearchMethod::generateTrialPointsFinal()
{

    // The frame center is only used to compute bounds, if they are not defined.
    // Use the first available point.
    auto barrier = getMegaIterationBarrier();
    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"UserSearchMethod (" + std::to_string(_id) +"): must have a MadsMegaIteration ancestor with a barrier");
    }
    auto frameCenter = barrier->getFirstPoint();

    OUTPUT_INFO_START
    AddOutputInfo("Generate point for " + getName() + "(" + std::to_string(_id) +")");
    OUTPUT_INFO_END

    auto mads = dynamic_cast<const NOMAD::Mads*>(_iterAncestor->getRootAlgorithm());


    bool success =false;
    if (_id == 1)
    {
        success = mads->runCallback(NOMAD::CallbackType::USER_METHOD_SEARCH, *this, _trialPoints);
    }
    else if (_id == 2)
    {
        success = mads->runCallback(NOMAD::CallbackType::USER_METHOD_SEARCH_2, *this, _trialPoints);
    }

    if (!success)
    {
        OUTPUT_INFO_START
        AddOutputInfo("User search (" + std::to_string(_id) +") cannot produce directions.");
        OUTPUT_INFO_END
        return;
    }

    // Insert the point. Projection on mesh and snap to bounds is done later
    for (auto tp : _trialPoints)
    {
        tp.setPointFrom(std::make_shared<NOMAD::EvalPoint>(*frameCenter), NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this)); // !!! Point from is a copy of frame center
        tp.addGenStep(getStepType());
        insertTrialPoint(tp);
    }
}

void NOMAD::UserSearchMethod::updateAtStepEnd()
{
    auto mads = dynamic_cast<const NOMAD::Mads*>(_iterAncestor->getRootAlgorithm());

    bool success = mads->runCallback(NOMAD::CallbackType::USER_METHOD_SEARCH_END, *this);

    if (!success)
    {
        OUTPUT_INFO_START
        AddOutputInfo("User search (" + std::to_string(_id) +") post evaluation is not working properly.");
        OUTPUT_INFO_END
        return;
    }
}
