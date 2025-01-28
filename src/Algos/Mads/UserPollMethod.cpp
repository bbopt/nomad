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
#include "../../Algos/Mads/UserPollMethod.hpp"

#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/SubproblemManager.hpp"

void NOMAD::UserPollMethod::init()
{
    setStepType(NOMAD::StepType::POLL_METHOD_USER);
    verifyParentNotNull();

    // Query the enabling parameter here
    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"UserPollMethod: must have an iteration ancestor");
    }
    auto mads = dynamic_cast<const NOMAD::Mads*>(_iterAncestor->getRootAlgorithm());
    if ( nullptr == mads || (!mads->hasUserPollMethod() && !mads->hasUserFreePollMethod()))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"UserPollMethod: the custom callback function for the user poll method must be added. This works only in library mode. See example in $NOMAD_HOME/examples/advanced/library/CustomPollMethod");
    }

}


// Generate poll directions
void NOMAD::UserPollMethod::generateUnitPollDirections(std::list<NOMAD::Direction> &directions, const size_t n) const
{

    directions.clear();

    auto mads = dynamic_cast<const NOMAD::Mads*>(_iterAncestor->getRootAlgorithm());

    // NOTE: cannot have both USER_METHOD_POLL and USER_METHOD_FREE_POLL callbacks provided.
    bool success;
    if (isFreePoll())
    {
        success = mads->runCallback(NOMAD::CallbackType::USER_METHOD_FREE_POLL, *this, directions, n);
    }
    else
    {
        success = mads->runCallback(NOMAD::CallbackType::USER_METHOD_POLL, *this, directions, n);
    }

    if (!success || directions.empty())
    {
        OUTPUT_INFO_START
        AddOutputInfo("User-defined poll method did not produced directions.");
        OUTPUT_INFO_END
        return;
    }

    // Free user poll method allows to provide any directions.
    // Regular user poll method must follow the same conditions as other poll method.
    if (!isFreePoll())
    {
        bool shownWarningMessage = false;
        for (const auto & dir: directions)
        {
            if (! shownWarningMessage && dir.squaredL2Norm() != 1)
            {
                OUTPUT_INFO_START
                AddOutputInfo("WARNING: User-defined poll method produces directions of L2 norm not equal to one. For proper scaling on the mesh in Mads, the directions should have norm 1.");
                OUTPUT_INFO_END
                shownWarningMessage = true;
            }
        }
    }
}


void NOMAD::UserPollMethod::updateEndUserPoll()
{
    auto mads = dynamic_cast<const NOMAD::Mads*>(_iterAncestor->getRootAlgorithm());

    bool success = mads->runCallback(NOMAD::CallbackType::USER_METHOD_FREE_POLL_END, *this);

    if (!success)
    {
        OUTPUT_INFO_START
        AddOutputInfo("User poll post evaluation function is not working properly.");
        OUTPUT_INFO_END
        return;
    }
}
