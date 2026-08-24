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

#include "../../Algos/CatMads/CatPollMethod.hpp"
#include "../../Algos/CatMads/CatMads.hpp"

void NOMAD::CatPollMethod::init()
{
    setStepType(NOMAD::StepType::POLL_METHOD_CAT);
    verifyParentNotNull();

    // Query the enabling parameter here
    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatPollMethod: must have an iteration ancestor");
    }
    auto catAlgo = dynamic_cast<const NOMAD::CatAlgoUtils*>(_iterAncestor->getRootAlgorithm());
    if ( nullptr == catAlgo)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatPollMethod: Cat poll only works with CatMads/CatAds.");
    }
}


// Generate poll directions
void NOMAD::CatPollMethod::generateUnitPollDirections(std::list<NOMAD::Direction> &directions, const size_t n) const
{
    directions.clear();

    // Feasible or infeasible solution might not exist: the poll is skipped
    auto frameCenterTmp = getFrameCenter();
    if (nullptr == frameCenterTmp)
    {
            OUTPUT_INFO_START
            AddOutputInfo("Cat poll method: no frame center, skipping poll.");
            OUTPUT_INFO_END
        return;
    }
    
    auto catAlgo = dynamic_cast<const NOMAD::CatAlgoUtils*>(_iterAncestor->getRootAlgorithm());
    
    std::vector<NOMAD::Direction> catPollDirections = catAlgo->runCatPollFree(frameCenterTmp);
    
    size_t dim = frameCenterTmp->size();
    NOMAD::Direction dir(dim, 0.0); // Initialize a direction vector of size n
    
    // -- Store directions -- //
    for (const auto& pd : catPollDirections)
    {
        if (pd.size() != dim)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"CatPoll direction not consistent with the subproblem dimension n");
        }
        OUTPUT_INFO_START
        std::string s = pd.display();
        AddOutputInfo("Cat poll method direction:"+s);
        OUTPUT_INFO_END
        directions.push_back(pd); // Store the constructed direction
    }
    
    if (directions.empty())
    {
        OUTPUT_INFO_START
        AddOutputInfo("Cat poll method did not produced directions.");
        OUTPUT_INFO_END
        return;
    }

    

}
