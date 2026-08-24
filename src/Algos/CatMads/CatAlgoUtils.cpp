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


#include "../../Algos/CatMads/CatAlgoUtils.hpp"

bool NOMAD::CatAlgoUtils::setBuildModelCallback(std::function<bool()> buildCatModel)
{
    _buildCatModel = buildCatModel;
    return true;
}


bool NOMAD::CatAlgoUtils::buildCatModel() const
{
    return _buildCatModel();
}

// For categorical free poll
bool NOMAD::CatAlgoUtils::setPollFreeCallback(std::function<std::vector<NOMAD::Direction>(NOMAD::EvalPointPtr)> catPollFree)
{
    _catPollFree = catPollFree;
    return true;
}

std::vector<NOMAD::Direction> NOMAD::CatAlgoUtils::runCatPollFree(EvalPointPtr frameCenter) const
{
    return _catPollFree(frameCenter);
}

// For categorical search
bool NOMAD::CatAlgoUtils::setCategoricalSearchCallback(std::function<std::vector<NOMAD::EvalPoint>(NOMAD::EvalPointPtr, NOMAD::EvalPointPtr)> catCategoricalSearch)
{
    _catCategoricalSearch = catCategoricalSearch;
    return true;
}

std::vector<NOMAD::EvalPoint> NOMAD::CatAlgoUtils::runCatCategoricalSearch(NOMAD::EvalPointPtr feasFrameCenter, NOMAD::EvalPointPtr infeasFrameCenter) const
{
    return _catCategoricalSearch(feasFrameCenter, infeasFrameCenter);
}
