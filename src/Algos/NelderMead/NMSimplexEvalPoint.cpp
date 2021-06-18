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

#include "../../Algos/EvcInterface.hpp" // To access EvalType
#include "../../Algos/NelderMead/NMSimplexEvalPoint.hpp"
#include "../../Eval/ComputeSuccessType.hpp"

bool NOMAD::NMSimplexEvalPointCompare::operator()(const NOMAD::EvalPoint& lhs,
                                                  const NOMAD::EvalPoint& rhs) const
{
    // Workaround to get EvalType for ComputeSuccessType.
    NOMAD::EvalType evalType = NOMAD::EvalType::BB;
    NOMAD::ComputeType computeType = NOMAD::ComputeType::STANDARD;
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (nullptr != evc)
    {
        evalType = evc->getEvalType();
        computeType = evc->getComputeType();
    }
    NOMAD::ComputeSuccessType computeSuccess(evalType, computeType);

    // This could be an inefficient function (see issue #392)
    NOMAD::SuccessType success = computeSuccess(std::make_shared<NOMAD::EvalPoint>(lhs),
                                                std::make_shared<NOMAD::EvalPoint>(rhs));

    if (success >= NOMAD::SuccessType::FULL_SUCCESS)
    {
        return true;
    }

    success = computeSuccess(std::make_shared<NOMAD::EvalPoint>(rhs),
                             std::make_shared<NOMAD::EvalPoint>(lhs));
    if (success >= NOMAD::SuccessType::FULL_SUCCESS)
    {
        return false;
    }

    // No dominance, compare h values.
    NOMAD::Double h1 = lhs.getH(evalType, computeType);
    NOMAD::Double h2 = rhs.getH(evalType, computeType);
    if (h1.isDefined() && h2.isDefined())
    {
        if (h1 < h2)
        {
            return true;
        }
        if (h2 < h1)
        {
            return false;
        }
    }
    else if (h1.isDefined() && !h2.isDefined())
    {
        return true;
    }
    else if (!h1.isDefined() && h2.isDefined())
    {
        return false;
    }

    // The "older" function from NM-Mads paper is implemented by comparing the tags
    return lhs.getTag() < rhs.getTag();
}

