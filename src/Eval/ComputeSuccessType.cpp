/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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
/**
 \file   ComputeSuccessType.cpp
 \brief  Comparison methods for EvalPoints (implementation)
 \author Viviane Rochon Montplaisir
 \date   March 2017 / September 2020
 \see    ComputeSuccessType.hpp
 */
#include "../Eval/ComputeSuccessType.hpp"

/*--------------------------*/
/* Class ComputeSuccessType */
/*--------------------------*/
NOMAD::SuccessType NOMAD::ComputeSuccessType::defaultComputeSuccessType(
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint1,
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint2,
                                const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    if (nullptr != evalPoint1)
    {
        if (nullptr == evalPoint2)
        {
            if (evalPoint1->getH(NOMAD::EvalType::BB) > hMax)
            {
                // Even if evalPoint2 is NULL, this case is still
                // not a success.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else
            {
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
        }
        else
        {
            success = NOMAD::Eval::defaultComputeSuccessType(evalPoint1->getEval(NOMAD::EvalType::BB),
                                                             evalPoint2->getEval(NOMAD::EvalType::BB),
                                                             hMax);
        }
    }

    return success;
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::computeSuccessTypePhaseOne(
                            const std::shared_ptr<NOMAD::EvalPoint>& evalPoint,
                            const std::shared_ptr<NOMAD::EvalPoint>& xInf,
                            const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    if (nullptr != evalPoint)
    {
        if (nullptr == xInf)
        {
            success = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else
        {
            success = NOMAD::Eval::computeSuccessTypePhaseOne(evalPoint->getEval(NOMAD::EvalType::BB),
                                                              xInf->getEval(NOMAD::EvalType::BB), hMax);
        }
    }

    return success;
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::computeSuccessTypeSgte(
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint1,
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint2,
                                const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;
    const NOMAD::EvalType evalTypeSgte = NOMAD::EvalType::SGTE;

    if (nullptr != evalPoint1)
    {
        if (evalPoint1->getH(evalTypeSgte) > hMax)
        {
            success = NOMAD::SuccessType::UNSUCCESSFUL;
        }
        else if (nullptr == evalPoint2)
        {
            success = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else
        {
            success = NOMAD::Eval::defaultComputeSuccessType(evalPoint1->getEval(evalTypeSgte),
                                                             evalPoint2->getEval(evalTypeSgte),
                                                             hMax);
        }
    }

    return success;
}


void NOMAD::ComputeSuccessType::setDefaultComputeSuccessTypeFunction(const NOMAD::EvalType& evalType)
{
    switch (evalType)
    {
        case NOMAD::EvalType::BB:
            setComputeSuccessTypeFunction(NOMAD::ComputeSuccessType::defaultComputeSuccessType);
            break;
        case NOMAD::EvalType::SGTE:
            setComputeSuccessTypeFunction(NOMAD::ComputeSuccessType::computeSuccessTypeSgte);
            break;
        case NOMAD::EvalType::UNDEFINED:
        default:
            break;
    }
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::operator()(const NOMAD::EvalPointPtr& p1,
                                                         const NOMAD::EvalPointPtr& p2,
                                                         const NOMAD::Double& hMax)
{
    return _computeSuccessType(p1, p2, hMax);
}


