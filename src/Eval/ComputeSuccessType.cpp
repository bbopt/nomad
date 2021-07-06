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
            if (evalPoint1->getH(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD) > hMax)
            {
                // Even if evalPoint2 is NULL, this case is still
                // not a success.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else if (evalPoint1->isFeasible(NOMAD::EvalType::BB))
            {
                // New feasible point: full success
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
            else
            {
                // New infeasible makes for partial success, not full success
                success = NOMAD::SuccessType::PARTIAL_SUCCESS;
            }
        }
        else
        {
            success = NOMAD::Eval::computeSuccessType(evalPoint1->getEval(NOMAD::EvalType::BB),
                                                      evalPoint2->getEval(NOMAD::EvalType::BB),
                                                      NOMAD::ComputeType::STANDARD,
                                                      hMax);
        }
    }

    return success;
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::computeSuccessTypePhaseOne(
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint1,
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint2,
                                const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    if (nullptr != evalPoint1)
    {
        if (nullptr == evalPoint2)
        {
            success = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else
        {
            success = NOMAD::Eval::computeSuccessType(evalPoint1->getEval(NOMAD::EvalType::BB),
                                                      evalPoint2->getEval(NOMAD::EvalType::BB),
                                                      NOMAD::ComputeType::PHASE_ONE,
                                                      hMax);
        }
    }

    return success;
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::computeSuccessTypePhaseOneSurrogate(
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint1,
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint2,
                                const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    if (nullptr != evalPoint1)
    {
        if (nullptr == evalPoint2)
        {
            success = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else
        {
            success = NOMAD::Eval::computeSuccessType(evalPoint1->getEval(NOMAD::EvalType::SURROGATE),
                                                      evalPoint2->getEval(NOMAD::EvalType::SURROGATE),
                                                      NOMAD::ComputeType::PHASE_ONE,
                                                      hMax);
        }
    }

    return success;
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::computeSuccessTypeModel(
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint1,
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint2,
                                const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;
    const NOMAD::EvalType evalTypeModel = NOMAD::EvalType::MODEL;

    if (nullptr != evalPoint1)
    {
        if (evalPoint1->getH(evalTypeModel, NOMAD::ComputeType::STANDARD) > hMax)
        {
            success = NOMAD::SuccessType::UNSUCCESSFUL;
        }
        else if (nullptr == evalPoint2)
        {
            success = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else
        {
            success = NOMAD::Eval::computeSuccessType(evalPoint1->getEval(evalTypeModel),
                                                      evalPoint2->getEval(evalTypeModel),
                                                      NOMAD::ComputeType::STANDARD,
                                                      hMax);
        }
    }

    return success;
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::computeSuccessTypeSurrogate(
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint1,
                                const std::shared_ptr<NOMAD::EvalPoint>& evalPoint2,
                                const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    if (nullptr != evalPoint1)
    {
        if (nullptr == evalPoint2)
        {
            if (evalPoint1->getH(NOMAD::EvalType::SURROGATE, NOMAD::ComputeType::STANDARD) > hMax)
            {
                // Even if evalPoint2 is NULL, this case is still
                // not a success.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else if (evalPoint1->isFeasible(NOMAD::EvalType::SURROGATE))
            {
                // New feasible point: full success
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
            else
            {
                // New infeasible makes for partial success, not full success
                success = NOMAD::SuccessType::PARTIAL_SUCCESS;
            }
        }
        else
        {
            success = NOMAD::Eval::computeSuccessType(evalPoint1->getEval(NOMAD::EvalType::SURROGATE),
                                                      evalPoint2->getEval(NOMAD::EvalType::SURROGATE),
                                                      NOMAD::ComputeType::STANDARD,
                                                      hMax);
        }
    }

    return success;
}


void NOMAD::ComputeSuccessType::setComputeSuccessTypeFunction(const NOMAD::EvalType& evalType,
                                                              const NOMAD::ComputeType& computeType)
{
    if (NOMAD::EvalType::BB == evalType)
    {
        if (NOMAD::ComputeType::STANDARD == computeType)
        {
            _computeSuccessType = NOMAD::ComputeSuccessType::defaultComputeSuccessType;
        }
        else if (NOMAD::ComputeType::PHASE_ONE == computeType)
        {
            _computeSuccessType = NOMAD::ComputeSuccessType::computeSuccessTypePhaseOne;
        }
        else
        {

        }
    }
    else if (NOMAD::EvalType::SURROGATE == evalType)
    {
        if (NOMAD::ComputeType::STANDARD == computeType)
        {
            _computeSuccessType = NOMAD::ComputeSuccessType::computeSuccessTypeSurrogate;
        }
        else if (NOMAD::ComputeType::PHASE_ONE == computeType)
        {
            _computeSuccessType = NOMAD::ComputeSuccessType::computeSuccessTypePhaseOneSurrogate;
        }
        else
        {
        }
    }
    else if (NOMAD::EvalType::MODEL == evalType)
    {
        _computeSuccessType = NOMAD::ComputeSuccessType::computeSuccessTypeModel;
    }
    else
    {
        std::string err = "No compute success type function available for " + NOMAD::evalTypeToString(evalType);
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }
}


NOMAD::SuccessType NOMAD::ComputeSuccessType::operator()(const NOMAD::EvalPointPtr& p1,
                                                         const NOMAD::EvalPointPtr& p2,
                                                         const NOMAD::Double& hMax)
{
    return _computeSuccessType(p1, p2, hMax);
}


