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


