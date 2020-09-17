
#include "../../Algos/EvcInterface.hpp" // To access EvalType
#include "../../Algos/NelderMead/NMSimplexEvalPoint.hpp"
#include "../../Eval/ComputeSuccessType.hpp"

bool NOMAD::NMSimplexEvalPointCompare::operator()(const NOMAD::EvalPoint& lhs,
                                                  const NOMAD::EvalPoint& rhs) const
{
    // Workaround to get EvalType for ComputeSuccessType.
    NOMAD::EvalType evalType = NOMAD::EvalType::BB;
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (nullptr != evc)
    {
        evalType = evc->getEvalType();
    }
    NOMAD::ComputeSuccessType computeSuccess(evalType);

    NOMAD::SuccessType success = computeSuccess(std::make_shared<NOMAD::EvalPoint>(lhs),
                                                std::make_shared<NOMAD::EvalPoint>(rhs));

    if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
    {
        return true;
    }

    success = computeSuccess(std::make_shared<NOMAD::EvalPoint>(rhs),
                             std::make_shared<NOMAD::EvalPoint>(lhs));
    if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
    {
        return false;
    }

    // The "older" function from NM-Mads paper is implemented by comparing the tags
    return lhs.getTag() < rhs.getTag();
}

