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
/**
 \file   Eval.cpp
 \brief  Evaluation at a point (implementation)
 \author Viviane Rochon Montplaisir
 \date   March 2017
 \see    Eval.hpp
 */
#include <utility>

#include "../Eval/Eval.hpp"
#include "../Type/EvalType.hpp"


/*---------------------------------------------------------------------*/
/*                            Constructor 1                            */
/*---------------------------------------------------------------------*/
// Note: NOMAD::Double() makes value = NOMAD::NaN; defined = false.
NOMAD::Eval::Eval()
  : _evalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _preEvalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _bbOutput(""),
    _bbOutputTypeList(),
    _bbOutputComplete(false)
{
    _moInfo = std::make_unique<MOInfo>();
}


/*---------------------------------------------------------------------*/
/*                            Constructor 2                            */
/*---------------------------------------------------------------------*/
NOMAD::Eval::Eval(const std::shared_ptr<NOMAD::EvalParameters>& params,
                  NOMAD::BBOutput bbOutput)
  : _evalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _preEvalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _bbOutput(bbOutput),
    _bbOutputTypeList(params->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"))
{
    _bbOutputComplete = _bbOutput.isComplete(_bbOutputTypeList);

    NOMAD::ArrayOfDouble f = _bbOutput.getObjectives(_bbOutputTypeList);
    if (_bbOutput.getEvalOk() && f.isComplete())
    {
        _evalStatus = NOMAD::EvalStatusType::EVAL_OK;
    }
    else
    {
        _evalStatus = NOMAD::EvalStatusType::EVAL_FAILED;
    }
    _moInfo = std::make_unique<MOInfo>();
}


/*---------------------------------------------------------------------*/
/*                           Copy Constructor                          */
/*---------------------------------------------------------------------*/
NOMAD::Eval::Eval(const NOMAD::Eval &eval)
  : _evalStatus(eval._evalStatus),
    _preEvalStatus(eval._preEvalStatus),
    _bbOutput(eval._bbOutput),
    _bbOutputTypeList(eval._bbOutputTypeList),
    _bbOutputComplete(eval._bbOutputComplete)
{
    _moInfo = std::make_unique<NOMAD::MOInfo>(*eval._moInfo);
}

/*-----------------------------------------------------------*/
/*                     Affectation operator                  */
/*-----------------------------------------------------------*/
NOMAD::Eval& NOMAD::Eval::operator=(const NOMAD::Eval& eval)
{
    if (this == &eval)
    {
        return *this;
    }
    _evalStatus = eval._evalStatus;
    _bbOutput = eval._bbOutput;
    _bbOutputTypeList = eval._bbOutputTypeList;
    _bbOutputComplete = eval._bbOutputComplete;

    // Deep copy
    _moInfo = std::make_unique<NOMAD::MOInfo>(*eval._moInfo);

    return *this;
}

/*-----------------------*/
/*     Other methods     */
/*-----------------------*/
bool NOMAD::Eval::isFeasible(const NOMAD::FHComputeTypeS& fhComputeType) const
{
    if (NOMAD::EvalStatusType::EVAL_OK != _evalStatus)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Eval::isFeasible: Needs status type EVAL_OK");
    }
    NOMAD::Double h = getH(fhComputeType);
    return (h.isDefined() && h.todouble() < NOMAD::Double::getEpsilon() + NOMAD::Double::getHMin());
}

// This Eval has a status that permits its point to be re-evaluated.
bool NOMAD::Eval::canBeReEvaluated() const
{
    bool reEval = false;
    if (   _evalStatus == NOMAD::EvalStatusType::EVAL_OK
        || _evalStatus == NOMAD::EvalStatusType::EVAL_NOT_STARTED
        || _preEvalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
        || _evalStatus == NOMAD::EvalStatusType::EVAL_ERROR)
    {
        reEval = true;
    }

    return reEval;

}


// This Eval has evaluation information that is useful to save to cache,
// either by its computed values or by its eval status.
bool NOMAD::Eval::goodForCacheFile() const
{
    bool goodForCache = false;
    if (_evalStatus == NOMAD::EvalStatusType::EVAL_OK
        || _evalStatus == NOMAD::EvalStatusType::EVAL_FAILED
        || _evalStatus == NOMAD::EvalStatusType::EVAL_ERROR
        || _preEvalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
        || _preEvalStatus == NOMAD::EvalStatusType::EVAL_USER_ACCEPTED)
    {
        goodForCache = true;
    }
    return goodForCache;
}


/*------------------------------------*/
/*      Get f. Always recomputed.     */
/*------------------------------------*/
NOMAD::Double NOMAD::Eval::getF(const NOMAD::FHComputeTypeS& fhComputeType) const
{
    NOMAD::Double f;

    if (NOMAD::EvalStatusType::EVAL_OK != _evalStatus)
    {
        return NOMAD::INF;
    }
    switch (fhComputeType.computeType)
    {
        case NOMAD::ComputeType::STANDARD:
            f = _bbOutput.getObjective(_bbOutputTypeList);
            break;
        case NOMAD::ComputeType::DMULTI_COMBINE_F:
            if (_moInfo->fvalues.isEmpty())
            {
                _moInfo->fvalues = _bbOutput.getObjectives(_bbOutputTypeList);
            }
            if (!_moInfo->combineFValue.isDefined())
            {
                _moInfo->combineFValue = fhComputeType.singleObjectiveCompute(_bbOutputTypeList, _bbOutput);
            }
            f = _moInfo->combineFValue;
            break;
        case NOMAD::ComputeType::PHASE_ONE:
            f = computeFPhaseOne(fhComputeType.hNormType);
            break;
        case NOMAD::ComputeType::USER:
            f = fhComputeType.singleObjectiveCompute(_bbOutputTypeList, _bbOutput);
            break;
        default:
            throw NOMAD::Exception(__FILE__,__LINE__,"getF(): ComputeType not supported");
    }

    return f;
}


const NOMAD::ArrayOfDouble& NOMAD::Eval::getFs(const NOMAD::FHComputeTypeS& fhComputeType) const
{
    if (NOMAD::EvalStatusType::EVAL_OK != _evalStatus)
    {
        _moInfo->intermediateVal.resize(1);
        _moInfo->intermediateVal[0] = NOMAD::INF;
        return _moInfo->intermediateVal;
    }
    switch (fhComputeType.computeType)
    {
        case NOMAD::ComputeType::STANDARD:
            if (_moInfo->fvalues.isEmpty())
            {
                _moInfo->fvalues = _bbOutput.getObjectives(_bbOutputTypeList);
            }
            return _moInfo->fvalues;
        case NOMAD::ComputeType::USER:
            if (_moInfo->fvalues.isEmpty())
            {
                _moInfo->fvalues.resize(1);
                _moInfo->fvalues[0] = fhComputeType.singleObjectiveCompute(_bbOutputTypeList, _bbOutput);
            }
            return _moInfo->fvalues;
        case NOMAD::ComputeType::DMULTI_COMBINE_F:
            _moInfo->intermediateVal.resize(1);
            if (_moInfo->fvalues.isEmpty())
            {
                _moInfo->fvalues = _bbOutput.getObjectives(_bbOutputTypeList);
            }
            if (!_moInfo->combineFValue.isDefined())
            {
                _moInfo->combineFValue = fhComputeType.singleObjectiveCompute(_bbOutputTypeList, _bbOutput);
            }
            _moInfo->intermediateVal[0] = _moInfo->combineFValue;
            return _moInfo->intermediateVal;
        case NOMAD::ComputeType::PHASE_ONE:
            _moInfo->intermediateVal.resize(1);
            _moInfo->intermediateVal[0] = computeFPhaseOne(fhComputeType.hNormType);
            return _moInfo->intermediateVal;
        default:
            throw NOMAD::Exception(__FILE__,__LINE__,"getFs(): ComputeType not supported");
    }

}


/*-------------------------------------*/
/*      Get h. Always recomputed.      */
/*-------------------------------------*/
NOMAD::Double NOMAD::Eval::getH(const NOMAD::FHComputeTypeS& fhComputeType) const
{
    NOMAD::Double h;

    if (NOMAD::EvalStatusType::EVAL_OK != _evalStatus)
    {
        return NOMAD::INF;
    }
    switch (fhComputeType.computeType)
    {
        case NOMAD::ComputeType::STANDARD:
        case NOMAD::ComputeType::DMULTI_COMBINE_F:
            h = computeHStandard(fhComputeType.hNormType);
            break;
        case NOMAD::ComputeType::PHASE_ONE:
            h = 0.0;
            break;
        case NOMAD::ComputeType::USER:
            h = fhComputeType.infeasHCompute(_bbOutputTypeList, _bbOutput);
            break;
        default:
            throw NOMAD::Exception(__FILE__,__LINE__,"getH(): ComputeType not supported");
    }

    return h;
}


NOMAD::Double NOMAD::Eval::computeHStandard(NOMAD::HNormType hNormType) const
{
    NOMAD::Double h = 0.0;
    bool hPos = false;

    const NOMAD::ArrayOfDouble bboArray = _bbOutput.getBBOAsArrayOfDouble();
    size_t bboIndex = 0;
    for (const auto & bbOutputType : _bbOutputTypeList)
    {
        const NOMAD::Double& bboI = bboArray[bboIndex];
        bboIndex++;
        if (!bbOutputType.isConstraint())
        {
            continue;
        }
        else if (!bboI.isDefined())
        {
            h = NOMAD::Double();    // h is undefined
            break;
        }
        else if (bboI > 0.0)
        {
            hPos = true;
            NOMAD::Double hTemp = 0.0;
            if (bbOutputType == NOMAD::BBOutputType::Type::EB)
            {
                hTemp = NOMAD::INF;
            }
            else if (bbOutputType == NOMAD::BBOutputType::Type::PB || bbOutputType == NOMAD::BBOutputType::Type::RPB )
            {
                switch (hNormType)
                {
                    case NOMAD::HNormType::L2:
                        hTemp = bboI * bboI;
                        break;
                    case NOMAD::HNormType::L1:
                        hTemp = bboI;
                        break;
                    case NOMAD::HNormType::Linf:
                        if ( bboI > h )
                            h = bboI;
                        break;
                    default:
                        break;
                }

            }

            // Violated Extreme Barrier constraint:
            // Set h to infinity and break.
            if (NOMAD::INF == hTemp)
            {
                h = NOMAD::INF;
                break;
            }
            h += hTemp;
        }
    }

    // Failsafe: If at least one PB constraint is positive, h must be set
    // to at least epsilon so that the Eval is recognized as infeasible.
    // Catch cases such as constraint violated by 1e-8, which gives h = 1e-16
    // which is considered as 0.
    if (hPos && h.isDefined() && (0 == h))
    {
        h = NOMAD::Double::getEpsilon();
    }

    return h;
}


NOMAD::Double NOMAD::Eval::computeFPhaseOne( NOMAD::HNormType hNormType) const
{
    NOMAD::Double f ;
    const NOMAD::ArrayOfDouble bboArray = _bbOutput.getBBOAsArrayOfDouble();
    bool fPos = false;

    if (NOMAD::EvalStatusType::EVAL_OK == _evalStatus)
    {
        f=0.0;
        size_t bboIndex = 0;
        for (const auto & bbOutputType : _bbOutputTypeList)
        {
            const NOMAD::Double& bboI = bboArray[bboIndex];
            bboIndex++;
            if (bbOutputType != NOMAD::BBOutputType::Type::EB)
            {
                continue;
            }
            else if (!bboI.isDefined())
            {
                f = NOMAD::Double();    // Undefined;
            }
            else if (bboI > 0.0)
            {
                fPos = true;
                switch (hNormType)
                {
                    case NOMAD::HNormType::L2:
                        f += bboI * bboI;
                        break;
                    case NOMAD::HNormType::L1:
                        f += bboI;
                        break;
                    case NOMAD::HNormType::Linf:
                        if ( bboI > f )
                            f = bboI;
                        break;
                    default:
                        break;
                }
            }
        }
    }
    else
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"computeFPhaseOne(): EvalStatusType not EVAL_OK");
    }

    // Failsafe: If at least one EB constraint is positive, f is set
    // to at least epsilon.
    if (fPos && f.isDefined() && (0 == f))
    {
        f = NOMAD::Double::getEpsilon();
    }

    return f;
}


/*-----------------------*/
/*      Set BBOutput     */
/*-----------------------*/
void NOMAD::Eval::setBBO(const std::string &bbo,
                         const NOMAD::BBOutputTypeList &bbOutputTypeList,
                         const bool evalOk)
{
    _bbOutput = NOMAD::BBOutput(bbo, evalOk);
    _bbOutputTypeList = bbOutputTypeList;
    _moInfo = std::make_unique<NOMAD::MOInfo>();

    // Revealed constraint are not set by evaluator. They are updated later by a callback.
    // Need to set a default value to pass the following tests.
    updateForRevealedConstraints();

    if (_bbOutputTypeList.empty())
    {
        // Assume it will be set later.
    }
    else if (!_bbOutput.checkSizeMatch(bbOutputTypeList))
    {
        _evalStatus = NOMAD::EvalStatusType::EVAL_ERROR;
        _bbOutputComplete = false;
    }
    else
    {
        _bbOutputComplete = _bbOutput.isComplete(_bbOutputTypeList);
        _evalStatus = _bbOutput.getObjectives(_bbOutputTypeList).isComplete() ? NOMAD::EvalStatusType::EVAL_OK : NOMAD::EvalStatusType::EVAL_FAILED;
    }

}

// Called by setBBO or EvaluatorControl (when setBBO(bbo) is used instead of setBBO(bbo,bbot_list)).
// Currently used only for DiscoMads revealed RPB
void NOMAD::Eval::updateForRevealedConstraints()
{
    if (_bbOutputTypeList.empty())
    {
        return;
    }

    auto bboAOD = _bbOutput.getBBOAsArrayOfDouble();
    // Just ONE RPB constraint can be present
    size_t diffSize = _bbOutputTypeList.size() - bboAOD.size();
    auto it =  std::find(_bbOutputTypeList.begin(),_bbOutputTypeList.end(),NOMAD::BBOutputType::RPB);
    if (diffSize == 1 &&  it!= _bbOutputTypeList.end())
    {
        // Update RPB constraint with a feasible default value.
        _bbOutput = NOMAD::BBOutput(_bbOutput.getBBO()+" -1.0", _bbOutput.getEvalOk());
        _bbOutputComplete = _bbOutput.isComplete(_bbOutputTypeList);
    }


}


/*-----------------------------------------------------------*/
/*                           operator ==                     */
/*-----------------------------------------------------------*/
bool NOMAD::Eval::operator==(const NOMAD::Eval &e) const
{
    // Ignore eval status for comparison.
    // Compare f (single objective) and h.

    bool equal = false;
    NOMAD::Double f1;
    NOMAD::Double f2;
    if (NOMAD::EvalStatusType::EVAL_OK == _evalStatus)
    {
        f1 = getF(NOMAD::defaultFHComputeTypeS);
    }
    if (NOMAD::EvalStatusType::EVAL_OK == e._evalStatus)
    {
        f2 = e.getF(NOMAD::defaultFHComputeTypeS);
    }

    if (this == &e)
    {
        equal = true;
    }
    else if (!f1.isDefined() || !f2.isDefined())
    {
        // If either value is undefined, consider Evals not equal - even if both are undefined.
        equal = false;
    }
    else
    {
        // General case
        NOMAD::Double h1 = getH(NOMAD::defaultFHComputeTypeS);
        NOMAD::Double h2 = e.getH(NOMAD::defaultFHComputeTypeS);
        // As for f, if either h value is undefined, consider Evals not equal - even if both are undefined.
        if (!h1.isDefined() || !h2.isDefined())
        {
            equal = false;
        }
        else
        {
            // Warning. Maybe we need a flag to do f1.todouble() == f2.todouble() && (h1.todouble() == h2.todouble() )
            equal = ( (f1 == f2) && (h1 == h2) );
        }
    }

    return equal;
}

/*--------------------------------------------------------------------------*/
/* Note: This operator must not be used to find/store EvalPoints in the     */
/* cache or in any set.                                                     */
/*--------------------------------------------------------------------------*/
bool NOMAD::Eval::operator<(const NOMAD::Eval &eval) const
{
    return dominates(eval, defaultFHComputeTypeS);
}


/*--------------------------------------------------------------------------*/
/* Dominance as defined by Definition 12.3 of DFBO ("The Book").            */
/* The feasible point x < the feasible point y                              */
/*     when f(x) < f(y).                                                    */
/* The infeasible point x < the infeasible point y                          */
/*     when f(x) <= f(y) and h(x) <= h(y),                                  */
/*     with at least one strict inequality.                                 */
/* Otherwise, return false.                                                 */
/*                                                                          */
/*--------------------------------------------------------------------------*/
bool NOMAD::Eval::dominates(const NOMAD::Eval &eval, const NOMAD::FHComputeTypeS& computeType) const
{

// Original version for testing dominance. Now use the more general compMO.
//    bool dom = false;
//    double f1 = getF(computeType).todouble();
//    NOMAD::Double h1 = getH(computeType);
//    double f2 = eval.getF(computeType).todouble();
//    NOMAD::Double h2 = eval.getH(computeType);
//
//    if (isFeasible(computeType) && eval.isFeasible(computeType))
//    {
//        dom = (f1 < f2);
//    }
//    else if (!isFeasible(computeType) && !eval.isFeasible(computeType))
//    {
//        if (h1 != NOMAD::INF)
//        {
//            dom = (f1 <= f2) && (h1 <= h2) && ((f1 < f2) || (h1 < h2));
//        }
//    }
//    // else - comparing a feasible point with an unfeasible point.
//    // Always false. Do nothing.


    NOMAD::CompareType compare = compMO(eval, computeType, false /* false: compare f and h*/);
    // comparing a feasible point with an unfeasible point --> UNDEFINED -> false.
    return (NOMAD::CompareType::DOMINATING == compare);
}

NOMAD::CompareType NOMAD::Eval::compMO(const NOMAD::Eval &eval,
                                       const NOMAD::FHComputeTypeS& computeType,
                                       const bool onlyfvalues) const
{
    NOMAD::CompareType compareFlag = NOMAD::CompareType::UNDEFINED;

    const NOMAD::ArrayOfDouble& f1 = getFs(computeType);
    const NOMAD::Double h1 = getH(computeType);
    const NOMAD::ArrayOfDouble& f2 = eval.getFs(computeType);
    const NOMAD::Double h2 = eval.getH(computeType);

    // Comparing objective vectors of different size is undefined
    if (f1.size() != f2.size())
    {
        return compareFlag;
    }

    // The comparison code has been adapted from
    // Jaszkiewicz, A., & Lust, T. (2018).
    // ND-tree-based update: a fast algorithm for the dynamic nondominance problem.
    // IEEE Transactions on Evolutionary Computation, 22(5), 778-791.
    if (isFeasible(computeType) && eval.isFeasible(computeType))
    {
        bool isbetter = false;
        bool isworse = false;
        for (size_t i = 0; i < f1.size(); ++i)
        {
            if (f1[i].todouble() < f2[i].todouble())
            {
                isbetter = true;
            }
            if (f2[i].todouble() < f1[i].todouble())
            {
                isworse = true;
            }
            if (isworse && isbetter)
            {
                break;
            }
        }
        if (isworse)
        {
            compareFlag = isbetter ? NOMAD::CompareType::INDIFFERENT : NOMAD::CompareType::DOMINATED;
        }
        else
        {
            compareFlag = isbetter ? NOMAD::CompareType::DOMINATING : NOMAD::CompareType::EQUAL;
        }
    }
    else if (!isFeasible(computeType) && !eval.isFeasible(computeType))
    {
        if (h1 != NOMAD::INF)
        {
            bool isbetter = false;
            bool isworse = false;
            for (size_t i = 0; i < f1.size(); ++i)
            {
                if (f1[i].todouble() < f2[i].todouble())
                {
                    isbetter = true;
                }
                if (f2[i].todouble() < f1[i].todouble())
                {
                    isworse = true;
                }
                if (isworse && isbetter)
                {
                    break;
                }
            }
            if (!(isworse && isbetter) && !onlyfvalues)
            {
                if (h1 < h2)
                {
                    isbetter = true;
                }
                if (h2 < h1)
                {
                    isworse = true;
                }
            }
            if (isworse)
            {
                compareFlag = isbetter ? NOMAD::CompareType::INDIFFERENT : NOMAD::CompareType::DOMINATED;
            }
            else
            {
                compareFlag = isbetter ? NOMAD::CompareType::DOMINATING : NOMAD::CompareType::EQUAL;
            }
        }
    }

    // Comparing an infeasible objective vector with a feasible objective vector is always UNDEFINED.
    return compareFlag;
}


// Comparison function used for Cache's findBest functions.
bool NOMAD::Eval::compEvalFindBest(const NOMAD::Eval &eval1, const NOMAD::Eval &eval2, const NOMAD::FHComputeTypeS& computeType)
{
    // Success is PARTIAL_SUCCESS or FULL_SUCCESS
    // if eval1 is better than eval2.
    // hMax is ignored (set to NOMAD::INF).
    NOMAD::SuccessType success = computeSuccessType(&eval1, &eval2, computeType, NOMAD::INF);

    return (success >= NOMAD::SuccessType::PARTIAL_SUCCESS);
}



NOMAD::SuccessType NOMAD::Eval::computeSuccessType(const NOMAD::Eval* eval1,
                                                   const NOMAD::Eval* eval2,
                                                   const NOMAD::FHComputeTypeS& computeType,
                                                   const NOMAD::Double& hMax)
{
    // NOT_EVALUATED,      // Not evaluated yet
    // UNSUCCESSFUL,       // Failure
    // PARTIAL_SUCCESS,    // Partial success (improving). Found an infeasible
    //                        solution with a better h. f is worse.
    // FULL_SUCCESS        // Full success (dominating)
    NOMAD::SuccessType success = NOMAD::SuccessType::UNDEFINED; /// Will trigger exception if not changed

    if (nullptr != eval1)
    {
        if (nullptr == eval2)
        {
            NOMAD::Double h = eval1->getH(computeType);
            if (!h.isDefined() || h > hMax || h == NOMAD::INF)
            {
                // Even if eval2 is NULL, this case is not successful.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else
            {
                // A new infeasible point, without prior infeasible point, is partial success,
                // not a full success.
                if (eval1->isFeasible(computeType))
                {
                    success = NOMAD::SuccessType::FULL_SUCCESS;
                }
                else
                {
                    success = NOMAD::SuccessType::PARTIAL_SUCCESS;
                }
            }
        }
        else
        {
            if (eval1->getH(computeType) <= hMax
                && eval2->getH(computeType) <= hMax
                && eval1->dominates(*eval2, computeType))
            {
                // Whether eval1 and eval2 are both feasible, or both
                // infeasible, dominance means FULL_SUCCESS, when their h are below hMax.
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
            else if (eval1->isFeasible(computeType) && eval2->isFeasible(computeType))
            {
                // Eval1 and eval2 are both feasible, but eval1 does
                // not dominate eval2.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else if (!eval1->isFeasible(computeType) && !eval2->isFeasible(computeType))
            {
                // Comparing two infeasible points
                if (eval1->getH(computeType) <= hMax
                    && eval1->getH(computeType) < eval2->getH(computeType))
                {
                    // Partial success (improving). Found an infeasible
                    // solution with a better h. f is worse or equivalent (as the dominance relation is false).
                    success = NOMAD::SuccessType::PARTIAL_SUCCESS;
                }
                else
                {
                    success = NOMAD::SuccessType::UNSUCCESSFUL;
                }
            }
        }
    }

    return success;
}


/*--------------------------------------------------*/
/*                      display                     */
/*--------------------------------------------------*/
std::string NOMAD::Eval::display(const NOMAD::FHComputeTypeS & computeType, const int prec) const
{
    std::string s;

    s += NOMAD::enumStr(_evalStatus);
    s += "\t ";

    try
    {
        if (computeType.computeType != NOMAD::ComputeType::STANDARD)
        {
            s += NOMAD::computeTypeToString(computeType.computeType) + " :";
        }
        const NOMAD::ArrayOfDouble& f = getFs(computeType);
        NOMAD::Double h = getH(computeType);
        if (f.isDefined())
        {
            s += "f = ";
            s += f.display(NOMAD::ArrayOfDouble(f.size(),prec));
        }
        else
        {
            s += "Undefined f";
        }
        s += "\t ";
        if (h.isDefined())
        {
            s += "h = ";
            s += h.display(prec);
            s += " (" + NOMAD::hNormTypeToString(computeType.hNormType) + ")";
        }
        else
        {
            s += "Undefined h";
        }
    }
    catch (NOMAD::Exception&)
    {
        // Could not compute f and h. Show raw bbo instead.
        s += getBBO();
    }

    return s;
}

NOMAD::ArrayOfDouble NOMAD::Eval::getBBOutputByType( const BBOutputType & bboType )
{
    std::vector<double> bbo;

    if (NOMAD::EvalStatusType::EVAL_OK == _evalStatus)
    {
        const NOMAD::ArrayOfDouble & allBBO = _bbOutput.getBBOAsArrayOfDouble();

        for (size_t i = 0 ; i < allBBO.size() ; i++)
        {
            if (_bbOutputTypeList[i] == bboType)
            {
                bbo.push_back(allBBO[i].todouble());
            }
        }
    }
    return NOMAD::ArrayOfDouble(bbo);
}

/*-------------------------------------*/
/* Convert an eval status to a string. */
/*-------------------------------------*/
std::string NOMAD::enumStr(const NOMAD::EvalStatusType evalStatus)
{
    std::string str;

    switch (evalStatus)
    {
        case NOMAD::EvalStatusType::EVAL_NOT_STARTED:
            str = "Evaluation not started";
            break;
        case NOMAD::EvalStatusType::EVAL_FAILED:
            str = "Evaluation failed";
            break;
        case NOMAD::EvalStatusType::EVAL_ERROR:
            str = "Evaluation error";
            break;
        case NOMAD::EvalStatusType::EVAL_USER_REJECTED:
            str = "Evaluation rejected by user (pre-eval; may be submitted again)";
            break;
        case NOMAD::EvalStatusType::EVAL_USER_ACCEPTED:
            str = "Evaluation accepted by user (pre-evaluation)";
            break;
        case NOMAD::EvalStatusType::EVAL_OK:
            str = "Evaluation OK";
            break;
        case NOMAD::EvalStatusType::EVAL_IN_PROGRESS:
            str = "Evaluation in progress";
            break;
        case NOMAD::EvalStatusType::EVAL_WAIT:
            str = "Waiting for evaluation in progress";
            break;
        case NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED:
            str = "Undefined evaluation status";
            break;
        default:
            str = "Unrecognized evaluation status";
            throw NOMAD::Exception(__FILE__, __LINE__, str);
    }

    return str;
}


// Output raw eval status
// Does not do the same as enumStr.
std::ostream& NOMAD::operator<<(std::ostream& out, const NOMAD::EvalStatusType &evalStatus)
{
    switch (evalStatus)
    {
        case NOMAD::EvalStatusType::EVAL_NOT_STARTED:
            out << "EVAL_NOT_STARTED";
            break;
        case NOMAD::EvalStatusType::EVAL_FAILED:
            out << "EVAL_FAILED";
            break;
        case NOMAD::EvalStatusType::EVAL_ERROR:
            out << "EVAL_ERROR";
            break;
        case NOMAD::EvalStatusType::EVAL_USER_REJECTED:
            out << "EVAL_USER_REJECTED";
            break;
        case NOMAD::EvalStatusType::EVAL_USER_ACCEPTED:
            out << "EVAL_USER_ACCEPTED";
            break;
        case NOMAD::EvalStatusType::EVAL_OK:
            out << "EVAL_OK";
            break;
        case NOMAD::EvalStatusType::EVAL_IN_PROGRESS:
            out << "EVAL_IN_PROGRESS";
            break;
        case NOMAD::EvalStatusType::EVAL_WAIT:
            out << "EVAL_WAIT";
            break;
        case NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED:
            out << "EVAL_STATUS_UNDEFINED";
            break;
        default:
            // Do not throw.
            std::cout << "Warning: Unknown eval status type" << std::endl;
            // We do not want a small mistake to disrupt the flow.
            break;
    }

    return out;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::EvalStatusType &evalStatus)
{
    std::string s;
    is >> s;

    // Remove BB_/SURROGATE_/MODEL_ (this may be used in cache file to show the different eval types)
    size_t indSep = s.find('_');
    if (indSep < std::string::npos && NOMAD::stringToEvalType(s.substr(0, indSep), true /* true: do not trigger exception */) != NOMAD::EvalType::UNDEFINED)
    {
        s.erase(0,indSep+1);
    }


    if ("EVAL_NOT_STARTED" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_NOT_STARTED;
    }
    else if ("EVAL_FAILED" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_FAILED;
    }
    else if ("EVAL_ERROR" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_ERROR;
    }
    else if ("EVAL_USER_REJECTED" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_USER_REJECTED;
    }
    else if ("EVAL_USER_ACCEPTED" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_USER_ACCEPTED;
    }
    else if ("EVAL_OK" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_OK;
    }
    else if ("EVAL_IN_PROGRESS" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_IN_PROGRESS;
    }
    else if ("EVAL_WAIT" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_WAIT;
    }
    else if ("EVAL_STATUS_UNDEFINED" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED;
    }
    else
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED;

        // Put back s to istream.
        for (unsigned i = 0; i < s.size(); i++)
        {
            is.unget();
        }
    }

    return is;

}
