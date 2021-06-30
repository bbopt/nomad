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
 \file   Eval.cpp
 \brief  Evaluation at a point (implementation)
 \author Viviane Rochon Montplaisir
 \date   March 2017
 \see    Eval.hpp
 */
#include "../Eval/Eval.hpp"


/*---------------------------------------------------------------------*/
/*                            Constructor 1                            */
/*---------------------------------------------------------------------*/
// Note: NOMAD::Double() makes value = NOMAD::NaN; defined = false.
NOMAD::Eval::Eval()
  : _evalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _bbOutput(""),
    _bbOutputTypeList(),
    _bbOutputComplete(false)
{
}


/*---------------------------------------------------------------------*/
/*                            Constructor 2                            */
/*---------------------------------------------------------------------*/
NOMAD::Eval::Eval(std::shared_ptr<NOMAD::EvalParameters> params,
                  const NOMAD::BBOutput &bbOutput)
  : _evalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _bbOutput(bbOutput),
    _bbOutputTypeList(params->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"))
{
    _bbOutputComplete = _bbOutput.isComplete(_bbOutputTypeList);

    NOMAD::Double f = _bbOutput.getObjective(_bbOutputTypeList);
    if (_bbOutput.getEvalOk() && f.isDefined())
    {
        _evalStatus = NOMAD::EvalStatusType::EVAL_OK;
    }
    else
    {
        _evalStatus = NOMAD::EvalStatusType::EVAL_FAILED;
    }
}


/*---------------------------------------------------------------------*/
/*                           Copy Constructor                          */
/*---------------------------------------------------------------------*/
NOMAD::Eval::Eval(const NOMAD::Eval &eval)
  : _evalStatus(eval._evalStatus),
    _bbOutput(eval._bbOutput),
    _bbOutputTypeList(eval._bbOutputTypeList),
    _bbOutputComplete(eval._bbOutputComplete)
{
}


/*-----------------------*/
/*     Other methods     */
/*-----------------------*/
bool NOMAD::Eval::isFeasible(const NOMAD::ComputeType& computeType) const
{
    if (NOMAD::EvalStatusType::EVAL_OK != _evalStatus)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Eval::isFeasible: Needs status type EVAL_OK");
    }
    NOMAD::Double h = getH(computeType);
    return (h.isDefined() && h.todouble() < NOMAD::Double::getEpsilon());
}


// This Eval has a status that permits its point to be re-evaluated.
bool NOMAD::Eval::canBeReEvaluated() const
{
    bool reEval = false;
    if (   _evalStatus == NOMAD::EvalStatusType::EVAL_OK
        || _evalStatus == NOMAD::EvalStatusType::EVAL_NOT_STARTED
        || _evalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
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
        || _evalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
        || _evalStatus == NOMAD::EvalStatusType::EVAL_ERROR)
    {
        goodForCache = true;
    }
    return goodForCache;
}


/*------------------------------------*/
/*      Get f. Always recomputed.     */
/*------------------------------------*/
NOMAD::Double NOMAD::Eval::getF(const NOMAD::ComputeType& computeType) const
{
    NOMAD::Double f;

    if (NOMAD::EvalStatusType::EVAL_OK != _evalStatus)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"getF(): EvalStatusType not EVAL_OK");
    }
    switch (computeType)
    {
        case NOMAD::ComputeType::STANDARD:
            f = _bbOutput.getObjective(_bbOutputTypeList);
            break;
        case NOMAD::ComputeType::PHASE_ONE:
            f = computeFPhaseOne();
            break;
        case NOMAD::ComputeType::USER:
            break;
        default:
            throw NOMAD::Exception(__FILE__,__LINE__,"getF(): ComputeType not supported");
    }

    return f;
}


/*-------------------------------------*/
/*      Get h. Always recomputed.      */
/*-------------------------------------*/
NOMAD::Double NOMAD::Eval::getH(const NOMAD::ComputeType& computeType) const
{
    NOMAD::Double h;

    if (NOMAD::EvalStatusType::EVAL_OK != _evalStatus)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"getH(): EvalStatusType not EVAL_OK: " + NOMAD::enumStr(_evalStatus));
    }
    switch (computeType)
    {
        case NOMAD::ComputeType::STANDARD:
            h = computeHStandard();
            break;
        case NOMAD::ComputeType::PHASE_ONE:
            h = 0.0;
            break;
        case NOMAD::ComputeType::USER:
            break;
        default:
            throw NOMAD::Exception(__FILE__,__LINE__,"getH(): ComputeType not supported");
    }

    return h;
}


NOMAD::Double NOMAD::Eval::computeHStandard() const
{
    NOMAD::Double h = 0.0;
    bool hPos = false;

    const NOMAD::ArrayOfDouble bboArray = _bbOutput.getBBOAsArrayOfDouble();
    size_t bboIndex = 0;
    for (auto bbOutputType : _bbOutputTypeList)
    {
        NOMAD::Double bboI = bboArray[bboIndex];
        bboIndex++;
        if (!NOMAD::BBOutputTypeIsConstraint(bbOutputType))
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
            if (NOMAD::BBOutputType::EB == bbOutputType)
            {
                hTemp = NOMAD::INF;
            }
            else if (NOMAD::BBOutputType::PB == bbOutputType)
            {
                hTemp = bboI * bboI;
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


NOMAD::Double NOMAD::Eval::computeFPhaseOne() const
{
    NOMAD::Double f ;
    const NOMAD::ArrayOfDouble bboArray = _bbOutput.getBBOAsArrayOfDouble();
    bool fPos = false;

    if (NOMAD::EvalStatusType::EVAL_OK == _evalStatus)
    {
        f=0.0;
        size_t bboIndex = 0;
        for (auto bbOutputType : _bbOutputTypeList)
        {
            NOMAD::Double bboI = bboArray[bboIndex];
            bboIndex++;
            if (NOMAD::BBOutputType::EB != bbOutputType)
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
                f += bboI * bboI;
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

    if (bbOutputTypeList.empty())
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
        _evalStatus = _bbOutput.getObjective(_bbOutputTypeList).isDefined() ? NOMAD::EvalStatusType::EVAL_OK : NOMAD::EvalStatusType::EVAL_FAILED;
    }

}


/*-----------------------------------------------------------*/
/*                           operator ==                     */
/*-----------------------------------------------------------*/
bool NOMAD::Eval::operator==(const NOMAD::Eval &e) const
{
    // Ignore eval status for comparison.
    // Compare f and h.

    bool equal = false;
    NOMAD::Double f1;
    NOMAD::Double f2;
    if (NOMAD::EvalStatusType::EVAL_OK == _evalStatus)
    {
        f1 = getF(NOMAD::ComputeType::STANDARD);
    }
    if (NOMAD::EvalStatusType::EVAL_OK == e._evalStatus)
    {
        f2 = e.getF(NOMAD::ComputeType::STANDARD);
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
        NOMAD::Double h1 = getH(NOMAD::ComputeType::STANDARD);
        NOMAD::Double h2 = e.getH(NOMAD::ComputeType::STANDARD);
        // As for f, if either h value is undefined, consider Evals not equal - even if both are undefined.
        if (!h1.isDefined() || !h2.isDefined())
        {
            equal = false;
        }
        else
        {
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
    return dominates(eval, NOMAD::ComputeType::STANDARD);
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
bool NOMAD::Eval::dominates(const NOMAD::Eval &eval, const NOMAD::ComputeType &computeType) const
{
    bool dom = false;

    NOMAD::Double f1 = getF(computeType);
    NOMAD::Double h1 = getH(computeType);
    NOMAD::Double f2 = eval.getF(computeType);
    NOMAD::Double h2 = eval.getH(computeType);

    if (isFeasible(computeType) && eval.isFeasible(computeType))
    {
        dom = (f1 < f2);
    }
    else if (!isFeasible(computeType) && !eval.isFeasible(computeType))
    {
        if (h1 != NOMAD::INF)
        {
            dom = (f1 <= f2) && (h1 <= h2) && ((f1 < f2) || (h1 < h2));
        }
    }
    // else - comparing a feasible point with an unfeasible point.
    // Always false. Do nothing.

    return dom;
}


// Comparison function used for Cache's findBest functions.
bool NOMAD::Eval::compEvalFindBest(const NOMAD::Eval &eval1, const NOMAD::Eval &eval2, const NOMAD::ComputeType& computeType)
{
    // Success is PARTIAL_SUCCESS or FULL_SUCCESS
    // if eval1 is better than eval2.
    // hMax is ignored (set to NOMAD::INF).
    NOMAD::SuccessType success = computeSuccessType(&eval1, &eval2, computeType, NOMAD::INF);

    return (success >= NOMAD::SuccessType::PARTIAL_SUCCESS);
}


/*
NOMAD::SuccessType NOMAD::Eval::computeSuccessType(const NOMAD::Eval* eval1,
                                                   const NOMAD::Eval* eval2,
                                                   const NOMAD::ComputeType& computeType,
                                                   const NOMAD::Double& hMax)
*/
// Comparison function used for Cache's findBest functions.
bool NOMAD::Eval::compInsertInBarrier(const NOMAD::Eval &eval1,
                                      const NOMAD::Eval &eval2,
                                      const NOMAD::ComputeType& computeType,
                                      NOMAD::SuccessType successType,
                                      bool strictEqual)
{
    // Success is PARTIAL_SUCCESS or FULL_SUCCESS
    // if eval1 is better than eval2.
    // hMax is ignored (set to NOMAD::INF).
    NOMAD::SuccessType success = computeSuccessType(&eval1, &eval2, computeType, NOMAD::INF);
    if (strictEqual)
    {
        return (success == successType);
    }
    else
    {
        return (success >= successType);
    }
}


// Comparison function used for Barrier update function.
bool NOMAD::Eval::compEvalBarrier(const NOMAD::Eval &eval1, const NOMAD::Eval &eval2)
{
    // - If eval1 domnates eval2, return true.
    // - Else, return true if eval1's f is better than eval2's.
    bool isBetter = false;
    if (eval1.dominates(eval2))
    {
        isBetter = true;
    }
    else if (eval2.dominates(eval1))
    {
        isBetter = false;
    }
    else if (eval1.getF() < eval2.getF())
    {
        isBetter = true;
    }
    else if (eval2.getF() < eval1.getF())
    {
        isBetter = false;
    }

    return isBetter;
}


NOMAD::SuccessType NOMAD::Eval::computeSuccessType(const NOMAD::Eval* eval1,
                                                   const NOMAD::Eval* eval2,
                                                   const NOMAD::ComputeType& computeType,
                                                   const NOMAD::Double& hMax)
{
    // NOT_EVALUATED,      // Not evaluated yet
    // UNSUCCESSFUL,       // Failure
    // PARTIAL_SUCCESS,    // Partial success (improving). Found an infeasible
    //                        solution with a better h. f is worse.
    // FULL_SUCCESS        // Full success (dominating)
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    if (nullptr != eval1)
    {
        if (nullptr == eval2)
        {
            if (eval1->getH(computeType) > hMax)
            {
                // Even if eval2 is NULL, this case is not successful.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else
            {
                // A new infeasible point, without prior infeasible point, is partial success,
                // not a full success.
                if (eval1->isFeasible())
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
            if (eval1->dominates(*eval2, computeType))
            {
                // Whether eval1 and eval2 are both feasible, or both
                // infeasible, dominance means FULL_SUCCESS.
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
                    && eval1->getH(computeType) < eval2->getH(computeType)
                    && eval1->getF(computeType) > eval2->getF(computeType))
                {
                    // Partial success (improving). Found an infeasible
                    // solution with a better h. f is worse.
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
std::string NOMAD::Eval::display(const NOMAD::ComputeType& computeType, const int prec) const
{
    std::string s;

    s += NOMAD::enumStr(_evalStatus);
    s += "\t ";

    try
    {
        NOMAD::Double f = getF(computeType);
        NOMAD::Double h = getH(computeType);
        if (f.isDefined())
        {
            s += "f = ";
            s += f.display(prec);
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
            str = "Evaluation rejected by user (may be submitted again)";
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
            std::cerr << "Warning: Unknown eval status type" << std::endl;
            // We do not want a small mistake to disrupt the flow.
            break;
    }

    return out;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::EvalStatusType &evalStatus)
{
    std::string s;
    is >> s;

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
        // Put back s to istream.
        for (unsigned i = 0; i < s.size(); i++)
        {
            is.unget();
        }
    }

    return is;

}

