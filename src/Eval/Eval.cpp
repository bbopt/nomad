/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

// static member initialization
std::function<NOMAD::SuccessType(const NOMAD::Eval* eval1, const NOMAD::Eval* eval2, const NOMAD::Double& hMax)> NOMAD::Eval::_computeSuccessType = NOMAD::Eval::defaultComputeSuccessType;

std::function<NOMAD::Double(const NOMAD::Eval& eval, const NOMAD::BBOutputTypeList &bbOutputTypeList)> NOMAD::Eval::_computeH = NOMAD::Eval::defaultComputeH;

std::function<NOMAD::Double(const NOMAD::BBOutputType &bbOutputType, size_t index, const NOMAD::Double& bbo)> NOMAD::Eval::_computeHComponent = NOMAD::Eval::defaultComputeHComponent;


/*---------------------------------------------------------------------*/
/*                            Constructor 1                            */
/*---------------------------------------------------------------------*/
// Note: NOMAD::Double() makes value = NOMAD::NaN; defined = false.
NOMAD::Eval::Eval()
  : _toBeRecomputed(false),
    _f(),
    _h(NOMAD::INF),
    _evalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _bbOutput("")
{
}


/*---------------------------------------------------------------------*/
/*                            Constructor 2                            */
/*---------------------------------------------------------------------*/
NOMAD::Eval::Eval(std::shared_ptr<NOMAD::EvalParameters> params,
                  const NOMAD::BBOutput &bbOutput)
  : _toBeRecomputed(true),
    _f(),
    _h(NOMAD::INF),
    _evalStatus(NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED),
    _bbOutput(bbOutput)
{
    auto bbOutputType = params->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    
    
    // Note: if bbOutput is not eval_ok, then _f and _h end up undefined Doubles.
    _f = computeF(bbOutputType);
    
    // Set H
    setH (_computeH(*this, bbOutputType));    
    _toBeRecomputed = false;
    
    if (_bbOutput.getEvalOk() && _f.isDefined())
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
  : _toBeRecomputed(eval._toBeRecomputed),
    _f(eval._f),
    _h(eval._h),
    _evalStatus(eval._evalStatus),
    _bbOutput(eval._bbOutput)
{
}


/*-----------------------*/
/*     Other methods     */
/*-----------------------*/
// Compute f
NOMAD::Double NOMAD::Eval::computeF(const NOMAD::BBOutputTypeList &bbOutputTypeList) const
{
    NOMAD::Double f;
    if (_bbOutput.getEvalOk())
    {
        f = _bbOutput.getObjective(bbOutputTypeList);
    }
    return f;
}


/// Compute h for a given BB Output Type list
NOMAD::Double NOMAD::Eval::defaultComputeH(const NOMAD::Eval& eval, const NOMAD::BBOutputTypeList &bbOutputTypeList)
{
    NOMAD::Double h = 0.0;
    const NOMAD::ArrayOfDouble bbo = eval.getBBOutput().getBBOAsArrayOfDouble();
    bool hPos = false;

    if (eval.getBBOutput().getEvalOk())
    {
        
        size_t bboIndex = 0;
        for (auto bbOutputType : bbOutputTypeList)
        {
            
            if ( NOMAD::BBOutputTypeIsConstraint(bbOutputType) )
            {
 
                /// The computeHComponent function has default version using cons*cons for PB constraint.
                /// The default function can be replaced by a custom one using the Eval::setComputeHComponent static function.
                NOMAD::Double hComp = _computeHComponent(bbOutputType , bboIndex, bbo[bboIndex]);
                
                /// Aggregate the H component for a given constraint to H.
                
                // Violated Extreme Barrier constraint:
                // Set h to infinity and break.
                if (NOMAD::INF == hComp)
                {
                    h = NOMAD::INF;
                    break;
                }
                if (bbo[bboIndex] > 0)
                {
                    hPos = true;
                    h += hComp;
                }
            }
            bboIndex++;
        }
    }

    // Failsafe: If at least one PB constraint is positive, h must be set
    // to at least epsilon so that the Eval is recognized as infeasible.
    // Catch cases such as constraint violated by 1e-8, which gives h = 1e-16 
    // which is considered as 0.
    if (hPos && (0 == h))
    {
        h = NOMAD::Double::getEpsilon();
    }


    return h;
}


/// Default computation of h component for a constraint (squared bbo for PB)
NOMAD::Double NOMAD::Eval::defaultComputeHComponent( const NOMAD::BBOutputType & bbOutputType , size_t index __attribute__((unused)), const NOMAD::Double &bbo )
{
    if ( ! NOMAD::BBOutputTypeIsConstraint(bbOutputType) )
    {
        std::string str = "H component must be computed from BB output that is a constraint";
        throw NOMAD::Exception(__FILE__, __LINE__, str);
    }
    
    NOMAD::Double h = 0.0;
    if (bbo > 0)
    {
        if (NOMAD::BBOutputType::EB == bbOutputType)
        {
            h = NOMAD::INF;
        }
        else if (NOMAD::BBOutputType::PB == bbOutputType)
        {
            h = bbo * bbo;
        }
    }
    
    return h;
}


NOMAD::Double NOMAD::Eval::computeHPB(const NOMAD::Eval& eval, const NOMAD::BBOutputTypeList &bbOutputTypeList)
{
    NOMAD::BBOutputTypeList bbOutputTypeListPB;
    for (auto bbOutputType : bbOutputTypeList)
    {
        if (NOMAD::BBOutputType::EB == bbOutputType)
        {
            // Replace EB by PB
            bbOutputTypeListPB.push_back(NOMAD::BBOutputType::PB);
        }
        else
        {
            bbOutputTypeListPB.push_back(bbOutputType);
        }
    }
    // Compute h with EB replaced by PB.
    // Required to revert to default compute H to have access to custom compute H component function.
    return defaultComputeH(eval, bbOutputTypeListPB);
}


bool NOMAD::Eval::isFeasible() const
{
    if (_toBeRecomputed)
    {
        std::cerr << "Warning: Eval::isFeasible() called on an Eval that needs to be recomputed." << std::endl;
    }

    // Comparison of NOMAD::Double accounts for epsilon
    return ( _h == 0.0 );
}


// This Eval has a status that permits its point to be re-evaluated.
bool NOMAD::Eval::canBeReEvaluated() const
{
    bool reEval = false;
    if (   _evalStatus == NOMAD::EvalStatusType::EVAL_OK
        || _evalStatus == NOMAD::EvalStatusType::EVAL_NOT_STARTED
        || _evalStatus == NOMAD::EvalStatusType::EVAL_USER_REJECTED
        || _evalStatus == NOMAD::EvalStatusType::EVAL_CONS_H_OVER
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
        || _evalStatus == NOMAD::EvalStatusType::EVAL_CONS_H_OVER
        || _evalStatus == NOMAD::EvalStatusType::EVAL_ERROR)
    {
        goodForCache = true;
    }
    return goodForCache;
}


/*----------------------------------------*/
/*      Get f. Warn if recompute needed.  */
/*----------------------------------------*/
NOMAD::Double NOMAD::Eval::getF() const
{
    if (_toBeRecomputed)
    {
        std::cerr << "Warning: Eval::getF() called on an Eval that needs to be recomputed." << std::endl;
    }
    return _f;
}


/*----------------------------------------*/
/*      Set f and update eval status      */
/*----------------------------------------*/
void NOMAD::Eval::setF(const NOMAD::Double &f)
{
    _f = f;

    _evalStatus = f.isDefined() ? NOMAD::EvalStatusType::EVAL_OK : NOMAD::EvalStatusType::EVAL_FAILED;
}


/*----------------------------------------*/
/*      Get h. Warn if recompute needed.  */
/*----------------------------------------*/
NOMAD::Double NOMAD::Eval::getH() const
{
    if (_toBeRecomputed)
    {
        std::cerr << "Warning: Eval::getH() called on an Eval that needs to be recomputed." << std::endl;
    }
    return _h;
}


/*--------------------------------------*/
/*      Set h and feasibility flag      */
/*--------------------------------------*/
void NOMAD::Eval::setH(const NOMAD::Double &h)
{
    if (h < 0)
    {
        std::string err = "Error: Trying to set a negative h (" + h.tostring() + ")";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    _h = h;
}


/*-----------------------*/
/*      Set BBOutput     */
/*-----------------------*/
void NOMAD::Eval::setBBOutput(const NOMAD::BBOutput &bbOutput)
{
    _bbOutput = bbOutput;
    toRecompute(true);
}


void NOMAD::Eval::setBBOutputAndRecompute(const NOMAD::BBOutput& bbOutput,
                                          const NOMAD::BBOutputTypeList& bbOutputType)
{
    setBBOutput(bbOutput);
    if (!bbOutput.checkSizeMatch(bbOutputType))
    {
        _evalStatus = NOMAD::EvalStatusType::EVAL_ERROR;
    }
    else
    {
        setF(computeF(bbOutputType));
        setH(_computeH(*this, bbOutputType));
    }
    toRecompute(false);
}


void NOMAD::Eval::setBBO(const std::string &bbo,
                         const NOMAD::BBOutputTypeList &bbOutputType,
                         const bool evalOk)
{
    _bbOutput.setBBO(bbo, evalOk);
    if (bbOutputType.size() > 0)
    {
        setF(computeF(bbOutputType));
        setH(_computeH(*this, bbOutputType));
        toRecompute(false);
    }
    else
    {
        toRecompute(true);
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

    if (this == &e)
    {
        equal = true;
    }
    else if (!_f.isDefined() || !e._f.isDefined())
    {
        // If either value is undefined, consider Evals not equal - even if both are undefined.
        equal = false;
    }
    else
    {
        // General case
        equal = ( (_f == e._f) && (_h == e._h) );
    }

    return equal;
}

/*--------------------------------------------------------------------------*/
/* Note: This operator must not be used to find/store EvalPoints in the     */
/* cache or in any set.                                                     */
/*--------------------------------------------------------------------------*/
bool NOMAD::Eval::operator<(const NOMAD::Eval &eval) const
{
    return dominates(eval);
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
bool NOMAD::Eval::dominates(const NOMAD::Eval &eval) const
{
    bool dom = false;

    NOMAD::Double f1 = getF();
    NOMAD::Double h1 = getH();
    NOMAD::Double f2 = eval.getF();
    NOMAD::Double h2 = eval.getH();

    if (isFeasible() && eval.isFeasible())
    {
        dom = (f1 < f2);
    }
    else if (!isFeasible() && !eval.isFeasible())
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
bool NOMAD::Eval::compEvalFindBest(const NOMAD::Eval &eval1, const NOMAD::Eval &eval2)
{
    // Success is PARTIAL_SUCCESS or FULL_SUCCESS
    // if eval1 is better than eval2.
    // hMax is ignored (set to NOMAD::INF).
    NOMAD::SuccessType success = _computeSuccessType(&eval1, &eval2, NOMAD::INF);

    return (success >= NOMAD::SuccessType::PARTIAL_SUCCESS);
}


NOMAD::SuccessType NOMAD::Eval::defaultComputeSuccessType(const Eval* eval1,
                                                          const Eval* eval2,
                                                          const Double& hMax)
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
            if (eval1->getH() > hMax)
            {
                // Even if eval2 is NULL, this case is not successful.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else
            {
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
        }
        else
        {
            if (eval1->dominates(*eval2))
            {
                // Whether eval1 and eval2 are both feasible, or both
                // infeasible, dominance means FULL_SUCCESS.
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
            else if (eval1->isFeasible() && eval2->isFeasible())
            {
                // Eval1 and eval2 are both feasible, but eval1 does
                // not dominate eval2.
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
            else if (!eval1->isFeasible() && !eval2->isFeasible())
            {
                // Comparing two infeasible points
                if (eval1->getH() <= hMax
                    && eval1->getH() < eval2->getH()
                    && eval1->getF() > eval2->getF())
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


NOMAD::SuccessType NOMAD::Eval::computeSuccessTypePhaseOne(const NOMAD::Eval* eval1,
                                                           const NOMAD::Eval* eval2,
                                                           const NOMAD::Double& hMax)
{
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    if (nullptr != eval1)
    {
        if (eval1->isFeasible())
        {
            success = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else if (nullptr == eval2)
        {
            success = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else
        {
            if (eval1->getH() < eval2->getH())
            {
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
            else
            {
                success = NOMAD::SuccessType::UNSUCCESSFUL;
            }
        }
    }

    return success;
}


/*--------------------------------------------------*/
/*                      display                     */
/*--------------------------------------------------*/
std::string NOMAD::Eval::display() const
{
    std::string s;

    s += NOMAD::enumStr(_evalStatus);
    s += "\t ";
    if (_f.isDefined())
    {
        s += "f = ";
        s += _f.tostring();
    }
    else
    {
        s += "Undefined f";
    }
    s += "\t ";
    if (_h.isDefined())
    {
        s += "h = ";
        s += _h.tostring();
    }
    else
    {
        s += "Undefined h";
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
        case NOMAD::EvalStatusType::EVAL_CONS_H_OVER:
            str = "Evaluation constraint violation is too high (may be submitted again)";
            break;
        case NOMAD::EvalStatusType::EVAL_OK:
            str = "Evaluation OK";
            break;
        case NOMAD::EvalStatusType::EVAL_IN_PROGRESS:
            str = "Evaluation in progress";
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
        case NOMAD::EvalStatusType::EVAL_CONS_H_OVER:
            out << "EVAL_CONS_H_OVER";
            break;
        case NOMAD::EvalStatusType::EVAL_OK:
            out << "EVAL_OK";
            break;
        case NOMAD::EvalStatusType::EVAL_IN_PROGRESS:
            out << "EVAL_IN_PROGRESS";
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
    else if ("EVAL_CONS_H_OVER" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_CONS_H_OVER;
    }
    else if ("EVAL_OK" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_OK;
    }
    else if ("EVAL_IN_PROGRESS" == s)
    {
        evalStatus = NOMAD::EvalStatusType::EVAL_IN_PROGRESS;
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

