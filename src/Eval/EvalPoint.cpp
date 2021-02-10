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
 \file   EvalPoint.cpp
 \brief  Evaluation point (implementation)
 \author Sebastien Le Digabel and Viviane Rochon Montplaisir
 \date   March 2017
 \see    EvalPoint.hpp
 */
#include "../Eval/EvalPoint.hpp"

size_t NOMAD::EvalPoint::_currentTag = 0;

/*---------------------------------------------------------------------*/
/*                            Constructor 1                            */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint ()
  : Point(),
    _eval(nullptr),
    _evalSgte(nullptr),
    _tag(0),
    _threadAlgo(NOMAD::getThreadNum()),
    _numberEval(0),
    _pointFrom(nullptr),
    _genStep("")
{
}


/*---------------------------------------------------------------------*/
/*                            Constructor 2                            */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint(size_t n)
  : NOMAD::Point(n),
    _eval(nullptr),
    _evalSgte(nullptr),
    _tag(0),
    _threadAlgo(NOMAD::getThreadNum()),
    _numberEval(0),
    _pointFrom(nullptr),
    _genStep("")
{
}


/*---------------------------------------------------------------------*/
/*                            Constructor 3                            */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint(const NOMAD::Point &x)
  : Point(x),
    _eval(nullptr),
    _evalSgte(nullptr),
    _tag(0),
    _threadAlgo(NOMAD::getThreadNum()),
    _numberEval(0),
    _pointFrom(nullptr),
    _genStep("")
{
}


/*---------------------------------------------------------------------*/
/*                           Copy Constructor                          */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint(const NOMAD::EvalPoint &evalPoint)
  : Point(evalPoint)
{
    copyMembers(evalPoint);
}


/*---------------------------------------------------------------------*/
/* Helper for copy constructor                                         */
/* Copy all members except Point.                                      */
/* Useful when converting from full to sub dimension and vice-versa.   */
/*---------------------------------------------------------------------*/
void NOMAD::EvalPoint::copyMembers(const NOMAD::EvalPoint &evalPoint)
{
    _tag = evalPoint._tag;
    _threadAlgo = evalPoint._threadAlgo;
    _numberEval = evalPoint._numberEval;

    _eval = nullptr;
    _evalSgte = nullptr;
    if (nullptr != evalPoint._eval)
    {
        // deep copy.
        _eval = NOMAD::EvalUPtr(new NOMAD::Eval(*evalPoint.getEval(NOMAD::EvalType::BB)));
    }
    if (nullptr != evalPoint._evalSgte)
    {
        // deep copy.
        _evalSgte = NOMAD::EvalUPtr(new NOMAD::Eval(*evalPoint.getEval(NOMAD::EvalType::SGTE)));
    }

    // shallow copy
    _pointFrom = evalPoint.getPointFrom();
    _genStep = evalPoint.getGenStep();
}


/*-----------------------------------------------------------*/
/*                     Affectation operator                  */
/*-----------------------------------------------------------*/
NOMAD::EvalPoint & NOMAD::EvalPoint::operator=(const NOMAD::EvalPoint &evalPoint)
{
    if (this == &evalPoint)
    {
        return *this;
    }

    Point::operator=(evalPoint);

    _tag = evalPoint._tag;
    _threadAlgo = evalPoint._threadAlgo;
    _numberEval = evalPoint._numberEval;

    _pointFrom = evalPoint._pointFrom;
    _genStep = evalPoint._genStep;

    // Do NOT delete _eval. Since it is a smart ptr, it will take care
    // of itself. Releasing the smart ptr here causes a memory leak.

    if (nullptr == evalPoint._eval)
    {
        _eval = nullptr;
    }
    else
    {
        // deep copy.
        _eval = NOMAD::EvalUPtr(new NOMAD::Eval(*evalPoint._eval));
    }

    if (nullptr == evalPoint._evalSgte)
    {
        _evalSgte = nullptr;
    }
    else
    {
        // deep copy.
        _evalSgte = NOMAD::EvalUPtr(new NOMAD::Eval(*evalPoint._evalSgte));
    }

    return *this;
}


/*---------------------------------------------------------------------*/
/*                               Destructor                            */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::~EvalPoint ()
{
    // Do NOT delete _eval. Since it is a smart ptr, it will take care
    // of itself. Releasing the smart ptr here causes a memory leak.
}


/*-----------------------------------------------------------*/
/*                           operator ==                     */
/*-----------------------------------------------------------*/
// Note: Considering both eval (bb) and evalSgte.
bool NOMAD::EvalPoint::operator== (const NOMAD::EvalPoint &evalPoint) const
{
    // First compare Points.
    bool equal = Point::operator==(evalPoint);

    // Ignore tag.
    // Ignore numberEval.
    // Ignore pointFrom and genStep.

    if (equal)
    {
        auto eval = getEval(NOMAD::EvalType::BB);
        auto eval2 = evalPoint.getEval(NOMAD::EvalType::BB);
        // Verify that evals are not to recompute. Otherwise, throw an exception.
        if (nullptr != eval && eval->toBeRecomputed())
        {
            std::string err = "Need to recompute f and h for this EvalPoint: ";
            err += this->display();
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }
        if (nullptr != eval2 && eval2->toBeRecomputed())
        {
            std::string err = "Need to recompute f and h for this EvalPoint: ";
            err += this->display();
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        // Compare Evals (bb).
        if (nullptr == eval && nullptr == eval2)
        {
            // Both Evals are NULL.
            equal = true;
        }
        else if (nullptr == eval || nullptr == eval2)
        {
            // One Eval is NULL, but not both.
            equal = false;
        }
        else
        {
            // General case
            equal = ( *eval == *(eval2) );
        }
    }
    if (equal)
    {
        // Compare Evals (sgte).
        auto eval = getEval(NOMAD::EvalType::SGTE);
        auto eval2 = evalPoint.getEval(NOMAD::EvalType::SGTE);
        if (nullptr == eval && nullptr == eval2)
        {
            // Both Evals are NULL.
            equal = true;
        }
        else if (nullptr == eval || nullptr == eval2)
        {
            // One Eval is NULL, but not both.
            equal = false;
        }
        else
        {
            // General case
            equal = ( *eval == *(eval2) );
        }
    }

    return equal;
}


/*--------------------------------*/
/*  Comparison operator '<':      */
/*  Compare eval values           */
/*      (f and h).                */
/*  To be used for filter         */
/*  Warning: only BB eval is      */
/*      considered                */
/*--------------------------------*/
bool NOMAD::EvalPoint::operator<(const NOMAD::EvalPoint & ep) const
{
    return this->dominates(ep, NOMAD::EvalType::BB);
}


/*---------------------*/
/* Other class methods */
/*---------------------*/
bool NOMAD::EvalPoint::isEvalOk(const NOMAD::EvalType& evalType) const
{
    bool ret = false;

    auto eval = getEval(evalType);
    if (eval)
    {
        ret = (NOMAD::EvalStatusType::EVAL_OK == eval->getEvalStatus());
    }

    return ret;
}


/*---------*/
/* Get/Set */
/*---------*/
NOMAD::Eval* NOMAD::EvalPoint::getEval(const NOMAD::EvalType& evalType) const
{
    NOMAD::Eval* eval = nullptr;

    switch (evalType)
    {
        case NOMAD::EvalType::SGTE:
            eval = _evalSgte.get();
            break;
        case NOMAD::EvalType::BB:
            eval = _eval.get();
            break;
        case NOMAD::EvalType::UNDEFINED:
        default:
            break;
    }

    return eval;
}


void NOMAD::EvalPoint::setEval(const NOMAD::Eval& eval,
                               const NOMAD::EvalType& evalType)
{
    // Do not release _eval before assigning it a new value.
    // It is a smart ptr, it will take care of itself.

    switch (evalType)
    {
        case NOMAD::EvalType::SGTE:
            _evalSgte = NOMAD::EvalUPtr(new NOMAD::Eval(eval));
            break;
        case NOMAD::EvalType::BB:
        default:
            _eval = NOMAD::EvalUPtr(new NOMAD::Eval(eval));
            break;
    }

}


NOMAD::Double NOMAD::EvalPoint::getF(const NOMAD::EvalType& evalType) const
{
    NOMAD::Double f;

    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        f = eval->getF();
    }

    return f;
}


void NOMAD::EvalPoint::setF(const NOMAD::Double f, const NOMAD::EvalType& evalType)
{
    auto eval = getEval(evalType);

    if (nullptr == eval)
    {
        std::string err = "Error: setting f on a null eval";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    eval->setF(f);
}


NOMAD::Double NOMAD::EvalPoint::getH(const NOMAD::EvalType& evalType) const
{
    NOMAD::Double h;

    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        h = eval->getH();
    }

    return h;
}


void NOMAD::EvalPoint::setH(const NOMAD::Double &h, const NOMAD::EvalType& evalType)
{
    auto eval = getEval(evalType);
    if (nullptr == eval)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error: setting h on an EvalPoint that has no f.");
    }

    eval->setH(h);
}


std::string NOMAD::EvalPoint::getBBO(const NOMAD::EvalType& evalType) const
{
    std::string bbo;
    auto eval = getEval(evalType);

    if (nullptr != eval)
    {
        bbo = eval->getBBOutput().getBBO();
    }

    return bbo;
}


void NOMAD::EvalPoint::setBBO(const std::string &bbo,
                              const NOMAD::BBOutputTypeList &bboutputtypes,
                              const NOMAD::EvalType& evalType,
                              const bool evalOk)
{
    auto eval = getEval(evalType);

    if (nullptr == eval)
    {
        switch (evalType)
        {
            // Create new Eval
            case NOMAD::EvalType::SGTE:
                _evalSgte = NOMAD::EvalUPtr(new NOMAD::Eval());
                break;
            case NOMAD::EvalType::BB:
            default:
                _eval = NOMAD::EvalUPtr(new NOMAD::Eval());
                break;
        }
        eval = getEval(evalType);
    }

    if (nullptr == eval)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "EvalPoint::setBBO: Could not create new Eval");
    }
    else
    {
        eval->setBBO(bbo, bboutputtypes, evalOk);
    }

}


void NOMAD::EvalPoint::setBBO(const std::string &bbo,
                              const std::string &sBBOutputTypes,
                              const NOMAD::EvalType& evalType,
                              const bool evalOk)
{
    NOMAD::BBOutputTypeList bboutputtypes = NOMAD::stringToBBOutputTypeList(sBBOutputTypes);
    setBBO(bbo, bboutputtypes, evalType, evalOk);
}


void NOMAD::EvalPoint::setBBO(const NOMAD::BBOutput& bbo,
                              const NOMAD::EvalType& evalType,
                              const bool evalOk)
{
    auto eval = getEval(evalType);

    if (nullptr == eval)
    {
        switch (evalType)
        {
            case NOMAD::EvalType::SGTE:
                _evalSgte = NOMAD::EvalUPtr(new NOMAD::Eval());
                break;
            case NOMAD::EvalType::BB:
            default:
                _eval = NOMAD::EvalUPtr(new NOMAD::Eval());
                break;
        }
        eval = getEval(evalType);
    }
    if (nullptr == eval)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "EvalPoint::setBBO: Could not create new Eval");
    }
    else
    {
        eval->setBBOutput(bbo);
    }
}


NOMAD::EvalStatusType NOMAD::EvalPoint::getEvalStatus(const NOMAD::EvalType& evalType) const
{
    NOMAD::EvalStatusType evalStatus = NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED;

    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        evalStatus = eval->getEvalStatus();
    }

    return evalStatus;
}


void NOMAD::EvalPoint::setEvalStatus(const NOMAD::EvalStatusType &evalStatus,
                                     const NOMAD::EvalType& evalType)
{
    auto eval = getEval(evalType);

    if (nullptr == eval)
    {
        switch (evalType)
        {
            case NOMAD::EvalType::SGTE:
                _evalSgte = NOMAD::EvalUPtr(new NOMAD::Eval());
                break;
            case NOMAD::EvalType::BB:
            default:
                _eval = NOMAD::EvalUPtr(new NOMAD::Eval());
                break;
        }
        eval = getEval(evalType);
    }

    if (nullptr == eval)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "EvalPoint::setEvalStatus: Could not create new Eval");
    }
    else
    {
        eval->setEvalStatus(evalStatus);
    }
}


// This method is declared const so we can use it inside a const method.
void NOMAD::EvalPoint::updateTag() const
{
    if (_tag==0)
    {
        _currentTag++;
        _tag = _currentTag;
    }
}


void NOMAD::EvalPoint::resetCurrentTag()
{
    _currentTag = 0;
}


const std::shared_ptr<NOMAD::Point> NOMAD::EvalPoint::getPointFrom(const NOMAD::Point& fixedVariable) const
{
    auto pointFrom = _pointFrom;
    if (nullptr != pointFrom)
    {
        pointFrom = std::make_shared<NOMAD::Point>(pointFrom->projectPointToSubspace(fixedVariable));
    }

    return pointFrom;
}


void NOMAD::EvalPoint::setPointFrom(const std::shared_ptr<NOMAD::Point> pointFrom, const NOMAD::Point& fixedVariable)
{
    auto pointFromFull = pointFrom;
    if (pointFromFull->size() < fixedVariable.size())
    {
        // pointFrom must always be in full dimension. Convert if needed.
        pointFromFull = std::make_shared<NOMAD::Point>(pointFromFull->makeFullSpacePointFromFixed(fixedVariable));
    }

    _pointFrom = pointFromFull;
}


void NOMAD::EvalPoint::setGenStep(const std::string& genStep)
{
    if (!_genStep.empty() && _genStep != genStep)
    {
        // Prepend genStep for more information. For instance:
        // MegaSearchPoll - LH Search Method
        // instead of only MegaSearchPoll.
        _genStep = genStep + " - " + _genStep;
    }
    else
    {
        _genStep = genStep;
    }
}


bool NOMAD::EvalPoint::isFeasible(const NOMAD::EvalType& evalType) const
{
    bool feas = false;

    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        feas = eval->isFeasible();
    }

    return feas;
}


void NOMAD::EvalPoint::recomputeFH(const NOMAD::BBOutputTypeList &bbOutputType)
{
    // Recompute evals for all EvalTypes.

    // Recompute for blackbox
    auto eval = getEval(NOMAD::EvalType::BB);
    if (nullptr != eval)
    {
        auto bbo = eval->getBBOutput();
        eval->setBBOutputAndRecompute(bbo, bbOutputType);
    }

    // Recompute for SGTE
    eval = getEval(NOMAD::EvalType::SGTE);
    if (nullptr != eval)
    {
        auto bbo = eval->getBBOutput();
        eval->setBBOutputAndRecompute(bbo, bbOutputType);
    }
}


NOMAD::EvalPoint NOMAD::EvalPoint::makeFullSpacePointFromFixed(const NOMAD::Point &fixedVariable) const
{
    NOMAD::EvalPoint fullSpaceEvalPoint(getX()->makeFullSpacePointFromFixed(fixedVariable));
    fullSpaceEvalPoint.copyMembers(*this);

    return fullSpaceEvalPoint;
}


NOMAD::EvalPoint NOMAD::EvalPoint::makeSubSpacePointFromFixed(const NOMAD::Point &fixedVariable) const
{
    NOMAD::EvalPoint subSpaceEvalPoint(getX()->makeSubSpacePointFromFixed(fixedVariable));
    subSpaceEvalPoint.copyMembers(*this);

    return subSpaceEvalPoint;
}


// Should we evaluate (possibly re-evaluate) this point?
bool NOMAD::EvalPoint::toEval(short maxPointEval, const NOMAD::EvalType& evalType) const
{

    bool reEval = false;

    auto eval = getEval(evalType);
    if (nullptr == eval)
    {
        // No eval, return true.
        reEval = true;
    }
    else if (_numberEval >= maxPointEval)
    {
        // Too many evaluations, return false.
        reEval = false;
    }
    else if (NOMAD::EvalType::SGTE == evalType)
    {
        // If using sgte, never allow re-evaluation.
        reEval = false;
    }
    else if (_numberEval >= 1 && NOMAD::EvalStatusType::EVAL_OK == eval->getEvalStatus())
    {
        // For now, we will not re-evaluate an EvalPoint that is EVAL_OK.
        reEval = false;
    }
    else
    {
        reEval = eval->canBeReEvaluated();
    }

    return reEval;
}


// Not displaying evalSgte, only bb eval
std::string NOMAD::EvalPoint::display(const NOMAD::ArrayOfDouble &format) const
{
    std::string s = "#" + std::to_string(_tag) + " ";
    s += NOMAD::Point::display(format);
    if (nullptr != _eval)
    {
        s += "\t";
        s += _eval->display();
    }
    return s;
}


// Show both eval and evalSgte. For debugging purposes.
std::string NOMAD::EvalPoint::displayAll() const
{
    std::string s = "#" + std::to_string(_tag) + " ";
    s += NOMAD::Point::display();
    if (nullptr != _eval)
    {
        s += "\t";
        s += "(BB - ";
        s += _eval->display();
        s += ")";
    }
    if (nullptr != _evalSgte)
    {
        s += "\t";
        s += "(SGTE - ";
        s += _evalSgte->display();
        s += ")";
    }
    return s;
}


// Determine if an evalpoint has a sgte eval.
bool NOMAD::EvalPoint::hasSgteEval(const NOMAD::EvalPoint& evalPoint)
{
    return (nullptr != evalPoint.getEval(NOMAD::EvalType::SGTE));
}


// Determine if an evalpoint has a bb (regular) eval.
bool NOMAD::EvalPoint::hasBbEval(const NOMAD::EvalPoint& evalPoint)
{
    return (nullptr != evalPoint.getEval(NOMAD::EvalType::BB));
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::EvalPoint &evalPoint)
{
    // Example:
    // ( 1.7 2.99 -2.42 2.09 -36 2.33 ) EVAL_FAILED NaN 0

    NOMAD::Point p = *(evalPoint.getX());
    // Since this operator is used to write cache, we need full precision on point.
    os << p.display(NOMAD::ArrayOfDouble(evalPoint.size(), NOMAD::DISPLAY_PRECISION_FULL));

    // Never use sgte in input/output stream operators
    const NOMAD::Eval* eval = evalPoint.getEval(NOMAD::EvalType::BB);
    if (nullptr != eval)
    {
        os << " " << eval->getEvalStatus();     // Raw, ex. "EVAL_OK"
        os << " " << NOMAD::BBOutput::bboStart << " " << eval->getBBO();
        os << " " << NOMAD::BBOutput::bboEnd;
    }


    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::EvalPoint &evalPoint)
{
    // Set up structures to gather member info
    NOMAD::Point point;
    NOMAD::EvalStatusType evalStatus = NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED;
    bool skip = false;

    std::string s;
    is >> s;

    if (s.empty() || !is.good() || is.eof())
    {
        skip = true;
    }

    if (!skip && NOMAD::ArrayOfDouble::pStart == s)
    {
        // Found start of point.
        is.unget();
        is >> point;

        evalPoint = NOMAD::EvalPoint(point);

        // Read Eval - if following field is an EvalStatus.
        is >> evalStatus;
        if (NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED != evalStatus)
        {
            // Never use sgte in input/output stream operators
            evalPoint.setEvalStatus(evalStatus, NOMAD::EvalType::BB);

            // Read BBOutput
            NOMAD::BBOutput bbo("");
            is >> bbo;

            evalPoint.setBBO(bbo, NOMAD::EvalType::BB);  // BBO is set but f and h need to be recomputed

            // For now, set numEval to 1 if Eval exists. Currently,
            // only 1 Eval is correctly supported.
            evalPoint.setNumberEval(1);
        }
    }
    else if (!skip)
    {
        is.setstate(std::ios::failbit);
        std::string err = "Expecting \"" + NOMAD::ArrayOfDouble::pStart + "\", got \"" + s + "\"";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }


    return is;
}


bool NOMAD::findInList(const NOMAD::Point& point,
                       const std::vector<NOMAD::EvalPoint>& evalPointList,
                       NOMAD::EvalPoint& foundEvalPoint)
{
    bool found = false;

    for (auto evalPoint : evalPointList)
    {
        if (point == *evalPoint.getX())
        {
            foundEvalPoint = evalPoint;
            found = true;
            break;
        }
    }

    return found;
}


void NOMAD::convertPointListToSub(NOMAD::EvalPointList &evalPointList, const NOMAD::Point& fixedVariable)
{
    if (fixedVariable.isEmpty())
    {
        std::string s = "Error: Fixed variable of dimension 0";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
    for (size_t i = 0; i < evalPointList.size(); i++)
    {
        if (evalPointList[i].size() == fixedVariable.size())
        {
            evalPointList[i] = evalPointList[i].makeSubSpacePointFromFixed(fixedVariable);
        }
    }
}


void NOMAD::convertPointListToFull(NOMAD::EvalPointList &evalPointList, const NOMAD::Point& fixedVariable)
{
    for (size_t i = 0; i < evalPointList.size(); i++)
    {
        if (evalPointList[i].size() == fixedVariable.size() - fixedVariable.nbDefined())
        {
            evalPointList[i] = evalPointList[i].makeFullSpacePointFromFixed(fixedVariable);
        }
    }
}


#ifdef USE_UNORDEREDSET
/// Used for unordered set.
/// Template specialization of std::hash<class T> to std::hash<NOMAD@::EvalPoint>
size_t std::hash<NOMAD::EvalPoint>::operator()(const NOMAD::EvalPoint& evalPoint) const
{
    double eps = NOMAD::Double::getEpsilon();
    double sizeMax = SIZE_MAX;
    size_t hashKey = 0;
    double hashKeyDouble = hashKey; // Compute hashKey in double.

    for (size_t i = 0; i < evalPoint.size(); i++)
    {
        hashKeyDouble *= 10.0;
        double t = std::trunc(std::fabs(evalPoint[i].todouble() / eps));
        hashKeyDouble += t;

        // avoid size_t overflow, which would set hashKey to 0.
        while (hashKeyDouble > sizeMax)
        {
            hashKeyDouble -= sizeMax;
        }
    }

    // convert to size_t
    hashKey = size_t(hashKeyDouble);

    return hashKey;
}


/// Used for unordered set.
/// Template specialization of std::equal_to<class T> to equal_to<NOMAD::EvalPoint>
bool std::equal_to<NOMAD::EvalPoint>::operator()(const NOMAD::EvalPoint& lhs, const NOMAD::EvalPoint& rhs) const
{

    return (lhs == rhs);
}
#endif // USE_UNORDERED_SET


bool NOMAD::EvalPoint::dominates(const NOMAD::EvalPoint &ep,
                                 const NOMAD::EvalType& evalType) const
{
    bool dom = false;
    if (this != &ep && nullptr != getEval(evalType) && nullptr != ep.getEval(evalType))
    {
        dom = getEval(evalType)->dominates(*ep.getEval(evalType));
    }

    return dom;
}


