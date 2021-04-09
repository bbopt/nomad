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

int NOMAD::EvalPoint::_currentTag = -1;

/*---------------------------------------------------------------------*/
/*                            Constructor 1                            */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint ()
  : NOMAD::Point(),
    _eval(nullptr),
    _evalModel(nullptr),
    _tag(-1),
    _threadAlgo(NOMAD::getThreadNum()),
    _numberEval(0),
    _pointFrom(nullptr),
    _genStep(""),
    _direction(nullptr),
    _angle()
{
}


/*---------------------------------------------------------------------*/
/*                            Constructor 2                            */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint(size_t n)
  : NOMAD::Point(n),
    _eval(nullptr),
    _evalModel(nullptr),
    _tag(-1),
    _threadAlgo(NOMAD::getThreadNum()),
    _numberEval(0),
    _pointFrom(nullptr),
    _genStep(""),
    _direction(nullptr),
    _angle()
{
}


/*---------------------------------------------------------------------*/
/*                            Constructor 3                            */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint(const NOMAD::Point &x)
  : NOMAD::Point(x),
    _eval(nullptr),
    _evalModel(nullptr),
    _tag(-1),
    _threadAlgo(NOMAD::getThreadNum()),
    _numberEval(0),
    _pointFrom(nullptr),
    _genStep(""),
    _direction(nullptr),
    _angle()
{
}


/*---------------------------------------------------------------------*/
/*                           Copy Constructor                          */
/*---------------------------------------------------------------------*/
NOMAD::EvalPoint::EvalPoint(const NOMAD::EvalPoint &evalPoint)
  : NOMAD::Point(evalPoint)
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
    _evalModel = nullptr;
    if (nullptr != evalPoint._eval)
    {
        // deep copy.
        _eval = NOMAD::EvalUPtr(new NOMAD::Eval(*evalPoint.getEval(NOMAD::EvalType::BB)));
    }
    if (nullptr != evalPoint._evalModel)
    {
        // deep copy.
        _evalModel = NOMAD::EvalUPtr(new NOMAD::Eval(*evalPoint.getEval(NOMAD::EvalType::MODEL)));
    }

    // shallow copy
    _pointFrom = evalPoint.getPointFrom();
    _genStep = evalPoint.getGenStep();
    _direction = evalPoint.getDirection();
    _angle = evalPoint.getAngle();
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

    NOMAD::Point::operator=(evalPoint);

    _tag = evalPoint._tag;
    _threadAlgo = evalPoint._threadAlgo;
    _numberEval = evalPoint._numberEval;

    _pointFrom = evalPoint._pointFrom;
    _genStep = evalPoint._genStep;
    _direction = evalPoint._direction;
    _angle = evalPoint._angle;

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

    if (nullptr == evalPoint._evalModel)
    {
        _evalModel = nullptr;
    }
    else
    {
        // deep copy.
        _evalModel = NOMAD::EvalUPtr(new NOMAD::Eval(*evalPoint._evalModel));
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
// Note: Considering both eval (bb) and evalModel.
bool NOMAD::EvalPoint::operator== (const NOMAD::EvalPoint &evalPoint) const
{
    // First compare Points.
    bool equal = NOMAD::Point::operator==(evalPoint);

    // Ignore tag.
    // Ignore numberEval.
    // Ignore pointFrom and genStep.

    if (equal)
    {
        auto eval = getEval(NOMAD::EvalType::BB);
        auto eval2 = evalPoint.getEval(NOMAD::EvalType::BB);

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
        // Compare Evals (model).
        auto eval = getEval(NOMAD::EvalType::MODEL);
        auto eval2 = evalPoint.getEval(NOMAD::EvalType::MODEL);
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
    return this->dominates(ep, NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD);
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
        case NOMAD::EvalType::MODEL:
            eval = _evalModel.get();
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
        case NOMAD::EvalType::MODEL:
            _evalModel = NOMAD::EvalUPtr(new NOMAD::Eval(eval));
            break;
        case NOMAD::EvalType::BB:
        default:
            _eval = NOMAD::EvalUPtr(new NOMAD::Eval(eval));
            break;
    }

}


NOMAD::Double NOMAD::EvalPoint::getF(const NOMAD::EvalType& evalType,
                                     const NOMAD::ComputeType& computeType) const
{
    NOMAD::Double f;

    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        f = eval->getF(computeType);
    }

    return f;
}


NOMAD::Double NOMAD::EvalPoint::getH(const NOMAD::EvalType& evalType,
                                     const NOMAD::ComputeType& computeType) const
{
    NOMAD::Double h;

    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        h = eval->getH(computeType);
    }

    return h;
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
                              const NOMAD::BBOutputTypeList &bbOutputTypeList,
                              const NOMAD::EvalType& evalType,
                              const bool evalOk)
{
    auto eval = getEval(evalType);

    if (nullptr == eval)
    {
        switch (evalType)
        {
            // Create new Eval
            case NOMAD::EvalType::MODEL:
                _evalModel = NOMAD::EvalUPtr(new NOMAD::Eval());
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
        eval->setBBO(bbo, bbOutputTypeList, evalOk);
    }

}


void NOMAD::EvalPoint::setBBO(const std::string &bbo,
                              const std::string &sBBOutputTypes,
                              const NOMAD::EvalType& evalType,
                              const bool evalOk)
{
    NOMAD::BBOutputTypeList bbOutputTypeList = NOMAD::stringToBBOutputTypeList(sBBOutputTypes);
    setBBO(bbo, bbOutputTypeList, evalType, evalOk);
}


void NOMAD::EvalPoint::setBBOutputType(const NOMAD::BBOutputTypeList& bbOutputType,
                                       const NOMAD::EvalType& evalType)
{
    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        eval->setBBOutputTypeList(bbOutputType);
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
            case NOMAD::EvalType::MODEL:
                _evalModel = NOMAD::EvalUPtr(new NOMAD::Eval());
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
    if (-1 == _tag)
    {
        _currentTag++;
        _tag = _currentTag;
    }
}


void NOMAD::EvalPoint::resetCurrentTag()
{
    _currentTag = -1;
}


const std::shared_ptr<NOMAD::EvalPoint> NOMAD::EvalPoint::getPointFrom(const NOMAD::Point& fixedVariable) const
{
    auto pointFrom = _pointFrom;
    if (nullptr != pointFrom)
    {
        pointFrom = std::make_shared<NOMAD::EvalPoint>(pointFrom->projectPointToSubspace(fixedVariable));
    }

    return pointFrom;
}


void NOMAD::EvalPoint::setPointFrom(const std::shared_ptr<NOMAD::EvalPoint> pointFrom, const NOMAD::Point& fixedVariable)
{
    auto pointFromFull = pointFrom;
    if (pointFromFull->size() < fixedVariable.size())
    {
        // pointFrom must always be in full dimension. Convert if needed.
        pointFromFull = std::make_shared<NOMAD::EvalPoint>(pointFromFull->makeFullSpacePointFromFixed(fixedVariable));
    }

    _pointFrom = pointFromFull;

    // Also set Direction.
    if (nullptr != pointFromFull)
    {
        NOMAD::Point pointFull(*getX());
        if (pointFull.size() < fixedVariable.size())
        {
            pointFull = pointFull.makeFullSpacePointFromFixed(fixedVariable);
        }
        _direction = std::make_shared<NOMAD::Direction>(NOMAD::Point::vectorize(*pointFromFull, pointFull));
    }
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


bool NOMAD::EvalPoint::isFeasible(const NOMAD::EvalType& evalType, const NOMAD::ComputeType& computeType) const
{
    bool feas = false;

    auto eval = getEval(evalType);
    if (nullptr != eval)
    {
        feas = eval->isFeasible(computeType);
    }

    return feas;
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
    else if (NOMAD::EvalType::MODEL == evalType)
    {
        // If using model, never allow re-evaluation.
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


// Displaying only bb eval
std::string NOMAD::EvalPoint::display(const NOMAD::ComputeType& computeType,
                                      const NOMAD::ArrayOfDouble &pointFormat,
                                      const int &solFormat) const
{
    std::string s;
    if (_tag >= 0)
    {
        s = "#" + std::to_string(_tag) + " ";
    }
    s += NOMAD::Point::display(pointFormat);
    if (nullptr != _eval)
    {
        s += "\t";
        s += _eval->display(computeType, solFormat);
    }
    return s;
}


std::string NOMAD::EvalPoint::display(const NOMAD::ArrayOfDouble &pointFormat,
                                      const int &solFormat) const
{
    return display(NOMAD::ComputeType::STANDARD, pointFormat, solFormat);
}


// Show all evals. For debugging purposes.
std::string NOMAD::EvalPoint::displayAll(const NOMAD::ComputeType& computeType) const
{
    std::string s;
    if (_tag >= 0)
    {
        s = "#" + std::to_string(_tag) + " ";
    }
    s += NOMAD::Point::display();
    if (nullptr != _eval)
    {
        s += "\t";
        s += "(BB - ";
        s += _eval->display(computeType);
        s += ")";
    }
    if (nullptr != _evalModel)
    {
        s += "\t";
        s += "(MODEL - ";
        s += _evalModel->display(computeType);
        s += ")";
    }
    return s;
}


// Determine if an evalpoint has a model eval.
bool NOMAD::EvalPoint::hasModelEval(const NOMAD::EvalPoint& evalPoint)
{
    return (nullptr != evalPoint.getEval(NOMAD::EvalType::MODEL));
}


// Determine if an evalpoint has a bb (regular) eval.
bool NOMAD::EvalPoint::hasBbEval(const NOMAD::EvalPoint& evalPoint)
{
    return (nullptr != evalPoint.getEval(NOMAD::EvalType::BB));
}


bool NOMAD::EvalPoint::isPhaseOneSolution(const NOMAD::EvalPoint& evalPoint)
{
    bool issol = false;

    auto eval = evalPoint.getEval(NOMAD::EvalType::BB);
    if (nullptr != eval)
    {
        issol = (0.0 == eval->getF(NOMAD::ComputeType::PHASE_ONE).todouble());
    }

    return issol;
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::EvalPoint &evalPoint)
{
    // Example:
    // ( 1.7 2.99 -2.42 2.09 -36 2.33 ) EVAL_FAILED ( NaN 0 -20 )

    NOMAD::Point p = *(evalPoint.getX());
    // Since this operator is used to write cache, we need full precision on point.
    os << p.display(NOMAD::ArrayOfDouble(evalPoint.size(), NOMAD::DISPLAY_PRECISION_FULL));

    // Never use model eval in input/output stream operators
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
            // Never use model eval in input/output stream operators
            evalPoint.setEvalStatus(evalStatus, NOMAD::EvalType::BB);

            // Read BBOutput
            NOMAD::BBOutput bbo("");
            is >> bbo;

            evalPoint.setBBO(bbo.getBBO(), NOMAD::BBOutputTypeList(), NOMAD::EvalType::BB);

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
                                 const NOMAD::EvalType& evalType,
                                 const NOMAD::ComputeType& computeType) const
{
    bool dom = false;
    if (this != &ep && nullptr != getEval(evalType) && nullptr != ep.getEval(evalType))
    {
        dom = getEval(evalType)->dominates(*ep.getEval(evalType), computeType);
    }

    return dom;
}


