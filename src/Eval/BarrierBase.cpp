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

#include "../Cache/CacheBase.hpp"
#include "../Eval/BarrierBase.hpp"
#include "../Eval/ComputeSuccessType.hpp"
#include "../Output/OutputQueue.hpp"

void NOMAD::BarrierBase::setN()
{
    bool isSet = false;
    std::string s;

    for (const auto &evalPoint : getAllPointsPtr())
    {
        if (!isSet)
        {
            _n = evalPoint->size();
            isSet = true;
        }
        else if (evalPoint->size() != _n)
        {
            s = "Barrier has points of size " + std::to_string(_n) + " and of size ";
            s += std::to_string(evalPoint->size());
            throw NOMAD::Exception(__FILE__, __LINE__, s);
        }
    }
    if (!isSet)
    {
        s = "Barrier could not set point size";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }


}


// Verify there is a Cache instantiated.
void NOMAD::BarrierBase::checkCache()
{
    try
    {
        NOMAD::CacheBase::getInstance();
    }
    catch (NOMAD::Exception&)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Cache must be instantiated before initializing Barrier.");
    }
}

void NOMAD::BarrierBase::checkHMax()
{
    if (!_hMax.isDefined())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Barrier: hMax is not defined.");
    }
    if (_hMax < NOMAD::Double::getEpsilon())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Barrier: hMax must be positive. Value: " + _hMax.display());
    }
}



void NOMAD::BarrierBase::clearXFeas()
{
    _xFeas.clear();
}


void NOMAD::BarrierBase::clearXInf()
{
    _xInf.clear();
    _xIncInf.clear();
}

std::vector<NOMAD::EvalPoint> NOMAD::BarrierBase::getAllPoints() const
{
    std::vector<NOMAD::EvalPoint> allPoints;

    for (const auto & p: _xFeas )
    {
        allPoints.push_back(*p);
    }
    for (const auto & p: _xInf )
    {
        allPoints.push_back(*p);
    }
    return allPoints;
}

std::vector<NOMAD::EvalPointPtr> NOMAD::BarrierBase::getAllPointsPtr() const
{
    std::vector<NOMAD::EvalPointPtr> allPoints;

    allPoints.reserve(_xFeas.size() + _xInf.size()); // preallocate memory
    allPoints.insert(allPoints.end(), _xFeas.begin(), _xFeas.end());
    allPoints.insert(allPoints.end(), _xInf.begin(), _xInf.end());

    return allPoints;
}

const NOMAD::EvalPointPtr NOMAD::BarrierBase::getFirstPoint() const
{
    if (_xIncFeas.size() > 0)
    {
        return _xIncFeas[0];
    }
    else if (_xFeas.size() > 0)
    {
        return _xFeas[0];
    }
    else if (_xIncInf.size() > 0)
    {
        return _xIncInf[0];
    }
    else if (_xInf.size() > 0 )
    {
        return _xInf[0];
    }
    else
    {
        return nullptr;
    }
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::BarrierBase& barrier)
{
    auto allXFeas = barrier.getAllXFeas();
    auto allXInf  = barrier.getAllXInf();

    for (const auto & xFeas : allXFeas)
    {
        os << "X_FEAS " << *xFeas << std::endl;
    }
    for (const auto & xInf : allXInf)
    {
        os << "X_INF " << *xInf << std::endl;
    }
    os << "H_MAX " << barrier.getHMax() << std::endl;

    return os;
}

std::istream& NOMAD::operator>>(std::istream& is, NOMAD::BarrierBase& barrier)
{
    // Set up structures to gather member info
    NOMAD::EvalPoint xFeas, xInf;
    NOMAD::Double hMax = NOMAD::INF;
    barrier.clearXFeas();
    barrier.clearXInf();
    
    std::vector<NOMAD::EvalPoint> evalPointList;

    // Read line by line
    // We suppose that we have blackbox evaluation
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
        if ("X_FEAS" == name)
        {
            is >> xFeas;
            // Looking for xFeas in cache will ensure its f and h are updated
            NOMAD::CacheBase::getInstance()->find(xFeas, xFeas);
            // EvalType undefined: No check will be done on the feasibility
            evalPointList.push_back(xFeas);
            
        }
        else if ("X_INF" == name)
        {
            is >> xInf;
            // Looking for xInf in cache will ensure its f and h are updated
            NOMAD::CacheBase::getInstance()->find(xInf, xInf);
            evalPointList.push_back(xInf);
        }
        else if ("H_MAX" == name)
        {
            is >> hMax;
        }
        else
        {
            for (unsigned i = 0; i < name.size(); i++)
            {
                is.unget();
            }
            break;
        }
    }
    
    barrier.updateWithPoints(evalPointList,NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD, false, true /* true: update incumbents and hMax */);
    return is;
}


std::vector<NOMAD::EvalPointPtr>::iterator NOMAD::BarrierBase::findEvalPoint(std::vector<NOMAD::EvalPointPtr>::iterator first, std::vector<NOMAD::EvalPointPtr>::iterator last, const NOMAD::EvalPoint & p  )
{
    while (first!=last)
    {
        if (**first==p)
            return first;
        ++first;
    }
    return last;
}


bool NOMAD::BarrierBase::findPoint(const NOMAD::Point& point,
                                   NOMAD::EvalPoint& foundEvalPoint)
{
    bool found = false;

    auto evalPointList = getAllPointsPtr();
    for (const auto &evalPoint : evalPointList)
    {
        if (point.size() != evalPoint->size())
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Error: Eval points have different dimensions");
        }
        if (point == *evalPoint->getX())
        {
            foundEvalPoint = *evalPoint;
            found = true;
            break;
        }
    }

    return found;
}

void NOMAD::BarrierBase::checkXInf(const NOMAD::EvalPoint &xInf, NOMAD::EvalType evalType)
{
    // If evalType is UNDEFINED, skip this check.
    if (NOMAD::EvalType::UNDEFINED != evalType)
    {
        if (nullptr == xInf.getEval(evalType))
        {
            throw NOMAD::Exception(__FILE__, __LINE__,
                                   "Barrier: xInf must be evaluated before being set.");
        }
    }
}


void NOMAD::BarrierBase::checkXFeas(const NOMAD::EvalPoint &xFeas,
                                    NOMAD::EvalType  evalType,
                                    NOMAD::ComputeType computeType)
{
    // If evalType is UNDEFINED, skip this check.
    if (NOMAD::EvalType::UNDEFINED != evalType)
    {
        if (nullptr == xFeas.getEval(evalType))
        {
            throw NOMAD::Exception(__FILE__, __LINE__,
                                "Barrier: xFeas must be evaluated before being set.");
        }
        checkXFeasIsFeas(xFeas, evalType, computeType);
    }
}


void NOMAD::BarrierBase::checkXFeasIsFeas(const NOMAD::EvalPoint &xFeas,
                                          NOMAD::EvalType  evalType,
                                          NOMAD::ComputeType computeType)
{
    // If evalType is UNDEFINED, skip this check.
    if (NOMAD::EvalType::UNDEFINED != evalType)
    {
        auto eval = xFeas.getEval(evalType);
        if (nullptr != eval && NOMAD::EvalStatusType::EVAL_OK == eval->getEvalStatus())
        {
            NOMAD::Double h = eval->getH(computeType);
            if (!h.isDefined() || 0.0 != h)
            {
                std::string err = "Error: Barrier: xFeas' h value must be 0.0, got: " + h.display();
                throw NOMAD::Exception(__FILE__,__LINE__,err);
            }
        }
    }
}
