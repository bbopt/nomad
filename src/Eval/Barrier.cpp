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

#include "../Cache/CacheBase.hpp"
#include "../Eval/Barrier.hpp"
#include "../Eval/ComputeSuccessType.hpp"
#include "../Output/OutputQueue.hpp"

void NOMAD::Barrier::init(const NOMAD::Point& fixedVariable,
                          const NOMAD::EvalType& evalType,
                          const std::vector<NOMAD::EvalPoint>& evalPointList,
                          const NOMAD::ComputeType& computeType,
                          bool barrierInitializedFromCache)
{
    std::vector<NOMAD::EvalPoint> cachePoints;

    if (fixedVariable.isEmpty())
    {
        std::string s = "Error: Fixed variable of dimension 0";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    if (barrierInitializedFromCache)
    {
        checkCache();

        // Get best feasible and infeasible solutions from cache.
        // Points from cache are in full dimension. Convert them
        // to subproblem dimension.
        if (NOMAD::CacheBase::getInstance()->findBestFeas(cachePoints, fixedVariable, evalType, computeType, nullptr) > 0)
        {
            for (auto evalPoint : cachePoints)
            {
                NOMAD::EvalPoint evalPointSub = evalPoint.makeSubSpacePointFromFixed(fixedVariable);
                _xFeas.push_back(evalPointSub);
            }
            cachePoints.clear();
        }
        if (NOMAD::CacheBase::getInstance()->findBestInf(cachePoints, _hMax, fixedVariable, evalType, computeType, nullptr) > 0)
        {
            for (auto evalPoint : cachePoints)
            {
                NOMAD::EvalPoint evalPointSub = evalPoint.makeSubSpacePointFromFixed(fixedVariable);
                _xInf.push_back(evalPointSub);
            }
            cachePoints.clear();
        }
    }

    // Constructor's call to update should not update ref best points.
    updateWithPoints(evalPointList, evalType, computeType, true);    // true: keep all points

    setN();


    // Check: xFeas or xInf could be non-evaluated, but not both.
    auto xFeas = getFirstXFeas();
    auto xInf = getFirstXInf();
    if (   (nullptr == xFeas || nullptr == xFeas->getEval(evalType))
        && (nullptr == xInf  || nullptr == xInf->getEval(evalType)))
    {
        std::string s = "Barrier constructor: xFeas or xInf must be evaluated.\n";
        if (nullptr != xFeas)
        {
            s += "There are " + std::to_string(_xFeas.size()) + " xFeas, the first one is:\n";
            s += xFeas->displayAll();
        }
        if (nullptr != xInf)
        {
            s += "There are " + std::to_string(_xInf.size()) + " xInf, the first one is:\n";
            s += xInf->displayAll();
        }
        if (_xFeas.size() == 0 && _xInf.size() == 0)
        {
            s += "There are no xFeas and no xInf defined.";
        }
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    checkHMax();
}


void NOMAD::Barrier::setN()
{
    bool isSet = (0 != _n);
    std::string s;

    for (auto evalPoint : getAllPoints())
    {
        if (!isSet)
        {
            _n = evalPoint.size();
            isSet = true;
        }
        else if (evalPoint.size() != _n)
        {
            s = "Barrier has points of size " + std::to_string(_n) + " and of size ";
            s += std::to_string(evalPoint.size());
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
void NOMAD::Barrier::checkCache()
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


void NOMAD::Barrier::checkXFeas(const NOMAD::EvalPoint &xFeas,
                                const NOMAD::EvalType& evalType,
                                const NOMAD::ComputeType& computeType)
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


void NOMAD::Barrier::checkXFeasIsFeas(const NOMAD::EvalPoint &xFeas,
                                      const NOMAD::EvalType& evalType,
                                      const NOMAD::ComputeType& computeType)
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


void NOMAD::Barrier::checkXInf(const NOMAD::EvalPoint &xInf, const NOMAD::EvalType& evalType)
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


void NOMAD::Barrier::checkHMax()
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


void NOMAD::Barrier::updateRefBests()
{
    _refBestFeas = getFirstXFeas();
    _refBestInf  = getFirstXInf();
}


NOMAD::EvalPointPtr NOMAD::Barrier::getFirstXFeas() const
{
    NOMAD::EvalPointPtr xFeas = nullptr;
    if (_xFeas.size() > 0)
    {
        xFeas = std::make_shared<NOMAD::EvalPoint>(_xFeas[0]);
    }

    return xFeas;
}


void NOMAD::Barrier::addXFeas(const NOMAD::EvalPoint &xFeas,
                              const NOMAD::EvalType& evalType,
                              const NOMAD::ComputeType& computeType)
{
    checkXFeas(xFeas, evalType, computeType);
    _xFeas.push_back(xFeas);
}


void NOMAD::Barrier::clearXFeas()
{
    _xFeas.clear();
}


NOMAD::EvalPointPtr NOMAD::Barrier::getFirstXInf() const
{
    NOMAD::EvalPointPtr xInf = nullptr;
    if (_xInf.size() > 0)
    {
        xInf = std::make_shared<NOMAD::EvalPoint>(_xInf[0]);
    }

    return xInf;
}


void NOMAD::Barrier::addXInf(const NOMAD::EvalPoint &xInf,
                             const NOMAD::EvalType& evalType)
{
    checkXInf(xInf, evalType);
    _xInf.push_back(xInf);
}

NOMAD::SuccessType NOMAD::Barrier::getSuccessTypeOfPoints(const std::shared_ptr<EvalPoint> & xFeas,
                                                          const std::shared_ptr<EvalPoint> & xInf,
                                                          const NOMAD::EvalType & evalType,
                                                          const NOMAD::ComputeType & computeType)
{
    NOMAD::SuccessType successType = SuccessType::UNSUCCESSFUL;
    NOMAD::SuccessType successType2 = SuccessType::UNSUCCESSFUL;

    NOMAD::EvalPointPtr newBestFeas,newBestInf;

    // Get the reference best points (should work for opportunistic or not)
    auto refBestFeas = getFirstXFeas();
    auto refBestInf = getFirstXInf();

    // Set the iter success based on best feasible and best infeasible points found compared to initial point.
    if (nullptr != refBestFeas || nullptr != refBestInf)
    {
        // Compute success
        // Get which of newBestFeas and newBestInf is improving
        // the solution. Check newBestFeas first.
        NOMAD::ComputeSuccessType computeSuccess(evalType, computeType);

        if (nullptr != refBestFeas)
        {
            successType = computeSuccess(xFeas, refBestFeas);
        }
        if (nullptr != refBestInf)
        {
            successType2 = computeSuccess(xInf, refBestInf);
        }
        if (successType2 > successType)
        {
            successType = successType2;
        }
    }
    return successType;
}


// Points from evalPointList are already in subproblem dimension.
// Parameter keepAllPoints:
// If the points have the same f and h values as the ones in the barrier,
// we want to keep them all. If the parameter is false, just keep one point
// and do not update it if the new point is equivalent.
// we want to keep only the one that is already in the barrier.
bool NOMAD::Barrier::updateWithPoints(const std::vector<NOMAD::EvalPoint>& evalPointList,
                                      const NOMAD::EvalType& evalType,
                                      const NOMAD::ComputeType& computeType,
                                      const bool keepAllPoints)
{
    bool (*comp)(const NOMAD::Eval&, const NOMAD::Eval&, const NOMAD::ComputeType&, NOMAD::SuccessType, bool) = NOMAD::Eval::compInsertInBarrier;
    bool updated = false;
    bool updatedFeas = false;
    bool updatedInf = false;
    std::string s;  // for output info

    // Temporary infeasible incumbent. For insertion of more than one improving/full success
    NOMAD::EvalPoint xInfTmp;

    OUTPUT_DEBUG_START
    s = "Updating barrier with " + std::to_string(evalPointList.size()) + " suggested points";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    s = "Current barrier: ";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    s = display(4);
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    // Do separate loop on evalPointList
    // First loop update the bestFeasible.
    // If a point is a full success, set updatedFeas = true. This flag is used in the second loop.

    for (auto evalPoint : evalPointList)
    {
        // All points should be of the same size inside a Barrier.
        if (0 == _n)
        {
            //Set _n
            _n = evalPoint.size();
        }
        else if (evalPoint.size() != _n)
        {
            s = "Barrier update: Barrier points are of dimension " + std::to_string(_n);
            s += ". Trying to add this point of dimension " + std::to_string(evalPoint.size());
            s += ": " + evalPoint.display();
            throw NOMAD::Exception(__FILE__, __LINE__, s);
        }

        auto eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            OUTPUT_DEBUG_START
            if (nullptr == eval)
            {
                s = "Eval is NULL, continue";
            }
            else if (NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
            {
                s = "Eval status is " + NOMAD::enumStr(eval->getEvalStatus());
                s += ", continue";
            }
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            NOMAD::OutputQueue::Flush();
            OUTPUT_DEBUG_END

            // Suggested point is not good.
            continue;
        }

        if (eval->isFeasible(computeType))
        {
            OUTPUT_DEBUG_START
            s = "Point suggested to update barrier (feasible): " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END

            // Ensure evalPoint is as good as previous points in xFeas
            if (_xFeas.empty())
            {
                // New point is first point
                _xFeas.push_back(evalPoint);
                updatedFeas = true;
                OUTPUT_DEBUG_START
                s = "New dominating xFeas: " + evalPoint.display();
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
            }
            else if (comp(*eval, *_xFeas[0].getEval(evalType), computeType, NOMAD::SuccessType::FULL_SUCCESS,true))
            {
                // New point is better
                _xFeas.clear();
                _xFeas.push_back(evalPoint);
                updatedFeas = true;
                OUTPUT_DEBUG_START
                s = "New dominating xFeas: " + evalPoint.display();
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
            }
            else if (keepAllPoints && !comp(*_xFeas[0].getEval(evalType), *eval, computeType, NOMAD::SuccessType::FULL_SUCCESS, true))
            {
                // Points are equivalent
                // If new point is not already there, add it
                if (std::find(_xFeas.begin(), _xFeas.end(), evalPoint) == _xFeas.end())
                {
                    _xFeas.push_back(evalPoint);
                    updated = true;
                    OUTPUT_DEBUG_START
                    s = "Adding xFeas: " + evalPoint.display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                }
                else
                {
                    OUTPUT_DEBUG_START
                    s = "xFeas is already there: " + evalPoint.display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                }
            }
        }
    }
    // Do separate loop on evalPointList
    // Second loop update the bestInfeasible.
    // Use the flag oneFeasEvalFullSuccess.
    // If the flag is true hmax will not change. A point improving the best infeasible should not replace it.

    for (auto evalPoint : evalPointList)
    {
        auto eval = evalPoint.getEval(evalType);

        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            // Suggested point is not good.
            continue;
        }

        if (!eval->isFeasible())
        {
            NOMAD::Double h = eval->getH(computeType);
            if (_hMax < NOMAD::INF && (!h.isDefined() || h > _hMax))
            {
                OUTPUT_DEBUG_START
                if (h.isDefined())
                {
                    s = "H is too large: ";
                    s += h.display(NOMAD::DISPLAY_PRECISION_FULL) + " > " + _hMax.display(NOMAD::DISPLAY_PRECISION_FULL) + ", continue.";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                }
                else
                {
                    s = "H is undefined";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                }
                OUTPUT_DEBUG_END
                continue;
            }

            if (_xInf.empty())
            {
                // New point is first point
                _xInf.push_back(evalPoint);
                updatedInf = true;
                OUTPUT_DEBUG_START
                s = "New dominating xInf: " + evalPoint.display();
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
            }
            else if (!updatedFeas && comp(*eval, *_xInf[0].getEval(evalType), computeType, NOMAD::SuccessType::PARTIAL_SUCCESS,false))
            {
                if (!xInfTmp.ArrayOfDouble::isDefined() || comp(*xInfTmp.getEval(evalType),*eval,computeType, NOMAD::SuccessType::PARTIAL_SUCCESS,true))
                {
                    // Keep this eval point for comparison with other eval points
                    xInfTmp = evalPoint;

                    updatedInf = true;

                    OUTPUT_DEBUG_START
                    s = "New dominating xInf: " + evalPoint.display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                }
            }
            else if (updatedFeas && comp(*eval, *_xInf[0].getEval(evalType), computeType, NOMAD::SuccessType::FULL_SUCCESS,true))
            {
                if (!xInfTmp.ArrayOfDouble::isDefined() || comp(*xInfTmp.getEval(evalType),*eval,computeType, NOMAD::SuccessType::PARTIAL_SUCCESS,true))
                {
                    // Keep this eval point for comparison with other eval points
                    xInfTmp = evalPoint;

                    updatedInf = true;

                    OUTPUT_DEBUG_START
                    s = "New dominating xInf: " + evalPoint.display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                }
            }
        }
    }

    // Keep the new best infeasible eval point
    if (xInfTmp.ArrayOfDouble::isDefined())
    {
        _xInf.clear();
        _xInf.push_back(xInfTmp);
    }

    // Perform another pass on eval points to keep points equivalent to xInf
    if (keepAllPoints && _xInf.size() > 0)
    {
        for (auto evalPoint : evalPointList)
        {
            auto eval = evalPoint.getEval(evalType);
            if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
            {
                continue;
            }

            if (!eval->isFeasible() && !comp(*_xInf[0].getEval(evalType), *eval, computeType, NOMAD::SuccessType::PARTIAL_SUCCESS, false))
            {
                // Points are equivalent
                // If new point is not already there, add it
                if (std::find(_xInf.begin(), _xInf.end(), evalPoint) == _xInf.end())
                {
                    _xInf.push_back(evalPoint);
                    updated = true;
                    OUTPUT_DEBUG_START
                    s = "Adding xInf: " + evalPoint.display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                }
                else
                {
                    OUTPUT_DEBUG_START
                    s = "xInf is already there: " + evalPoint.display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                }
            }
        }
    }
    updated = updated || updatedFeas || updatedInf;

    OUTPUT_DEBUG_START
    if (updated)
    {
        s = "Updated barrier";
        if (updatedFeas && updatedInf)
        {
            s += " with feasible and infeasible points";
        }
        else if (updatedFeas)
        {
            s += " with feasible points";
        }
        else if (updatedInf)
        {
            s += " with infeasible points";
        }
    }
    else
    {
        s = "Barrier not updated";
    }
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    return updated;
}


void NOMAD::Barrier::clearXInf()
{
    _xInf.clear();
}


std::vector<NOMAD::EvalPoint> NOMAD::Barrier::getAllPoints() const
{
    std::vector<NOMAD::EvalPoint> allPoints;

    allPoints.reserve(_xFeas.size() + _xInf.size()); // preallocate memory
    allPoints.insert(allPoints.end(), _xFeas.begin(), _xFeas.end());
    allPoints.insert(allPoints.end(), _xInf.begin(), _xInf.end());

    return allPoints;
}


const NOMAD::EvalPoint& NOMAD::Barrier::getFirstPoint() const
{
    if (_xFeas.size() > 0)
    {
        return _xFeas[0];
    }
    else
    {
        return _xInf[0];
    }
}


void NOMAD::Barrier::setHMax(const NOMAD::Double &hMax)
{
    _hMax = hMax;
    checkHMax();
}


// Nice formatted display.
std::string NOMAD::Barrier::display(const size_t max) const
{
    std::string s;

    auto allXFeas = getAllXFeas();
    auto allXInf  = getAllXInf();
    size_t nbXFeas = 0;
    size_t nbXInf = 0;

    for (auto xFeas : allXFeas)
    {
        s += "X_FEAS " + xFeas.displayAll() + "\n";
        nbXFeas++;
        if (nbXFeas >= max && allXFeas.size() > max)
        {
            s += "... (total " + std::to_string(allXFeas.size()) + ")\n";
            break;
        }
    }
    for (auto xInf : allXInf)
    {
        s += "X_INF " + xInf.displayAll() + "\n";
        nbXInf++;

        if (nbXInf >= max && allXInf.size() > max)
        {
            s += "... (total " + std::to_string(allXInf.size()) + ")\n";
            break;
        }
    }
    s += "H_MAX " + getHMax().display(NOMAD::DISPLAY_PRECISION_FULL) + "\n";
    s += "Ref Best Feasible:   " + (_refBestFeas ? _refBestFeas->displayAll() : "NULL") + "\n";
    s += "Ref Best Infeasible: " + (_refBestInf ? _refBestInf->displayAll() : "NULL") + "\n";

    return s;
}


// Raw display.
std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::Barrier& barrier)
{
    auto allXFeas = barrier.getAllXFeas();
    auto allXInf  = barrier.getAllXInf();

    for (auto xFeas : allXFeas)
    {
        os << "X_FEAS " << xFeas << std::endl;
    }
    for (auto xInf : allXInf)
    {
        os << "X_INF " << xInf << std::endl;
    }
    os << "H_MAX " << barrier.getHMax() << std::endl;

    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::Barrier& barrier)
{
    // Set up structures to gather member info
    NOMAD::EvalPoint xFeas, xInf;
    NOMAD::Double hMax;
    barrier.clearXFeas();
    barrier.clearXInf();

    // Read line by line
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
        if ("X_FEAS" == name)
        {
            is >> xFeas;
            // Looking for xFeas in cache will ensure its f and h are updated
            NOMAD::CacheBase::getInstance()->find(xFeas, xFeas);
            // EvalType undefined: No check will be done on the feasibility
            barrier.addXFeas(xFeas, NOMAD::EvalType::UNDEFINED, NOMAD::ComputeType::STANDARD);
        }
        else if ("X_INF" == name)
        {
            is >> xInf;
            // Looking for xInf in cache will ensure its f and h are updated
            NOMAD::CacheBase::getInstance()->find(xInf, xInf);
            barrier.addXInf(xInf, NOMAD::EvalType::UNDEFINED);
        }
        else if ("H_MAX" == name)
        {
            is >> hMax;
            barrier.setHMax(hMax);
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

    return is;
}
