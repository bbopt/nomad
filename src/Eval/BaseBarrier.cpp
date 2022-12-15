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
#include "../Eval/BaseBarrier.hpp"
#include "../Eval/ComputeSuccessType.hpp"
#include "../Output/OutputQueue.hpp"

void NOMAD::BaseBarrier::init(const NOMAD::Point& fixedVariable,
                          NOMAD::EvalType  evalType,
                          NOMAD::ComputeType computeType,
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
            // Cache will return all infeasible points, even those with h=INF which should not be included in barrier
            for (const auto & evalPoint : cachePoints)
            {

                NOMAD::Double h = evalPoint.getEval(evalType)->getH(computeType);
                if (h.isDefined() && h != NOMAD::INF && (h <= _hMax))
                {
                    NOMAD::EvalPoint evalPointSub = evalPoint.makeSubSpacePointFromFixed(fixedVariable);
                    _xInf.push_back(evalPointSub);
                }
            }
            cachePoints.clear();
        }
    }
    if (_xFeas.size() > 0 || _xInf.size() > 0)
    {
        setN();
        checkHMax();
    }
}


void NOMAD::BaseBarrier::setN()
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
void NOMAD::BaseBarrier::checkCache()
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

void NOMAD::BaseBarrier::checkHMax()
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


void NOMAD::BaseBarrier::setHMax(const NOMAD::Double &hMax)
{
    _hMax = hMax;
    checkHMax();
}


void NOMAD::BaseBarrier::addXFeas(const NOMAD::EvalPoint &xFeas,
                              NOMAD::EvalType evalType,
                              NOMAD::ComputeType computeType)
{
    checkXFeas(xFeas, evalType, computeType);
    _xFeas.push_back(xFeas);
}

NOMAD::EvalPointPtr NOMAD::BaseBarrier::getFirstXFeas() const
{
    NOMAD::EvalPointPtr xFeas = nullptr;
    if (_xFeas.size() > 0)
    {
        xFeas = std::make_shared<NOMAD::EvalPoint>(_xFeas[0]);
    }

    return xFeas;
}


void NOMAD::BaseBarrier::clearXFeas()
{
    _xFeas.clear();
}


void NOMAD::BaseBarrier::addXInf(const NOMAD::EvalPoint &xInf,
                             NOMAD::EvalType evalType)
{
    checkXInf(xInf, evalType);
    _xInf.push_back(xInf);
}

NOMAD::EvalPointPtr NOMAD::BaseBarrier::getFirstXInf() const
{
    NOMAD::EvalPointPtr xInf = nullptr;
    if (_xInf.size() > 0)
    {
        xInf = std::make_shared<NOMAD::EvalPoint>(_xInf[0]);
    }

    return xInf;
}

void NOMAD::BaseBarrier::clearXInf()
{
    _xInf.clear();
}

void NOMAD::BaseBarrier::checkXInf(const NOMAD::EvalPoint &xInf, NOMAD::EvalType evalType)
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


std::vector<NOMAD::EvalPoint> NOMAD::BaseBarrier::getAllPoints() const
{
    std::vector<NOMAD::EvalPoint> allPoints;

    allPoints.reserve(_xFeas.size() + _xInf.size()); // preallocate memory
    allPoints.insert(allPoints.end(), _xFeas.begin(), _xFeas.end());
    allPoints.insert(allPoints.end(), _xInf.begin(), _xInf.end());

    return allPoints;
}

const NOMAD::EvalPoint& NOMAD::BaseBarrier::getFirstPoint() const
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
