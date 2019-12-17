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

#include "../Cache/CacheBase.hpp"
#include "../Eval/Barrier.hpp"

void NOMAD::Barrier::init(const NOMAD::Point& fixedVariable,
                          const NOMAD::EvalType& evalType)
{
    std::vector<NOMAD::EvalPoint> evalPointList;

    checkCache();

    // Get best feasible and infeasible solutions from cache.
    // Points from cache are in full dimension. Convert them
    // to subproblem dimension.
    if (NOMAD::CacheBase::getInstance()->findBestFeas(evalPointList, fixedVariable, evalType) > 0)
    {
        for (auto evalPoint : evalPointList)
        {
            NOMAD::EvalPoint evalPointSub = evalPoint.makeSubSpacePointFromFixed(fixedVariable);
            _xFeas.push_back(std::make_shared<NOMAD::EvalPoint>(evalPointSub));
        }
        evalPointList.clear();
    }
    if (NOMAD::CacheBase::getInstance()->findBestInf(evalPointList, _hMax, fixedVariable, evalType) > 0)
    {
        for (auto evalPoint : evalPointList)
        {
            NOMAD::EvalPoint evalPointSub = evalPoint.makeSubSpacePointFromFixed(fixedVariable);
            _xInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPointSub));
        }
        evalPointList.clear();
    }

    // Check: xFeas or xInf could be non-evaluated, but not both.
    try
    {
        checkXFeas(evalType);
    }
    catch (NOMAD::Exception &exceptionFeas)
    {
        try
        {
            checkXInf();
        }
        catch (NOMAD::Exception &exceptionInf)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Barrier: Either xFeas or xInf must be evaluated before being set.");
        }
    }

    checkHMax();
}


// Verify there is a Cache instanciated.
void NOMAD::Barrier::checkCache()
{
    try
    {
        NOMAD::CacheBase::getInstance();
    }
    catch (NOMAD::Exception &e)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Cache must be instanciated before initializing Barrier.");
    }
}


void NOMAD::Barrier::checkXFeas(const NOMAD::EvalType& evalType)
{
    bool badXFeas = false;
    for (size_t i = 0; i < _xFeas.size(); i++)
    {
        if (nullptr == _xFeas[i])
        {
            badXFeas = true;
        }
    }
    if (0 == _xFeas.size() || badXFeas)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Barrier: xFeas must be evaluated before being set.");
    }

    checkXFeasIsFeas(evalType);
}


void NOMAD::Barrier::checkXFeasIsFeas(const NOMAD::EvalType& evalType)
{
    if (NOMAD::EvalType::UNDEFINED == evalType)
    {
        // Skip this check
    }
    else
    {
        for (size_t i = 0; i < _xFeas.size(); i++)
        {
            auto eval = _xFeas[i]->getEval(evalType);
            if (nullptr != eval && 0.0 != eval->getH())
            {
                std::string warn = "Warning: Barrier: xFeas' H value will be enforced to 0.0. xFeas input value for h was " + eval->getH().tostring();
                std::cerr << warn << std::endl;
                _xFeas[i]->setH(0.0, evalType);
            }
        }
    }
}


void NOMAD::Barrier::checkXInf()
{
    bool badXInf = false;
    for (size_t i = 0; i < _xInf.size(); i++)
    {
        if (nullptr == _xInf[i])
        {
            badXInf = true;
        }
    }
    if (0 == _xInf.size() || badXInf)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "Barrier: xInf must be evaluated before being set.");
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


NOMAD::EvalPointPtr NOMAD::Barrier::getFirstXFeas() const
{
    NOMAD::EvalPointPtr xFeas(nullptr);
    if (_xFeas.size() > 0)
    {
        xFeas = _xFeas[0];
    }

    return xFeas;
}


void NOMAD::Barrier::addXFeas(const NOMAD::EvalPointPtr &xFeas,
                              const NOMAD::EvalType& evalType)
{
    _xFeas.push_back(xFeas);
    checkXFeas(evalType);
}


void NOMAD::Barrier::clearXFeas()
{
    _xFeas.clear();
}


NOMAD::EvalPointPtr NOMAD::Barrier::getFirstXInf() const
{
    NOMAD::EvalPointPtr xInf(nullptr);
    if (_xInf.size() > 0)
    {
        xInf = _xInf[0];
    }

    return xInf;
}


void NOMAD::Barrier::addXInf(const NOMAD::EvalPointPtr &xInf)
{
    _xInf.push_back(xInf);
    checkXInf();
}


void NOMAD::Barrier::clearXInf()
{
    _xInf.clear();
}


std::vector<NOMAD::EvalPointPtr> NOMAD::Barrier::getAllPoints()
{
    std::vector<NOMAD::EvalPointPtr> allPoints;
    for (auto xFeas : _xFeas)
    {
        allPoints.push_back(xFeas);
    }
    for (auto xInf : _xInf)
    {
        allPoints.push_back(xInf);
    }

    return allPoints;
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
        s += "X_FEAS " + xFeas->display() + "\n";
        nbXFeas++; 
        if (nbXFeas >= max && allXFeas.size() > max)
        {
            s += "... (total " + std::to_string(allXFeas.size()) + ")\n";
            break;
        }
    }
    for (auto xInf : allXInf)
    {
        s += "X_INF " + xInf->display() + "\n"; 
        nbXInf++; 

        if (nbXInf >= max && allXInf.size() > max)
        {
            s += "... (total " + std::to_string(allXInf.size()) + ")\n";
            break;
        }
    }
    s += "H_MAX " + getHMax().tostring() + "\n";


    return s;
}


// Raw display.
std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::Barrier& barrier)
{
    auto allXFeas = barrier.getAllXFeas();
    auto allXInf  = barrier.getAllXInf();

    for (auto xFeas : allXFeas)
    {
        if (nullptr != xFeas)
        {
            os << "X_FEAS " << *xFeas << std::endl;
        }
    }
    for (auto xInf : allXInf)
    {
        if (nullptr != xInf)
        {
            os << "X_INF " << *xInf << std::endl;
        }
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
            barrier.addXFeas(std::make_shared<NOMAD::EvalPoint>(xFeas),
                             NOMAD::EvalType::UNDEFINED);
        }
        else if ("X_INF" == name)
        {
            is >> xInf;
            // Looking for xInf in cache will ensure its f and h are updated
            NOMAD::CacheBase::getInstance()->find(xInf, xInf);
            barrier.addXInf(std::make_shared<NOMAD::EvalPoint>(xInf));
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

