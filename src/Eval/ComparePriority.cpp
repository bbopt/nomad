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
#include "../Eval/ComparePriority.hpp"


/*-------------------------*/
/* Class OrderByDirection  */
/*-------------------------*/
bool NOMAD::OrderByDirection::comp(NOMAD::EvalQueuePointPtr& point1,
                                   NOMAD::EvalQueuePointPtr& point2) const
{
    std::string err;
    bool lowerPriority = false; // Sorting from less interesting to most interesting point, so return true if point1 is less interesting than point2.
    bool useTag = false;    // If tie, or anything preventing computation, use tag.

    if (nullptr == point1 || nullptr == point2)
    {
        // Should not happen!
        throw NOMAD::Exception(__FILE__, __LINE__, "OrderByDirection: Comparing a point that is NULL");
    }

    std::shared_ptr<NOMAD::Direction> lastSuccessfulDir1;
    std::shared_ptr<NOMAD::Direction> lastSuccessfulDir2;

    auto point1From = point1->getPointFrom();
    auto point2From = point2->getPointFrom();

    if (nullptr != point1From && NOMAD::EvalStatusType::EVAL_OK == point1From->getEvalStatus(NOMAD::EvalType::BB))
    {
        lastSuccessfulDir1 = (point1From->isFeasible(NOMAD::EvalType::BB)) ? _lastSuccessfulFeasDirs[point1->getThreadAlgo()]
                                                                           : _lastSuccessfulInfDirs[point1->getThreadAlgo()];
    }
    if (nullptr != point2From && NOMAD::EvalStatusType::EVAL_OK == point2From->getEvalStatus(NOMAD::EvalType::BB))
    {
        lastSuccessfulDir2 = (point2From->isFeasible(NOMAD::EvalType::BB)) ? _lastSuccessfulFeasDirs[point2->getThreadAlgo()]
                                                                           : _lastSuccessfulInfDirs[point2->getThreadAlgo()];
    }

    if (   nullptr == point1From || nullptr == point2From
        || nullptr == lastSuccessfulDir1 || nullptr == lastSuccessfulDir2
        || !lastSuccessfulDir1->isComplete() || !lastSuccessfulDir2->isComplete()
        || 0 == lastSuccessfulDir1->norm() || 0 == lastSuccessfulDir2->norm())
    {
        // If no valid last direction of success revert to tag ordering
        useTag = true;
    }
    else
    {
        // General case, both point1 and point2 have points from.
        NOMAD::Direction dir1 = *point1->getDirection();
        NOMAD::Direction dir2 = *point2->getDirection();
        if (   lastSuccessfulDir1->size() != dir1.size()
            || lastSuccessfulDir2->size() != dir2.size())
        {
            err = "Error: Last successful direction is not of the same dimension as points";
            std::cerr << err << std::endl;
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }
        else if (0 == dir1.norm() || 0 == dir2.norm())
        {
            useTag = true;
        }
        else
        {
            NOMAD::Double angle1 = point1->getAngle();
            if (!angle1.isDefined())
            {
                angle1 = NOMAD::Direction::angle(*lastSuccessfulDir1, dir1);
                point1->setAngle(angle1);
            }
            NOMAD::Double angle2 = point2->getAngle();
            if (!angle2.isDefined())
            {
                angle2 = NOMAD::Direction::angle(*lastSuccessfulDir2, dir2);
                point2->setAngle(angle2);
            }
            if (!angle1.isDefined() || !angle2.isDefined())
            {
                useTag = true;
            }
            else if (angle1 < angle2)
            {
                lowerPriority = false;
            }
            else if (angle2 < angle1)
            {
                lowerPriority = true;
            }
            else
            {
                useTag = true;
            }
        }
    }

    if (useTag)
    {
        lowerPriority = (point1->getTag() > point2->getTag());
    }

    return lowerPriority;
}


/*------------------*/
/* Class RandomComp */
/*------------------*/
NOMAD::RandomComp::RandomComp(const size_t n)
  : _randomPickup(n),
    _tagToRank()
{
    setName("Random");
}


bool NOMAD::RandomComp::comp(NOMAD::EvalQueuePointPtr& point1,
                             NOMAD::EvalQueuePointPtr& point2) const
{
    size_t tag1 = point1->getTag();
    size_t tag2 = point2->getTag();

    if (_tagToRank.end() == _tagToRank.find(tag1))
    {
        _tagToRank[tag1] = _randomPickup.pickup();
    }
    if (_tagToRank.end() == _tagToRank.find(tag2))
    {
        _tagToRank[tag2] = _randomPickup.pickup();
    }

    return (_tagToRank.at(tag1) < _tagToRank.at(tag2));
}


/*-----------------*/
/* Class BasicComp */
/*-----------------*/
// Currently only compares iteration number k.
bool NOMAD::BasicComp::comp(NOMAD::EvalQueuePointPtr& point1,
                            NOMAD::EvalQueuePointPtr& point2) const
{
    return (point1->getK() < point2->getK());
}


/*------------------------*/
/* Class OrderBySurrogate */
/*------------------------*/
bool NOMAD::OrderBySurrogate::comp(NOMAD::EvalQueuePointPtr& point1,
                            NOMAD::EvalQueuePointPtr& point2) const
{
    bool useTag = false;    // If tie, or anything preventing computation, use tag.
    bool lowerPriority = false; // Sorting from less interesting to most interesting point, so return true if point1 is less interesting than point2.

    auto eval1 = point1->getEval(NOMAD::EvalType::SURROGATE);
    auto eval2 = point2->getEval(NOMAD::EvalType::SURROGATE);
    if (nullptr == eval1)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "OrderBySurrogate: SURROGATE evaluation missing for point " + point1->displayAll());
    }
    else if (nullptr == eval2)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "OrderBySurrogate: SURROGATE evaluation missing for point " + point2->displayAll());
    }

    if (eval1->dominates(*eval2))
    {
        lowerPriority = false;
    }
    else if (eval2->dominates(*eval1))
    {
        lowerPriority = true;
    }
    else
    {
        // Revert to tag ordering
        useTag = true;
    }

    if (useTag)
    {
        lowerPriority = (point1->getTag() > point2->getTag());
    }

    return lowerPriority;
}


/*------------------------*/
/* Class ComparePriority  */
/*------------------------*/
bool NOMAD::ComparePriority::operator()(NOMAD::EvalQueuePointPtr& point1,
                                        NOMAD::EvalQueuePointPtr& point2)
{
    bool ret = false;
    try
    {
        if (nullptr != _compMethod)
        {
            ret = _compMethod->comp(point1, point2);
        }
    }
    catch (NOMAD::Exception &e)
    {
        std::string compMethodName = _compMethod->getName();
        std::string err = "Error: ComparePriority: Comparison ";
        if (!compMethodName.empty())
        {
            err += "with method ";
            err += _compMethod->getName() + " ";
        }
        err += "failed for point1 = ";
        err += point1->display() + ", point2 = " + point2->display();
        err += " " + std::string(e.what());
        std::cerr << err << std::endl;
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    return ret;
}


