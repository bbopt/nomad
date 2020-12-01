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
#include "../Eval/ComparePriority.hpp"


/*-------------------------*/
/* Class OrderByDirection  */
/*-------------------------*/
bool NOMAD::OrderByDirection::comp(NOMAD::EvalQueuePointPtr& point1,
                                   NOMAD::EvalQueuePointPtr& point2) const
{
    std::string err;
    bool lowerPriority = false;

    std::shared_ptr<NOMAD::Direction> lastSuccessfulDir1 = _lastSuccessfulDirs[point1->getThreadAlgo()];
    std::shared_ptr<NOMAD::Direction> lastSuccessfulDir2 = _lastSuccessfulDirs[point2->getThreadAlgo()];

    if (nullptr == lastSuccessfulDir1 || nullptr == lastSuccessfulDir2
        || !lastSuccessfulDir1->isDefined() || !lastSuccessfulDir2->isDefined()
        || 0 == lastSuccessfulDir1->norm() || 0 == lastSuccessfulDir2->norm())
    {
        lowerPriority = false;
    }
    else if (nullptr == point1 || nullptr == point2)
    {
        lowerPriority = false;
    }
    else if (!point1->getPointFrom())
    {
        lowerPriority = false;
    }
    else if (!point2->getPointFrom())
    {
        lowerPriority = true;
    }
    else
    {
        // General case, both point1 and point2 have points from.
        NOMAD::Direction dir1 = NOMAD::Point::vectorize(*point1->getPointFrom(), *point1);
        NOMAD::Direction dir2 = NOMAD::Point::vectorize(*point2->getPointFrom(), *point2);
        if (   lastSuccessfulDir1->size() != dir1.size()
            || lastSuccessfulDir2->size() != dir2.size())
        {
            err = "Error: Last successful direction is not of the same dimension as points";
            std::cerr << err << std::endl;
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }
        else if (0 == dir1.norm())
        {
            lowerPriority = false;
        }
        else if (0 == dir2.norm())
        {
            lowerPriority = true;
        }
        else
        {
            NOMAD::Double val1 = 1;
            NOMAD::Double val2 = 1;
            val1 = NOMAD::Direction::cos(*lastSuccessfulDir1, dir1);
            val2 = NOMAD::Direction::cos(*lastSuccessfulDir2, dir2);

            // The point farthest from lastSuccessfulDir gets lower priority.
            if (val1 < val2)
            {
                lowerPriority = true;
            }
        }
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


