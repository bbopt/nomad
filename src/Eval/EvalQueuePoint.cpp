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
#include "../Eval/EvalQueuePoint.hpp"

// Initialize static variables
NOMAD::Direction NOMAD::OrderByDirection::_lastSuccessfulDir = NOMAD::Direction();
std::function<bool(NOMAD::EvalQueuePointPtr &p1, NOMAD::EvalQueuePointPtr &p2)> NOMAD::ComparePriority::_comp = basicDefaultComp;

/*-------------------------*/
/* Class OrderByDirection  */
/*-------------------------*/
bool NOMAD::OrderByDirection::comp(NOMAD::EvalQueuePointPtr& point1,
                                   NOMAD::EvalQueuePointPtr& point2)
{
    std::string err;
    bool lowerPriority = false;

    if (!_lastSuccessfulDir.isDefined() || nullptr == point1 || nullptr == point2)
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

        NOMAD::Double val1 = 1;
        NOMAD::Double val2 = 1;
        if (   _lastSuccessfulDir.size() != dir1.size()
            || _lastSuccessfulDir.size() != dir2.size())
        {
            err = "Error: Last successful direction is not of the same dimension as points";
            std::cerr << err << std::endl;
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }
        else
        {
            if (0 == _lastSuccessfulDir.norm())
            {
                lowerPriority = false;
            }
            else
            {
                val1 = NOMAD::Direction::cos(_lastSuccessfulDir, dir1);
                val2 = NOMAD::Direction::cos(_lastSuccessfulDir, dir2);
            }
        }

        // The point farthest from _lastSuccessfulDir gets lower priority.
        if (val1 < val2)
        {
            lowerPriority = true;
        }
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
        ret = _comp(point1, point2);
    }
    catch (NOMAD::Exception &e)
    {
        std::string err = "ComparePriority: Comparison failed for point1 = ";
        err += point1->display() + ", point2 = " + point2->display();
        err += " " + std::string(e.what());
        std::cerr << err << std::endl;
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    return ret;
}


// Static function.
// Currently only compares iteration number k.
bool NOMAD::ComparePriority::basicDefaultComp(NOMAD::EvalQueuePointPtr& point1,
                                              NOMAD::EvalQueuePointPtr& point2)
{
    bool hasLowerPriority = false;

    hasLowerPriority = (point1->getK() < point2->getK());

    return hasLowerPriority;
}


// Static function.
bool NOMAD::ComparePriority::lastSuccessfulDirComp(NOMAD::EvalQueuePointPtr& point1,
                                                   NOMAD::EvalQueuePointPtr& point2)
{
    return NOMAD::OrderByDirection::comp(point1, point2);
}


/*------------------------*/
/* Class EvalQueuePoint   */
/*------------------------*/



