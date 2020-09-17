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



