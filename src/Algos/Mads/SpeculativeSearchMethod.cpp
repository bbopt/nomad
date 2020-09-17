/**
 \file   SpeculativeSearchMethod.cpp
 \brief  Speculative search (implementation)
 \author Christophe Tribes and Sebastien Le Digabel
 \date   2018-03-1
 */
#include "../../Algos/Mads/SpeculativeSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"

/*-------------------------------------------------------------*/
/*                     MADS speculative search                 */
/*-------------------------------------------------------------*/
/* Multiple points: i=1, ..., SPECULATIVE_SEARCH_MAX           */
/* d: direction of last success scaled to intersect the frame  */
/*  x_t = x_{k-1} + d * i                                      */
/*-------------------------------------------------------------*/

void NOMAD::SpeculativeSearchMethod::init()
{
    _name = "Speculative Search Method";

    //setComment("(SpecSearch)");

    auto enable = _runParams->getAttributeValue<bool>("SPECULATIVE_SEARCH");

    setEnabled(enable);
}


void NOMAD::SpeculativeSearchMethod::generateTrialPointsImp()
{
    bool canGenerate = true;
    std::shared_ptr<NOMAD::Point> pointFrom;

    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"SpeculativeSearchMethod: must have an iteration ancestor");
    }
    auto frameCenter = _iterAncestor->getFrameCenter();
    if (nullptr == frameCenter)
    {
        canGenerate = false;
    }
    else
    {
        // Test that the frame center has a valid generating direction
        pointFrom = frameCenter->getPointFrom(NOMAD::SubproblemManager::getSubFixedVariable(this));
        if (nullptr == pointFrom || *pointFrom == *frameCenter)
        {
            canGenerate = false;
        }
    }

    if (canGenerate)
    {
        auto dir = NOMAD::Point::vectorize(*pointFrom, *frameCenter);

        // Make the direction intersect the frame
        auto mesh = _iterAncestor->getMesh();
        if (nullptr == mesh)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"SpeculativeSearchMethod: must have a mesh");
        }
        NOMAD::ArrayOfDouble deltaFrameSize = mesh->getDeltaFrameSize();
        NOMAD::Double factor = NOMAD::INF;
        for (size_t i = 0; i < dir.size(); ++i)
        {
            if ( dir[i] != 0 )
            {
                factor = min(factor,deltaFrameSize[i]/dir[i].abs());
            }
        }

        if ( factor == NOMAD::INF )
        {
            std::string err("SpeculativeSearch: Cannot scale direction on frame");
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        OUTPUT_INFO_START
        AddOutputInfo("Direction before scaling: " + dir.display());
        OUTPUT_INFO_END
        auto nbSearches = _runParams->getAttributeValue<size_t>("SPECULATIVE_SEARCH_MAX");
        for (size_t i = 1; i <= nbSearches; i++)
        {
            auto diri = dir;
            diri *= factor * i;
            OUTPUT_INFO_START
            AddOutputInfo("Scaled direction on frame: " + diri.display());
            OUTPUT_INFO_END

            NOMAD::Point point = NOMAD::Point(*(frameCenter->getX()) + diri);

            // Insert the point
            insertTrialPoint(NOMAD::EvalPoint(point));

        }
    }
}
