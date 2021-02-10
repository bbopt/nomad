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
    std::shared_ptr<NOMAD::Point> pointFrom;

    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"SpeculativeSearchMethod: must have an iteration ancestor");
    }

    auto barrier = getMegaIterationBarrier();

    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"SpeculativeSearchMethod needs a barrier");
    }

    // Generate points starting from all points in the barrier.
    // If FRAME_CENTER_USE_CACHE is false (default), that is the same
    // as using the best feasible and best infeasible points.
    for (auto frameCenter : barrier->getAllPoints())
    {
        bool canGenerate = true;
        // Test that the frame center has a valid generating direction
        pointFrom = frameCenter.getPointFrom(NOMAD::SubproblemManager::getSubFixedVariable(this));
        if (nullptr == pointFrom || *pointFrom == frameCenter)
        {
            canGenerate = false;
        }

        if (canGenerate)
        {
            auto dir = NOMAD::Point::vectorize(*pointFrom, frameCenter);

            // Make the direction intersect the frame
            auto mesh = _iterAncestor->getMesh();
            if (nullptr == mesh)
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"SpeculativeSearchMethod: must have a mesh");
            }
            NOMAD::ArrayOfDouble deltaFrameSize = mesh->getDeltaFrameSize();
            NOMAD::ArrayOfDouble factor(dir.size(),0.0);
            for (size_t i = 0; i < dir.size(); ++i)
            {
                if ( dir[i] != 0 )
                {
                    // CT TODO determine the base factor. 3 gives better results on 18pbs with constraints than 1.
                    factor[i] = max(deltaFrameSize[i]/dir[i].abs()-1.0,3.0);
                }
            }
            OUTPUT_INFO_START
            AddOutputInfo("Direction before scaling: " + dir.display());
            OUTPUT_INFO_END

            auto nbSearches = _runParams->getAttributeValue<size_t>("SPECULATIVE_SEARCH_MAX");
            for (size_t i = 1; i <= nbSearches; i++)
            {
                auto diri = dir;
                for(size_t j = 0 ; j < dir.size(); j++)
                {
                    diri[j] *= factor[j] * i;
                }

                OUTPUT_INFO_START
                AddOutputInfo("Scaled direction : " + diri.display());
                OUTPUT_INFO_END

                // Generate
                NOMAD::Point point = NOMAD::Point(*(frameCenter.getX()) + diri);

                // Insert the point
                auto evalPoint = NOMAD::EvalPoint(point);
                evalPoint.setPointFrom(std::make_shared<NOMAD::EvalPoint>(frameCenter), NOMAD::SubproblemManager::getSubFixedVariable(this));
                evalPoint.setGenStep(getName());
                insertTrialPoint(evalPoint);
            }
        }
    }
}
