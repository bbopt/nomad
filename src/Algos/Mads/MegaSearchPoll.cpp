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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/MegaSearchPoll.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/Search.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::MegaSearchPoll::init()
{
    _name = "MegaSearchPoll";
    verifyParentNotNull();

    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>( _megaIterAncestor );
    if (nullptr == megaIter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"An instance of class MegaSearch must have a MadsMegaIteration among its ancestors");
    }

}


const std::shared_ptr<NOMAD::MadsIteration> NOMAD::MegaSearchPoll::getIterForPoint(const NOMAD::EvalPoint& point) const
{
    return _iterForPoint[point];
}


void NOMAD::MegaSearchPoll::startImp()
{

    // Generate trial points using poll and search and merge them
    generateTrialPoints();

}


bool NOMAD::MegaSearchPoll::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }

    // Update MegaIteration success type with best success found.
    _megaIterAncestor->setSuccessType(_success);

    return foundBetter;
}


void NOMAD::MegaSearchPoll::endImp()
{
    postProcessing(NOMAD::EvcInterface::getEvaluatorControl()->getEvalType());
}


// Generate new points to evaluate from Poll and Search
void NOMAD::MegaSearchPoll::generateTrialPoints()
{
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, true);
    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + _name, true, false);
    OUTPUT_INFO_END

    NOMAD::EvalPointSet trialPoints;

    // Generate points for all frame centers, all meshes, all Search and Poll strategies.
    for (size_t i = 0; i < _megaIterAncestor->getNbIterations() ; i++)
    {
        // downcast from Iteration to MadsIteration
        const std::shared_ptr<NOMAD::MadsIteration> & iter = std::dynamic_pointer_cast<NOMAD::MadsIteration> ( _megaIterAncestor->getIter(i));

        if ( iter == nullptr )
            throw NOMAD::Exception(__FILE__, __LINE__, "Cannot convert to MadsIteration shared pointer");

        // Generate trial points for Search (all enabled search methods) and Poll.
        // Note: Search and Poll generateTrialPoints() methods both
        // take care of verifying that the generated are on mesh, and also
        // update the "PointFrom" with the Iteration frame center.
        NOMAD::Search search(iter.get() );
        search.generateTrialPoints();
        auto trialPointsSearch = search.getTrialPoints() ;

        NOMAD::Poll poll(iter.get() );
        poll.generateTrialPoints();
        auto trialPointsPoll = poll.getTrialPoints();

        // Merge two sets and remove duplicates
        // Naive implementation. Easier to understand - I could not make std::merge,
        // std::unique or std::set_union work fine.
        // Caveat: Multiple EvalPoints copy.
        for (auto point : trialPointsSearch)
        {
            insertTrialPoint( point );
            // Remember which iteration generated these points
            auto pointIterPair = std::pair<NOMAD::EvalPoint, std::shared_ptr<NOMAD::MadsIteration>>(point, std::make_shared<NOMAD::MadsIteration>(*iter));
            _iterForPoint.insert(pointIterPair);
        }
        for (auto point : trialPointsPoll)
        {
            insertTrialPoint( point );

            // Remember which iteration generated these points
            auto pointIterPair = std::pair<NOMAD::EvalPoint, std::shared_ptr<NOMAD::MadsIteration>>(point, std::make_shared<NOMAD::MadsIteration>(*iter));
            _iterForPoint.insert(pointIterPair);
        }

    }

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + NOMAD::itos(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);
    OUTPUT_INFO_END

}
