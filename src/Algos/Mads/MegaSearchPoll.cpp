/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/MegaSearchPoll.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::MegaSearchPoll::init()
{
    setStepType(NOMAD::StepType::MEGA_SEARCH_POLL);
    verifyParentNotNull();

    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>( _megaIterAncestor );
    if (nullptr == megaIter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"An instance of class MegaSearch must have a MadsMegaIteration among its ancestors");
    }
    
    _poll = std::make_unique<NOMAD::Poll>(this);
    _search = std::make_unique<NOMAD::Search>(this);

}


void NOMAD::MegaSearchPoll::startImp()
{
    // Generate trial points using poll and search and merge them
    generateTrialPoints();
}


bool NOMAD::MegaSearchPoll::runImp()
{
    bool foundBetter = false;
    // Ensure no max lap for MegaSearchPoll. Also reset counter before evaluation.
    NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval(NOMAD::INF_SIZE_T);
    NOMAD::EvcInterface::getEvaluatorControl()->resetLapBbEval();

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }

    return foundBetter;
}


void NOMAD::MegaSearchPoll::endImp()
{
    postProcessing();
}


// Generate new points to evaluate from Poll and Search
void NOMAD::MegaSearchPoll::generateTrialPointsImp()
{
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, true);

    NOMAD::EvalPointSet trialPoints;

    // Generate trial points for Search (all enabled search methods) and Poll.
    // Note: Search and Poll generateTrialPoints() methods both
    // take care of verifying that the generated are on mesh, and also
    // update the "PointFrom" with the frame center.
    _search->generateTrialPoints();
    auto trialPointsSearch = _search->getTrialPoints();

    _poll->generateTrialPoints();
    _poll->generateTrialPointsSecondPass();

    // Add extra points to reach a given number of trial points
    // -> First: count the points that would need eval (check cache and barrier)
    _poll->countTrialPointsThatNeedEval(this);
    _poll->generateTrialPointsExtra();
    
    auto trialPointsPoll = _poll->getTrialPoints();

    // Merge two sets and remove duplicates
    // Naive implementation. Easier to understand - I could not make std::merge,
    // std::unique or std::set_union work fine.
    // Caveat: Multiple EvalPoints copy.
    for (auto point : trialPointsSearch)
    {
        insertTrialPoint(point);
    }
    for (auto point : trialPointsPoll)
    {
        insertTrialPoint(point);
    }
    
    // Complete trial points information for sorting before eval
    completeTrialPointsInformation();
    

}

