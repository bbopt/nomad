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
#include "../../Algos/Mads/Search.hpp"
#include "../../Algos/Mads/QPSolverAlgoSearchMethod.hpp"
#include "../../Algos/Mads/QuadSearchMethod.hpp"
#include "../../Algos/Mads/SgtelibSearchMethod.hpp"
#include "../../Algos/Mads/SimpleLineSearchMethod.hpp"
#include "../../Algos/Mads/SpeculativeSearchMethod.hpp"
#include "../../Algos/Mads/TemplateAlgoSearchMethod.hpp"
#include "../../Algos/Mads/LHSearchMethod.hpp"
#include "../../Algos/Mads/NMSearchMethod.hpp"
#include "../../Algos/Mads/UserSearchMethod.hpp"
#include "../../Algos/Mads/VNSmartSearchMethod.hpp"
#include "../../Algos/Mads/VNSSearchMethod.hpp"
#include "../../Output/OutputQueue.hpp"

#ifdef TIME_STATS
#include "../../Util/Clock.hpp"

// Initialize static variables
// 11 search methods are available
std::vector<double> NOMAD::Search::_searchTime(11, 0.0);
std::vector<double> NOMAD::Search::_searchEvalTime(11, 0.0);
#endif // TIME_STATS

void NOMAD::Search::init()
{
    setStepType(NOMAD::StepType::SEARCH);
    verifyParentNotNull();

    auto speculativeSearch      = std::make_shared<NOMAD::SpeculativeSearchMethod>(this);
    auto simpleLineSearch       = std::make_shared<NOMAD::SimpleLineSearchMethod>(this);

    auto qpsolverSearch         = std::make_shared<NOMAD::QPSolverAlgoSearchMethod>(this);
    auto quadSearch             = std::make_shared<NOMAD::QuadSearchMethod>(this);
    auto sgtelibSearch          = std::make_shared<NOMAD::SgtelibSearchMethod>(this);
    auto lhSearch               = std::make_shared<NOMAD::LHSearchMethod>(this);
    auto nmSearch               = std::make_shared<NOMAD::NMSearchMethod>(this);
    auto vnsmartSearch          = std::make_shared<NOMAD::VNSmartAlgoSearchMethod>(this);
    auto vnsSearch              = std::make_shared<NOMAD::VNSSearchMethod>(this);
    auto templateAlgoSearch     = std::make_shared<NOMAD::TemplateAlgoSearchMethod>(this);


    // The default search methods will be executed in the same order
    // as they are inserted.
    _searchMethods.push_back(speculativeSearch);    // 1. speculative search
    _searchMethods.push_back(simpleLineSearch);    // 1b. speculative search
    _searchMethods.push_back(qpsolverSearch);       // 5a. QP solver on quad Model
    _searchMethods.push_back(quadSearch);           // 5b. Mads solver Quad Model Search
    _searchMethods.push_back(sgtelibSearch);        // 5b. Sgtelib Model Search
    _searchMethods.push_back(vnsSearch);            // 6a. VNS Search
    _searchMethods.push_back(vnsmartSearch);        // 6b. VNSmart algo search
    _searchMethods.push_back(lhSearch);             // 7. Latin-Hypercube (LH) search
    _searchMethods.push_back(nmSearch);             // 8. NelderMead (NM) search

    _searchMethods.push_back(templateAlgoSearch); // 10. Template algo (iteration) search (order of search method is important; a new search method copied from this template should be carefully positioned in the list)

    // Extra search method added dynamically to Mads are moved to Search
    auto *mads = getParentOfType<NOMAD::Mads*>(NOMAD::Step::_parentStep);
    if (nullptr !=mads )
    {
        auto & extraSearchMethods = mads->accessExtraSearchMethods();

        for (auto & sm: extraSearchMethods)
        {
            // Need to update the parent step and the corresponding mega iter ancestor to have access to barrier
            sm.second->setParentStep(this);
            insertSearchMethod(sm.first,sm.second);
        }
    }
}


void NOMAD::Search::startImp()
{
   // Sanity check.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    // Reset the current trial stats before run is called and increment the number of calls
    // This is managed by iteration utils when using generateTrialPoint instead of the (start, run, end) sequence.
    _trialPointStats.resetCurrentStats();
    _trialPointStats.incrementNbCalls();

}


bool NOMAD::Search::runImp()
{
    bool searchSuccessful = false;
    std::string s;

    // Sanity check. The runImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    if (!isEnabled())
    {
        // Early out - Return false: No new success found.
        OUTPUT_DEBUG_START
        AddOutputDebug("Search methods are all disabled.");
        OUTPUT_DEBUG_END
        return false;
    }

    // Go through all search methods until we get a success.
    OUTPUT_DEBUG_START
    s = "Going through all search methods until we get a success";
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
    for (size_t i = 0; !searchSuccessful && i < _searchMethods.size(); i++)
    {

        // A local user stop is requested. Do not perform remaining search methods. Stop type reset is done at the end of iteration/megaiteration and algorithm.
        if (_stopReasons->testIf(NOMAD::IterStopType::USER_ITER_STOP) || _stopReasons->testIf(NOMAD::IterStopType::USER_ALGO_STOP) ||
            _stopReasons->testIf(NOMAD::EvalGlobalStopType::CUSTOM_GLOBAL_STOP)) // C.T : I think we should test NOMAD::EvalMainThreadStopType::CUSTOM_OPPORTUNISTIC_STOP as the others are not yet triggered
        {
            break;
        }


        auto searchMethod = _searchMethods[i];
        bool enabled = searchMethod->isEnabled();
        OUTPUT_DEBUG_START
        s = "Search method " + NOMAD::stepTypeToString(searchMethod->getStepType()) + (enabled ? " is enabled" : " not enabled");



        AddOutputDebug(s);
        OUTPUT_DEBUG_END

        if (!enabled)
        {
            continue;
        }


#ifdef TIME_STATS
        double searchStartTime = NOMAD::Clock::getCPUTime();
        double searchEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS

        searchMethod->start();
        searchMethod->run();
        searchMethod->end();

#ifdef TIME_STATS
        _searchTime[i] += NOMAD::Clock::getCPUTime() - searchStartTime;
        _searchEvalTime[i] += NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime() - searchEvalStartTime;
#endif // TIME_STATS

        // Search is successful only if full success type.
        searchSuccessful = (searchMethod->getSuccessType() >= NOMAD::SuccessType::FULL_SUCCESS);
        if (searchSuccessful)
        {
            // Do not go through the other search methods if a search is
            // successful.
            OUTPUT_INFO_START
            s = searchMethod->getName();
            s += " is successful. Stop reason: ";
            s += _stopReasons->getStopReasonAsString() ;

            AddOutputInfo(s);
            OUTPUT_INFO_END
            break;
        }
    }


    return searchSuccessful;
}


void NOMAD::Search::endImp()
{
    // Sanity check. The endImp function should be called only when trial points are generated and evaluated for each search method separately.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

    if (!isEnabled())
    {
        // Early out
        return;
    }

    // Update the trial stats of the parent (Mads)
    // This is directly managed by iteration utils when using generateTrialPoint instead of the (start, run, end) sequence done here.
    _trialPointStats.updateParentStats();


    // Need to reset the EvalStopReason if a sub optimization is used during Search and the max bb is reached for this sub optimization
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (evc->testIf(NOMAD::EvalMainThreadStopType::LAP_MAX_BB_EVAL_REACHED))
    {
        evc->setStopReason(NOMAD::getThreadNum(), NOMAD::EvalMainThreadStopType::STARTED);
    }

}


void NOMAD::Search::generateTrialPointsImp()
{
    // Sanity check. The generateTrialPoints function should be called only when trial points are generated for all each search methods. Evaluations are delayed.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, true);
    for (const auto& searchMethod : _searchMethods)
    {
        if (searchMethod->isEnabled())
        {
            searchMethod->generateTrialPoints();

            // Aggregation of trial points from several search methods.
            // The trial points produced by a search method are already snapped on bounds and on mesh.
            const auto& searchMethodPoints = searchMethod->getTrialPoints();
            for (const auto& point : searchMethodPoints)
            {
                // NOTE trialPoints includes points from multiple SearchMethods.
                insertTrialPoint(point);
            }
        }
    }

    // NOTE: Trial points information is completed (MODEL or SURROGATE eval) after all poll and search points are produced
}


bool NOMAD::Search::isEnabled() const
{
    return std::any_of(_searchMethods.begin(), _searchMethods.end(),
                       [](const std::shared_ptr<NOMAD::SearchMethodBase>& searchMethod) { return searchMethod->isEnabled(); });
}
