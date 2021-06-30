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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/Search.hpp"
#include "../../Algos/Mads/QuadSearchMethod.hpp"
#include "../../Algos/Mads/SgtelibSearchMethod.hpp"
#include "../../Algos/Mads/SpeculativeSearchMethod.hpp"
#include "../../Algos/Mads/LHSearchMethod.hpp"
#include "../../Algos/Mads/NMSearchMethod.hpp"
#include "../../Algos/Mads/UserSearchMethod.hpp"
#include "../../Algos/Mads/VNSSearchMethod.hpp"
#include "../../Output/OutputQueue.hpp"

#ifdef TIME_STATS
#include "../../Util/Clock.hpp"

// Initialize static variables
// 5 search methods are available
std::vector<double> NOMAD::Search::_searchTime(5, 0.0);
std::vector<double> NOMAD::Search::_searchEvalTime(5, 0.0);
#endif // TIME_STATS

void NOMAD::Search::init()
{
    setStepType(NOMAD::StepType::SEARCH);
    verifyParentNotNull();

    auto speculativeSearch      = std::make_shared<NOMAD::SpeculativeSearchMethod>(this);
    auto userSearch             = std::make_shared<NOMAD::UserSearchMethod>(this);
    auto quadSearch             = std::make_shared<NOMAD::QuadSearchMethod>(this);
    auto sgtelibSearch          = std::make_shared<NOMAD::SgtelibSearchMethod>(this);
    auto lhSearch               = std::make_shared<NOMAD::LHSearchMethod>(this);
    auto nmSearch               = std::make_shared<NOMAD::NMSearchMethod>(this);
    auto vnsSearch              = std::make_shared<NOMAD::VNSSearchMethod>(this);


    // The search methods will be executed in the same order
    // as they are inserted.
    // This is the order for NOMAD 3:
    // 1. speculative search
    // 2. user search
    // 3. trend matrix basic line search
    // 4. cache search
    // 5. Model Searches
    // 6. VNS search
    // 7. Latin-Hypercube (LH) search
    // 8. NelderMead (NM) search

    _searchMethods.push_back(speculativeSearch);    // 1. speculative search
    _searchMethods.push_back(userSearch);           // 2. user search
    _searchMethods.push_back(quadSearch);           // 5a. Quad Model Searches
    _searchMethods.push_back(sgtelibSearch);        // 5b. Model Searches
    _searchMethods.push_back(vnsSearch);
    _searchMethods.push_back(lhSearch);             // 7. Latin-Hypercube (LH) search
    _searchMethods.push_back(nmSearch);             // 8. NelderMead (NM) search
}


void NOMAD::Search::startImp()
{
   // Sanity check.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, false);

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


    NOMAD::SuccessType bestSuccessYet = NOMAD::SuccessType::NOT_EVALUATED;
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    // Go through all search methods until we get a success.
    OUTPUT_DEBUG_START
    s = "Going through all search methods until we get a success";
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
    for (size_t i = 0; !searchSuccessful && i < _searchMethods.size(); i++)
    {
        auto searchMethod = _searchMethods[i];
        bool enabled = searchMethod->isEnabled();
        OUTPUT_DEBUG_START
        s = "Search method " + searchMethod->getName() + (enabled ? " is enabled" : " not enabled");
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
        if (!enabled) { continue; }
#ifdef TIME_STATS
        double searchStartTime = NOMAD::Clock::getCPUTime();
        double searchEvalStartTime = NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime();
#endif // TIME_STATS
        searchMethod->start();
        searchMethod->run();
        success = searchMethod->getSuccessType();
        searchSuccessful = (success >= NOMAD::SuccessType::FULL_SUCCESS);
        if (success > bestSuccessYet)
        {
            bestSuccessYet = success;
        }
        searchMethod->end();
#ifdef TIME_STATS
        _searchTime[i] += NOMAD::Clock::getCPUTime() - searchStartTime;
        _searchEvalTime[i] += NOMAD::EvcInterface::getEvaluatorControl()->getEvalTime() - searchEvalStartTime;
#endif // TIME_STATS


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

    setSuccessType(bestSuccessYet);

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

    // Need to reset the EvalStopReason if a sub optimization is used during Search and the max bb is reached for this sub optimization
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (evc->testIf(NOMAD::EvalMainThreadStopType::LAP_MAX_BB_EVAL_REACHED))
    {
        evc->setStopReason(NOMAD::getThreadNum(), NOMAD::EvalMainThreadStopType::STARTED);
    }

}


void NOMAD::Search::generateTrialPoints()
{
    // Sanity check. The generateTrialPoints function should be called only when trial points are generated for all each search method. After that all trials points evaluated at once.
    verifyGenerateAllPointsBeforeEval(NOMAD_PRETTY_FUNCTION, true);

    for (auto searchMethod : _searchMethods)
    {
        if (searchMethod->isEnabled())
        {
            searchMethod->generateTrialPoints();

            // Aggregation of trial points from several search methods.
            // The trial points produced by a search method are already snapped on bounds and on mesh.
            auto searchMethodPoints = searchMethod->getTrialPoints();
            for (auto point : searchMethodPoints)
            {
                // NOTE trialPoints includes points from multiple SearchMethods.
                insertTrialPoint(point);
            }
        }
    }
}


bool NOMAD::Search::isEnabled() const
{
    return std::any_of(_searchMethods.begin(), _searchMethods.end(),
                       [](std::shared_ptr<NOMAD::SearchMethodBase> searchMethod) { return searchMethod->isEnabled(); });
}
