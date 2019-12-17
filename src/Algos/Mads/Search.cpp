/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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

#include "../../Algos/Mads/Search.hpp"
#include "../../Algos/Mads/SgtelibSearchMethod.hpp"
#include "../../Algos/Mads/SpeculativeSearchMethod.hpp"
#include "../../Algos/Mads/LHSearchMethod.hpp"
#include "../../Algos/Mads/NMSearchMethod.hpp"
#include "../../Algos/Mads/UserSearchMethod.hpp"

void NOMAD::Search::init()
{
    _name = "Search";
    verifyParentNotNull();

    auto lhSearch           = std::make_shared<NOMAD::LHSearchMethod>(this);
    auto userSearch         = std::make_shared<NOMAD::UserSearchMethod>(this);
    auto modelSgtelibSearch = std::make_shared<NOMAD::SgtelibSearchMethod>(this);
    auto nmSearch           = std::make_shared<NOMAD::NMSearchMethod>(this);
    auto speculativeSearch  = std::make_shared<NOMAD::SpeculativeSearchMethod>(this);

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
    _searchMethods.push_back(modelSgtelibSearch);   // 5. Model Searches
    _searchMethods.push_back(lhSearch);             // 7. Latin-Hypercube (LH) search
    _searchMethods.push_back(nmSearch);             // 8. NelderMead (NM) search
}


void NOMAD::Search::startImp()
{
   // This is sanity check.
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    if (!isEnabled())
    {
        // Early out
        return;
    }

    // Generate the points from all the enabled search methods before starting evaluations
    if ( _runParams->getAttributeValue<bool>("GENERATE_ALL_POINTS_BEFORE_EVAL") )
    {
        generateTrialPoints();
    }

}


bool NOMAD::Search::runImp()
{
    bool foundBetter = false;
    std::string s;

    // This function should be called only when trial points are generated for each search method separately and evaluated.
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    if (!isEnabled())
    {
        // Early out --> no found better!
        AddOutputDebug("Search method is disabled. Early out.");
        return false;
    }


    NOMAD::SuccessType bestSuccessYet = NOMAD::SuccessType::NOT_EVALUATED;
    NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;

    // Go through all search methods until we get a success.
    s = "Going through all search methods until we get a success";
    AddOutputDebug(s);
    for (size_t i = 0; !foundBetter && i < _searchMethods.size(); i++)
    {
        auto searchMethod = _searchMethods[i];
        bool enabled = searchMethod->isEnabled();
        s = "Search method " + searchMethod->getName() + (enabled ? " is enabled" : " not enabled");
        AddOutputDebug(s);
        if (!enabled) { continue; }
        searchMethod->start();
        foundBetter = searchMethod->run();
        success = searchMethod->getSuccessType();
        if (success > bestSuccessYet)
        {
            bestSuccessYet = success;
        }
        searchMethod->end();

        if (foundBetter)
        {
            // Do not go through the other search methods if we found
            // an improving solution.
            s = searchMethod->getName();
            s += " found an improving solution. Stop reason: ";
            s += _stopReasons->getStopReasonAsString() ;

            AddOutputInfo(s);
            break;
        }
    }

    setSuccessType(bestSuccessYet);

    return foundBetter;
}


void NOMAD::Search::endImp()
{
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    if (!isEnabled())
    {
        // Early out
        return;
    }

    // Need to reset the EvalStopReason if a sub optimization is used during Search and the max bb is reached for this sub optimization
    if (_stopReasons->testIf(NOMAD::EvalStopType::LAP_MAX_BB_EVAL_REACHED))
    {
        _stopReasons->set(NOMAD::EvalStopType::STARTED);
    }

}


// Generate trial points for ALL enabled search methods.
// To be used only when parameter GENERATE_ALL_POINTS_BEFORE_EVAL is true.
void NOMAD::Search::generateTrialPoints()
{

    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, true);

    for (auto searchMethod : _searchMethods)
    {
        if (searchMethod->isEnabled())
        {
            searchMethod->generateTrialPoints();

            // NB verifyPointsAreOnMesh and updatePointsWithFrameCenter are
            // done here, so that if an user adds a new search method, they
            // will not have to think about adding these verifications.
            searchMethod->verifyPointsAreOnMesh(getName());
            searchMethod->updatePointsWithFrameCenter();

            auto searchMethodPoints = searchMethod->getTrialPoints();

            for (auto point : searchMethodPoints)
            {
                // NOTE trialPoints includes points from multiple SearchMethods.
                // TODO Probably we should use some move operators instead of
                // copying EvalPoints.
                insertTrialPoint(point);
            }
        }
    }
}


bool NOMAD::Search::isEnabled() const
{
    bool searchEnabled = false;
    if (_searchMethods.size() > 0)
    {
        for (auto searchMethod : _searchMethods)
        {
            if (searchMethod->isEnabled())
            {
                searchEnabled = true;
                break;
            }
        }
    }

    return searchEnabled;
}
