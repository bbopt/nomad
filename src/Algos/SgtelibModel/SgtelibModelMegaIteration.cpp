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

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelFilterCache.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelIteration.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelMegaIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::SgtelibModelMegaIteration::init()
{
    setStepType(NOMAD::StepType::MEGA_ITERATION);
}


NOMAD::SgtelibModelMegaIteration::~SgtelibModelMegaIteration()
{
    // Clear model info from cache.
    // Very important so we don't have false info in a later MegaIteration.
    NOMAD::CacheBase::getInstance()->clearModelEval(NOMAD::getThreadNum());
}


void NOMAD::SgtelibModelMegaIteration::startImp()
{
    // Create EvalPoints and send them to EvaluatorControl
    generateTrialPoints();

    if (0 == getTrialPointsCount())
    {
        auto sgtelibModelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        sgtelibModelStopReasons->set(NOMAD::ModelStopType::NOT_ENOUGH_POINTS);
    }

}


bool NOMAD::SgtelibModelMegaIteration::runImp()
{
    // Evaluate points here for BB. Compute success.
    // return true if we have a partial or full success.
    bool foundBetter = false;
    std::string s;

    if (_stopReasons->checkTerminate())
    {
        OUTPUT_DEBUG_START
        s = getName() + ": stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
    }
    else
    {
        foundBetter = evalTrialPoints(this);
    }

    if (!foundBetter)
    {
        // If no better points found, we should terminate, otherwise we will spin.
        auto sgtelibModelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        sgtelibModelStopReasons->set(NOMAD::ModelStopType::NO_NEW_POINTS_FOUND);
    }


    return foundBetter;
}


void NOMAD::SgtelibModelMegaIteration::endImp()
{
    postProcessing();

    // Clear model info from cache.
    // Very important so we don't have false info in a later MegaIteration.
    NOMAD::CacheBase::getInstance()->clearModelEval(NOMAD::getThreadNum());
    NOMAD::MegaIteration::endImp();
}


void NOMAD::SgtelibModelMegaIteration::generateIterations()
{
    // Create a single Iteration for this MegaIteration.
    // The X0s will be set to all barrier xfeas and xinf by setupPbParameters().
    size_t k = _k;  // Main iteration counter
    // Note: NOMAD 3 uses SGTELIB_MODEL_TRIALS only.
    size_t nbIter = _runParams->getAttributeValue<size_t>("MAX_ITERATION_PER_MEGAITERATION");
    nbIter = std::min(nbIter, _runParams->getAttributeValue<size_t>("SGTELIB_MODEL_SEARCH_TRIALS"));

    for (size_t iterCount = 0; iterCount < nbIter; iterCount++)
    {
        std::shared_ptr<NOMAD::SgtelibModelIteration> iteration = std::make_shared<NOMAD::SgtelibModelIteration>(this, k);
        _iterList.push_back(iteration);
        k++;
    }

    OUTPUT_INFO_START
    AddOutputInfo(getName() + " has " + NOMAD::itos(nbIter) + " iteration" + ((nbIter > 1)? "s" : "") + ".");
    OUTPUT_INFO_END

    OUTPUT_DEBUG_START
    AddOutputDebug("Iterations generated:");
    for (size_t i = 0; i < nbIter; i++)
    {
        AddOutputDebug(_iterList[i]->getName());
    }
    OUTPUT_DEBUG_END
}


void NOMAD::SgtelibModelMegaIteration::runIterationsAndSetTrialPoints()
{
    std::string s;

    if (_iterList.empty())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "No iterations to run");
    }

    for (size_t i = 0; i < _iterList.size(); i++)
    {
        if (_stopReasons->checkTerminate())
        {
            break;
        }
        auto iteration = _iterList[i];

        if (nullptr == iteration)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Invalid shared pointer cast");
        }

        // Run iteration. The points are not BB-evaluated in the Iterations,
        // so ignore returned success type.
        iteration->start();
        iteration->run();
        iteration->end();

        // Update MegaIteration's trial points with Iteration's oracle points
        NOMAD::EvalPointSet oraclePoints = iteration->getOraclePoints();
        size_t nbInserted = 0;
        auto modelAlgo = getParentOfType<NOMAD::SgtelibModel*>();
        auto lb = modelAlgo->getExtendedLowerBound();
        auto ub = modelAlgo->getExtendedUpperBound();
        for (auto oraclePoint : oraclePoints)
        {
            // New oracle point - optimized on sgte
            // To be evaluated by blackbox
            // Add it to the list.
            // Snap to bounds, but there is no useful mesh in the context.
            auto pointFrom = modelAlgo->getX0();
            // PointFrom needs to be updated. It could have been set inside sub-algorithm,
            // but that value is moot.
            oraclePoint.setPointFrom(pointFrom, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
            if (snapPointToBoundsAndProjectOnMesh(oraclePoint, lb, ub))
            {
                bool inserted = insertTrialPoint(oraclePoint);
                if (inserted)
                {
                    nbInserted++;
                }
                OUTPUT_INFO_START
                s = "Generated point";
                s += (inserted) ? ": " : " not inserted: ";
                s += oraclePoint.display();
                AddOutputInfo(s);
                OUTPUT_INFO_END
            }
        }

        // If this iteration failed to generate new points, end it here.
        if (0 == nbInserted)
        {
            auto sgtelibModelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
            sgtelibModelStopReasons->set(NOMAD::ModelStopType::NO_NEW_POINTS_FOUND);
        }

        // Update MegaIteration's stop reason
        if (_stopReasons->checkTerminate())
        {
            OUTPUT_DEBUG_START
            s = getName() + " stop reason set to: " + _stopReasons->getStopReasonAsString();
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }

        _k++;   // Count one more iteration.

        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }

    }
}


void NOMAD::SgtelibModelMegaIteration::generateTrialPoints()
{
    generateIterations();
    runIterationsAndSetTrialPoints();
    filterCache();
}


void NOMAD::SgtelibModelMegaIteration::filterCache()
{
    // Select additonal candidates out of the cache
    int nbCandidates = _runParams->getAttributeValue<int>("SGTELIB_MODEL_SEARCH_CANDIDATES_NB");
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams();

    if (nbCandidates < 0)
    {
        // Update nbCandidates.
        // Use the largest value: Either BB_MAX_BLOCK_SIZE, or 2 * DIMENSION.
        nbCandidates = static_cast<int>(std::max(
                            evcParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE"),
                            2 * _pbParams->getAttributeValue<size_t>("DIMENSION")));
    }

    // We already have a certain number of points.
    nbCandidates -= (int)getTrialPointsCount();

    if (nbCandidates > 0)
    {
        // _trialPoints already contains points found by optimizing models.
        // Filter cache to add some more.
        auto modelAlgo = getParentOfType<NOMAD::SgtelibModel*>();
        NOMAD::SgtelibModelFilterCache filter(modelAlgo, nbCandidates);
        // _trialPoints is updated by filter
        filter.start();
        bool filterOk = filter.run();
        filter.end();

        if (!filterOk)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Filter cache failed");
        }

        auto lb = modelAlgo->getExtendedLowerBound();
        auto ub = modelAlgo->getExtendedUpperBound();
        for (auto oraclePoint : filter.getOraclePoints())
        {
            // Snap to bounds. No useful mesh in the context.
            if (snapPointToBoundsAndProjectOnMesh(oraclePoint, lb, ub))
            {
                insertTrialPoint(oraclePoint);
            }
        }

    }
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::SgtelibModelMegaIteration& megaIteration)
{
    megaIteration.display ( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::SgtelibModelMegaIteration& megaIteration)
{

    megaIteration.read( is );
    return is;

}
