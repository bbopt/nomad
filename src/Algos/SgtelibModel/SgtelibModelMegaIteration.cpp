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

#include <sstream>

#include "../../Algos/SgtelibModel/SgtelibModelFilterCache.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelIteration.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelMegaIteration.hpp"

#include "../../Algos/EvcInterface.hpp"


void NOMAD::SgtelibModelMegaIteration::init()
{
    _name = getAlgoName() + NOMAD::MegaIteration::getName();
}


NOMAD::SgtelibModelMegaIteration::~SgtelibModelMegaIteration()
{
    // Clear sgte info from cache.
    // Very important so we don't have false info in a later MegaIteration.
    NOMAD::CacheBase::getInstance()->clearSgte();
}


void NOMAD::SgtelibModelMegaIteration::startImp()
{
    // Create EvalPoints and send them to EvaluatorControl
    generateTrialPoints();

    if (0 == getTrialPointsCount())
    {
        auto sgteStopReasons = NOMAD::AlgoStopReasons<NOMAD::SgtelibModelStopType>::get(_stopReasons);
        sgteStopReasons->set(NOMAD::SgtelibModelStopType::NO_POINTS);
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
        s = getName() + ": stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
    }
    else
    {
        // DEBUG - ensure OPPORTUNISM is off.
        NOMAD::EvcInterface evcInterface(this);
        auto evcParams = evcInterface.getEvaluatorControl()->getEvaluatorControlParams();
        auto previousOpportunism = evcParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL");
        if (previousOpportunism)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Parameter OPPORTUNISTIC_EVAL should be false");
        }

        foundBetter = evalTrialPoints(this);
    }

    if (!foundBetter)
    {
        // If no better points found, we should terminate, otherwise we will spin.
        auto sgteStopReasons = NOMAD::AlgoStopReasons<NOMAD::SgtelibModelStopType>::get(_stopReasons);
        sgteStopReasons->set(NOMAD::SgtelibModelStopType::NO_NEW_POINTS_FOUND);
    }


    return foundBetter;
}


void NOMAD::SgtelibModelMegaIteration::endImp()
{
    // Clear sgte info from cache.
    // Very important so we don't have false info in a later MegaIteration.
    NOMAD::CacheBase::getInstance()->clearSgte();
    NOMAD::MegaIteration::endImp();
}


void NOMAD::SgtelibModelMegaIteration::generateIterations()
{
    // Create a single Iteration for this MegaIteration.
    // The X0s will be set to all barrier xfeas and xinf by setupPbParameters().
    size_t k = _k;  // Main iteration counter
    // Note: NOMAD 3 uses SGTELIB_MODEL_TRIALS only.
    size_t nbIter = _runParams->getAttributeValue<size_t>("MAX_ITERATION_PER_MEGAITERATION");
    nbIter = std::min(nbIter, _runParams->getAttributeValue<size_t>("SGTELIB_MODEL_TRIALS"));

    for (size_t iterCount = 0; iterCount < nbIter; iterCount++)
    {
        std::shared_ptr<NOMAD::SgtelibModelIteration> iteration = std::make_shared<NOMAD::SgtelibModelIteration>(this, k);
        _iterList.push_back(iteration);
        k++;
    }

    AddOutputInfo(_name + " has " + NOMAD::itos(nbIter) + " iteration" + ((nbIter > 1)? "s" : "") + ".");

    AddOutputDebug("Iterations generated:");
    for (size_t i = 0; i < nbIter; i++)
    {
        AddOutputDebug(_iterList[i]->getName());
    }
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
        // downcast from Iteration to SgtelibModelIteration
        std::shared_ptr<NOMAD::SgtelibModelIteration> iteration = std::dynamic_pointer_cast<NOMAD::SgtelibModelIteration>(_iterList[i]);
        
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
        auto modelAlgo = dynamic_cast<const NOMAD::SgtelibModel*>(getParentOfType<NOMAD::SgtelibModel*>());
        auto lb = modelAlgo->getExtendedLowerBound();
        auto ub = modelAlgo->getExtendedUpperBound();
        for (auto oraclePoint : oraclePoints)
        {
            // New oracle point - optimized on sgte
            // To be evaluated by blackbox
            // Add it to the list.
            // Snap to bounds, but there is no useful mesh in the context.
            if (snapPointToBoundsAndProjectOnMesh(oraclePoint, lb, ub, nullptr, nullptr))
            {
                bool inserted = insertTrialPoint(oraclePoint);
                if (inserted)
                {
                    nbInserted++;
                }
                s = "Generated point";
                s += (inserted) ? ": " : " not inserted: ";
                s += oraclePoint.display();
                AddOutputInfo(s);
            }
        }

        // If this iteration failed to generate new points, end it here.
        if (0 == nbInserted)
        {
            auto sgteStopReasons = NOMAD::AlgoStopReasons<NOMAD::SgtelibModelStopType>::get(_stopReasons);
            sgteStopReasons->set(NOMAD::SgtelibModelStopType::NO_NEW_POINTS_FOUND);
        }

        // Update MegaIteration's stop reason
        if (_stopReasons->checkTerminate())
        {
            s = getName() + " stop reason set to: " + _stopReasons->getStopReasonAsString();
            AddOutputDebug(s);
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
    int nbCandidates = _runParams->getAttributeValue<int>("SGTELIB_MODEL_CANDIDATES_NB");
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlParams();

    if (nbCandidates < 0)
    {
        // Update nbCandidates.
        // Use the largest value: Either BB_MAX_BLOCK_SIZE, or 2 * DIMENSION.
        nbCandidates = static_cast<int>(std::max(
                            evcParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE"),
                            2 * _pbParams->getAttributeValue<size_t>("DIMENSION")));
    }

    // We already have a certain number of points.
    nbCandidates -= getTrialPointsCount();

    if (nbCandidates > 0)
    {
        // _trialPoints already contains points found by optimizing models.
        // Filter cache to add some more.
        auto modelAlgo = dynamic_cast<const NOMAD::SgtelibModel*>(getParentOfType<NOMAD::SgtelibModel*>());
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
            if (snapPointToBoundsAndProjectOnMesh(oraclePoint, lb, ub, nullptr, nullptr))
            {
                insertTrialPoint(oraclePoint);
            }
        }

    }
}


void NOMAD::SgtelibModelMegaIteration::display(std::ostream& os) const
{
    NOMAD::MegaIteration::display(os);
}


void NOMAD::SgtelibModelMegaIteration::read(std::istream& is)
{
    // Set up structures to gather member info
    size_t k=0;
    // Read line by line
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
        if ("ITERATION_COUNT" == name)
        {
            is >> k;
        }
        else if ("BARRIER" == name)
        {
            if (nullptr != _barrier)
            {
                is >> *_barrier;
            }
            else
            {
                std::string err = "Error: Reading a Barrier onto a NULL pointer";
                std::cerr << err;
            }
        }
        else
        {
            for (size_t i = 0; i < name.size(); i++)
            {
                is.unget();
            }
            break;
        }
    }
    
    setK(k);
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
