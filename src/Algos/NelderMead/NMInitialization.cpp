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
#include "../../Algos/NelderMead/NMInitialization.hpp"
#include "../../Algos/SubproblemManager.hpp"


void NOMAD::NMInitialization::init()
{
    setStepType(NOMAD::StepType::INITIALIZATION);

    _nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get( _stopReasons );
}


bool NOMAD::NMInitialization::runImp()
{
    bool doContinue = ! _stopReasons->checkTerminate();

    if (doContinue)
    {
        // For a standalone NM, evaluate the trial points generated during start (simplex is created later)
        // Otherwise, there are no trial points available
        evalTrialPoints(this);
        doContinue = ! _stopReasons->checkTerminate();
        if ( ! doContinue )
            _nmStopReason->set(NOMAD::NMStopType::INITIAL_FAILED);

    }
    return doContinue;
}

void NOMAD::NMInitialization::startImp()
{

    if ( ! _stopReasons->checkTerminate() )
    {
        // If needed, generate trial points and put them in cache to form simplex
        // For a standalone optimization (NM_OPTIMIZATION true), initial trial points must be generated to form a valid simplex around x0. Otherwise, the cache will be used to construct the simplex.
        auto nm_opt = _runParams->getAttributeValue<bool>("NM_OPTIMIZATION");
        if ( nm_opt && ! checkCacheCanFormSimplex() )
        {
            generateTrialPoints();
        }
    }

}


void NOMAD::NMInitialization::endImp()
{
    // Construct _barrier member with evaluated _trialPoints for future use
    // _trialPoints are already updated with Evals.
    if (_trialPoints.size() > 0)
    {
        std::vector<NOMAD::EvalPoint> evalPointList;
        std::copy(_trialPoints.begin(), _trialPoints.end(),
                          std::back_inserter(evalPointList));
        auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
        _barrier = std::make_shared<NOMAD::Barrier>(hMax,
                                NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                NOMAD::EvcInterface::getEvaluatorControl()->getEvalType(),
                                NOMAD::EvcInterface::getEvaluatorControl()->getComputeType(),
                                evalPointList);
    }
}


bool NOMAD::NMInitialization::checkCacheCanFormSimplex()
{
    // Complete this function: see Issue #393
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    if ( NOMAD::CacheBase::getInstance()->size() < n+1 )
        return false;
    return false;

}

// Generate trial points to form a simplex
void NOMAD::NMInitialization::generateTrialPoints()
{
    NOMAD::Point x0 = _pbParams->getAttributeValue<NOMAD::Point>("X0");
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    if (!x0.isComplete() || x0.size() != n)
    {
        std::string err = "Initialization: evalY0: Invalid X0 " + x0.display();
        size_t cacheSize = NOMAD::CacheBase::getInstance()->size();
        if (cacheSize > (size_t)0)
        {
            err += ". Hint: Try not setting X0 so that the cache is used (";
            err += std::to_string(cacheSize) + " points)";
        }
        else
        {
            err += ". Cache is empty.";
        }
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    NOMAD::EvalPoint evalPoint_x0(x0);
    insertTrialPoint(evalPoint_x0);
    OUTPUT_INFO_START
    AddOutputInfo("Using X0: " + evalPoint_x0.display());
    OUTPUT_INFO_END

    // Method to generate simplex points using X0 adapted from fminsearch (matlab)
    const NOMAD::Double usualDelta = 0.05;    //  x0 + 5 percent
    const NOMAD::Double zeroDelta = 0.00025;  //
    for ( size_t j = 0 ; j < n ; j++ )
    {
        NOMAD::EvalPoint trialPoint(x0);
        if ( trialPoint[j] != 0 )
            trialPoint[j] *= (1 + usualDelta );
        else
            trialPoint[j] = zeroDelta;

        insertTrialPoint(trialPoint);
    }

    OUTPUT_INFO_START
    NOMAD::OutputQueue::Flush();
    OUTPUT_INFO_END
}
