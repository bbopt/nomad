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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/Mads/QuadSearchMethod.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelInitialization.hpp"

void NOMAD::QuadModelInitialization::init()
{
    _name = getAlgoName() + "Initialization";
    
    _qmStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get( _stopReasons );

}


/*-------------------------*/
/*       Destructor        */
/*-------------------------*/
NOMAD::QuadModelInitialization::~QuadModelInitialization()
{
}


void NOMAD::QuadModelInitialization::startImp()
{
    
    if ( ! _stopReasons->checkTerminate() )
    {
        // For a standalone optimization (no Search Method), X0s points must be evaluated if available, otherwise, the cache can be used.
        // Do nothing if this is part of a Search Method.
        
        auto searchMethodConst = getParentOfType<NOMAD::QuadSearchMethod*>(false);

        if ( searchMethodConst == nullptr )
        {
            // The name generateTrialPoints is not well suited here because we use provided X0s and check provided cache.
            generateTrialPoints();
        }
    }
    
}


bool NOMAD::QuadModelInitialization::runImp()
{
    bool doContinue = ! _stopReasons->checkTerminate();
    
    // For a standalone optimization (no Search method), X0s points must be evaluated if available, otherwise, the cache can be used.
    // Do nothing if this is part of a sub-optimization.
    auto searchMethodConst = getParentOfType<NOMAD::QuadSearchMethod*>(false);

    if ( doContinue && searchMethodConst == nullptr )
    {
        // For a standalone quad model optimization, evaluate the X0s
        bool evalOk = eval_x0s();
        
        doContinue = ! _stopReasons->checkTerminate();
        if ( ! doContinue || ! evalOk )
            _qmStopReason->set(NOMAD::ModelStopType::X0_FAIL);

    }
    return doContinue;
}


// The name generateTrialPoints is not well suited here because we use provided X0s and check provided cache.
void NOMAD::QuadModelInitialization::generateTrialPoints()
{
    auto x0s = _pbParams->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    bool validX0available = false;
    std::string err;
    
    for (auto x0 : x0s )
    {
        if (!x0.isComplete() || x0.size() != n)
        {
            err += "Initialization: eval_x0s: Invalid X0 " + x0.display() + ".";
        }
        else
        {
            // New EvalPoint to be evaluated.
            // Add it to the list (local or in Search method).
            validX0available = insertTrialPoint(NOMAD::EvalPoint(x0));;
        }
        
    }

    if (validX0available)
    {
        if (!err.empty())
        {
            // Show invalid X0s
            AddOutputWarning(err);
        }
    }
    else
    {
        // No valid X0 available, no cache. Throw exception.
        size_t cacheSize = NOMAD::CacheBase::getInstance()->size();
        if (cacheSize > 0)
        {
            err += " Hint: Try not setting X0 so that the cache is used (";
            err += std::to_string(cacheSize) + " points)";
        }
        else
        {
            err += ". Cache is empty.";
        }
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

}


// Eval X0s, using blackbox.
// Method is copied from MadsInitialization.
bool NOMAD::QuadModelInitialization::eval_x0s()
{
    bool evalOk = false;

    // Add X0s that need evaluation to eval queue
    NOMAD::EvcInterface evcInterface(this);
    evcInterface.getEvaluatorControl()->lockQueue();

    // Enforce no opportunism.
    auto evcParams = evcInterface.getEvaluatorControl()->getEvaluatorControlParams();
    auto previousOpportunism = evcParams->getAttributeValue<bool>("OPPORTUNISTIC_EVAL");
    evcParams->setAttributeValue("OPPORTUNISTIC_EVAL", false);
    evcParams->checkAndComply();

    evcInterface.getEvaluatorControl()->unlockQueue(false); // false: do not sort eval queue

    // Evaluate all x0s. Ignore returned success type.
    // Note: EvaluatorControl would not be able to compare/compute success since there is no barrier.
    // Sanity check
    if (NOMAD::EvalType::BB != getEvalType())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Evaluation of X0 must be using blackbox");
    }
    
    evalOk = evalTrialPoints(this);
    
    // Reset opportunism to previous values.
    evcInterface.getEvaluatorControl()->lockQueue();
    evcParams->setAttributeValue("OPPORTUNISTIC_EVAL", previousOpportunism);
    evcParams->checkAndComply();
    evcInterface.getEvaluatorControl()->unlockQueue(false); // false: do not sort eval queue


    NOMAD::OutputQueue::Flush();

    return evalOk;
}
