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
#include "../../Algos/MainStep.hpp"
#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelEvaluator.hpp"
#include "../../Algos/QuadModel/QuadModelMegaIteration.hpp"
#include "../../Algos/QuadModel/QuadModelInitialization.hpp"

#include "../../Algos/QuadModel/QuadModelUpdate.hpp"

#include "../../../ext/sgtelib/src/Surrogate_Factory.hpp"
//

void NOMAD::QuadModelAlgo::init()
{
    setName("QuadModelAlgo");
    verifyParentNotNull();
   
    // Instanciate quad model initialization class
    _initialization = std::make_unique<NOMAD::QuadModelInitialization>(this);

}


/*-------------------------*/
/*       Destructor        */
/*-------------------------*/
NOMAD::QuadModelAlgo::~QuadModelAlgo()
{
}


// Start is executed when QuadModelAlgo is used as an algorithm on its own.
void NOMAD::QuadModelAlgo::startImp()
{
    // Default algorithm start. Manages initialization among other things.
    NOMAD::Algorithm::startImp();
    
    // Comment to appear at the end of stats lines
    NOMAD::MainStep::setAlgoComment("(QuadModelAlgo)");
    
}


bool NOMAD::QuadModelAlgo::runImp()
{
    bool success = false;

    NOMAD::SuccessType bestSuccessType = NOMAD::SuccessType::NOT_EVALUATED;
    
    size_t k = 0;   // Iteration number
    
    if (!_termination->terminate(k))
    {
        // Barrier constructor automatically finds the best points in the cache.
        // Barrier is used for MegaIteration management.
        
        auto barrier = _initialization->getBarrier();
        if (nullptr == barrier)
        {
            auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
            barrier = std::make_shared<NOMAD::Barrier>(hMax, getSubFixedVariable(), NOMAD::EvalType::BB);
        }
        
        NOMAD::SuccessType megaIterSuccessType = NOMAD::SuccessType::NOT_EVALUATED;
        
        // TODO fix this
        /*
        if (nullptr != _megaIteration)
        {
            // Case hot restart
            k       = _megaIteration->getK();
            barrier = _megaIteration->getBarrier();
            megaIterSuccessType = _megaIteration->getSuccessType();
        }
        */
        
        
        // A single megaiteration is done
        
        // Create an MegaIteration: manage multiple iterations around
        // different frame centers at the same time.
        NOMAD::QuadModelMegaIteration megaIteration(this, k, barrier, megaIterSuccessType);
        megaIteration.start();
        bool currentMegaIterSuccess = megaIteration.run();
        megaIteration.end();
        
        success = success || currentMegaIterSuccess;
        
        // Remember these values to construct the next MegaIteration.
        k       = megaIteration.getK();
        barrier = megaIteration.getBarrier();
        megaIterSuccessType = megaIteration.NOMAD::MegaIteration::getSuccessType();
        
        if ( megaIterSuccessType > bestSuccessType )
        {
            bestSuccessType = megaIterSuccessType;
        }
        
        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }
        
        // member _megaIteration is used for hot restart (read and write)
        // Update it here.
        _megaIteration = std::make_shared<NOMAD::QuadModelMegaIteration>(this, k++, barrier, megaIterSuccessType);
        
    }

    _termination->start();
    _termination->run();
    _termination->end();

    NOMAD::OutputQueue::Flush();
    
    return success;
}


void NOMAD::QuadModelAlgo::endImp()
{
    // Remove any remaining points from eval queue.
    EvcInterface::getEvaluatorControl()->clearQueue();

    NOMAD::MainStep::resetPreviousAlgoComment();
    NOMAD::Algorithm::endImp();
}
