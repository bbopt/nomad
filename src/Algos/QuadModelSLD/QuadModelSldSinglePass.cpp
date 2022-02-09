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

#include "../../Algos/QuadModelSLD/QuadModelSldSinglePass.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldOptimize.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldUpdate.hpp"
#include "../../Cache/CacheBase.hpp"


void NOMAD::QuadModelSldSinglePass::generateTrialPointsImp ()
{
    // Select the sample points to construct the model. Use a center pt and the cache
    NOMAD::QuadModelSldUpdate update(this,emptyEvalPointSet /* No trial points */);
    
    update.start();
    bool updateSuccess = update.run();
    update.end();

    // Model Update is handled in start().
    if (!_stopReasons->checkTerminate() && updateSuccess && getModel()->check() )
    {
        // Clear model value info from cache. For each pass we suppose we have a different quadratic model and MODEL value must be re-evaluated.
        NOMAD::CacheBase::getInstance()->clearModelEval(NOMAD::getThreadNum());
        
        // Optimize to generate oracle points on this model
        // Initialize optimize member - model optimizer on sgte
        NOMAD::QuadModelSldOptimize optimize (this, _pbParams );
        
        auto previousMaxOutputLevel =  NOMAD::OutputQueue::getInstance()->getMaxOutputLevel();
        NOMAD::OutputQueue::getInstance()->setMaxOutputLevel(NOMAD::OutputLevel::LEVEL_NORMAL);
        
        optimize.start();
        // No run, the trial points are evaluated somewhere else.
        optimize.end();
        
        NOMAD::OutputQueue::getInstance()->setMaxOutputLevel(previousMaxOutputLevel);
        
        auto trialPts = optimize.getTrialPoints();
        
        // Manage all trial points
        for ( const auto & pt : trialPts )
        {
            
            // Need to copy to a non const eval point
            NOMAD::EvalPoint scaledPt(pt);
            getModel()->unscale(scaledPt);
            insertTrialPoint(scaledPt);
            OUTPUT_DEBUG_START
            std::string s = "Unscaled xt: " + scaledPt.display();
            AddOutputInfo(s, OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
        }
    }

    // If everything is ok we set the stop reason.
    if (! _stopReasons->checkTerminate())
    {
        auto stopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );
        stopReason->set(NOMAD::ModelStopType::MODEL_SINGLE_PASS_COMPLETED);
    }

}
