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

#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Ads/AdsMegaIteration.hpp"
#include "../../Algos/Ads/AdsUpdate.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::AdsMegaIteration::observe(const std::vector<NOMAD::EvalPoint>& evalPointList)
{
    // Update cache with new points.
    NOMAD::EvalPoint evalPointFound;
    for (const auto& evalPoint : evalPointList)
    {
        if (NOMAD::CacheBase::getInstance()->find(evalPoint, evalPointFound))
        {
            // New eval for point already in cache
            NOMAD::CacheBase::getInstance()->update(evalPoint, NOMAD::EvalType::BB);
        }
        else
        {
            // Point is not in cache yet
            evalPoint.updateTag();
            NOMAD::CacheBase::getInstance()->smartInsert(evalPoint);
        }
    }

    // Update barrier with new points.
    _barrier->updateRefBests();
    _barrier->updateWithPoints(evalPointList, _runParams->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"), true /* true: update incumbents and hMax */);

    // Update main mesh using AdsUpdate
    NOMAD::AdsUpdate update(this);
    update.start();
    update.run();
    update.end();

    OUTPUT_DEBUG_START
    AddOutputDebug("MegaIteration generated: " + getName());
    NOMAD::ArrayOfDouble meshSize  = _mainMesh->getdeltaMeshSize();
    NOMAD::ArrayOfDouble frameSize = _mainMesh->getDeltaFrameSize();
    AddOutputDebug("Mesh size:  " + meshSize.display());
    AddOutputDebug("Frame size: " + frameSize.display());
    OUTPUT_DEBUG_END
}


void NOMAD::AdsMegaIteration::startImp()
{
    // Used for callback if user call enabled
    NOMAD::MegaIteration::startImp();
    
    // Update main mesh and barrier using AdsUpdate
    NOMAD::AdsUpdate update( this );
    update.start();
    update.run();
    update.end();

    // Verify mesh stop conditions.
    _mainMesh->checkMeshForStopping(_stopReasons);

    OUTPUT_DEBUG_START
    AddOutputDebug("Mesh Stop Reason: " + _stopReasons->getStopReasonAsString());
    OUTPUT_DEBUG_END
}
