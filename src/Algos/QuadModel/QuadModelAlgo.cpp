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


#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelMegaIteration.hpp"
#include "../../Algos/QuadModel/QuadModelInitialization.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Output/OutputQueue.hpp"

#include "../../../ext/sgtelib/src/Surrogate_Factory.hpp"
//

void NOMAD::QuadModelAlgo::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_QUAD_MODEL);
    verifyParentNotNull();

    // Instantiate quad model initialization class
    _initialization = std::make_unique<NOMAD::QuadModelInitialization>(this);

    _minTrustRegionRadius = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_OPTIMIZATION_MIN_TR_RADIUS");

}


/*-------------------------*/
/*       Destructor        */
/*-------------------------*/
NOMAD::QuadModelAlgo::~QuadModelAlgo() = default;

bool NOMAD::QuadModelAlgo::runImp()
{
    bool success = false;

    size_t k = 0;   // Iteration number

    if (!_termination->terminate(k))
    {
        // Barrier constructor automatically finds the best points in the cache.
        // Barrier is used for MegaIteration management.

        auto barrier = _initialization->getBarrier();
        if (nullptr == barrier)
        {
            auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");

            // Create a single objective progressive barrier
            barrier = std::make_shared<NOMAD::ProgressiveBarrier>(hMax,
                                                       NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                                       NOMAD::EvcInterface::getEvaluatorControl()->getCurrentEvalType(),
                                                       NOMAD::EvcInterface::getEvaluatorControl()->getFHComputeTypeS());
        }

        NOMAD::SuccessType megaIterSuccessType = NOMAD::SuccessType::UNDEFINED;

        auto gMeshStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();


        // Model based trust region (TR) uses the mesh to manage the TR radius (Delta).
        // Mesh refining/coarsening mechanics is used for changing the TR radius
        std::shared_ptr<NOMAD::MeshBase> mesh = std::make_shared<NOMAD::GMesh>(_pbParams,_runParams);

        // member _megaIteration is used for hot restart (read and write)
        // Update it here.
        _refMegaIteration = std::make_shared<NOMAD::QuadModelMegaIteration>(this, k, barrier, mesh, megaIterSuccessType);

        // Create an MegaIteration: manage multiple iterations around
        // different frame centers at the same time.
        NOMAD::QuadModelMegaIteration megaIteration(this, k, barrier, mesh, megaIterSuccessType);

        while (!_termination->terminate(megaIteration.getK()))
        {
            megaIteration.start();
            bool currentMegaIterSuccess = megaIteration.run();
            megaIteration.end();

            success = success || currentMegaIterSuccess;

            if (megaIteration.hasSufficientDecrease())
            {
                // The TR radius (frame size) is updated.
                mesh->enlargeDeltaFrameSize(NOMAD::Direction());
            }
            else
            {
                // Reset some stop type that should not trigger a stop
                auto stopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(megaIteration.getAllStopReasons());
                if (stopReason->testIf(NOMAD::ModelStopType::NO_NEW_POINTS_FOUND))
                {
                    megaIteration.getAllStopReasons()->setStarted();
                }
                // Refine the mesh
                mesh->refineDeltaFrameSize();

                // Check for stopping. Mesh is not anisotropic. Compare on a single value.
                if (mesh->getDeltaFrameSize()[0] < _minTrustRegionRadius)
                {
                    auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );
                    qmsStopReason->set( NOMAD::ModelStopType::MIN_TRUST_REGION_RADIUS);
                }

            }

            if (getUserInterrupt())
            {
                hotRestartOnUserInterrupt();
            }
        }
    }

    _termination->start();
    _termination->run();
    _termination->end();

    NOMAD::OutputQueue::Flush();

    return success;
}
