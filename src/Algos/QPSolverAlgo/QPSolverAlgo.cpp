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


#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Util/fileutils.hpp"

// QPSolver algo specifics
#include "../../Algos/QPSolverAlgo/QPSolverAlgo.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverAlgoMegaIteration.hpp"

// QuadModel specifics
#include "../../Algos/QuadModel/QuadModelInitialization.hpp"


void NOMAD::QPSolverAlgo::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_QPSOLVER);

    bool qpsolverAlgoOpt = _runParams->getAttributeValue<bool>("QP_OPTIMIZATION"); // true if standalone

    if (!qpsolverAlgoOpt)
    {
        throw NOMAD::InvalidParameter(__FILE__,__LINE__,"QP algo is intended for standalone optimization. Set QP_OPTIMIZATION true.");
    }
    _initialization = std::make_unique<NOMAD::QuadModelInitialization>(this);

    _minTrustRegionRadius = _runParams->getAttributeValue<NOMAD::Double>("QP_OPTIMIZATION_MIN_TR_RADIUS");
}

bool NOMAD::QPSolverAlgo::runImp()
{
    _algoSuccessful = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        size_t k = 0;   // Iteration number

        // Barrier was computed by Initialization (if X0 provided).
        // Barrier is used for MegaIteration management.
        std::shared_ptr<NOMAD::BarrierBase> barrier = _initialization->getBarrier();
        if (nullptr == barrier)
        {
            // Barrier constructor automatically finds the best points in the cache.

            auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
            auto hNormType = _runParams->getAttributeValue<NOMAD::HNormType>("H_NORM");

            // Compute type for this optim
            FHComputeTypeS computeType; // Default from struct initializer
            computeType.hNormType = hNormType; // REM: No PhaseOne search for this algo!

            // Eval type for this optim
            auto evalType = NOMAD::EvcInterface::getEvaluatorControl()->getCurrentEvalType();

            // Create a single objective progressive barrier
            barrier = std::make_shared<NOMAD::ProgressiveBarrier>(hMax,
                                                       NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                                       evalType,
                                                       computeType);
        }

        // Model based trust region (TR) uses the mesh to manage the TR radius (Delta).
        // Mesh refining/coarsening mechanics is used for changing the TR radius
        std::shared_ptr<NOMAD::MeshBase> mesh = std::make_shared<NOMAD::GMesh>(_pbParams,_runParams);

        // Create a single MegaIteration: manage multiple iterations.
        NOMAD::QPSolverAlgoMegaIteration megaIteration(this, k, barrier, mesh, NOMAD::SuccessType::UNDEFINED);
        while (!_termination->terminate(megaIteration.getK()))
        {
            megaIteration.start();
            bool currentMegaIterSuccess = megaIteration.run();
            megaIteration.end();

            _algoSuccessful = _algoSuccessful || currentMegaIterSuccess;

            if ( megaIteration.hasSufficientDecrease())
            {
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

                // Check for stopping. Mesh is NOT anisotropic. Compare on a single value.
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

        _termination->start();
        _termination->run();
        _termination->end();
    }

    return _algoSuccessful;
}

void NOMAD::QPSolverAlgo::readInformationForHotRestart()
{
    if (_runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES"))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"QPSolverAlgo does not currently support hot restart.");
    }
}
