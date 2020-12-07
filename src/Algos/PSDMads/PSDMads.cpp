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
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Algos/PSDMads/PSDMads.hpp"
#include "../../Algos/PSDMads/PSDMadsMegaIteration.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"

#include <unistd.h> // For usleep

#ifdef _OPENMP
// Initialize static lock variable
omp_lock_t NOMAD::PSDMads::_psdMadsLock;
#endif


void NOMAD::PSDMads::init(const std::shared_ptr<NOMAD::Evaluator>& evaluator,
                          const std::shared_ptr<NOMAD::EvaluatorControlParameters>& evalContParams)
{
    _name = "PSD-MADS";
    verifyParentNotNull();

    // Instanciate MadsInitialization member
    _initialization = std::make_unique<NOMAD::MadsInitialization>(this);

    // Initialize all the main threads we will need.
    // The main threads will be the ones with thread numbers 0-(nbMainThreads-1).
    // Main thread 0 is already added to EvaluatorControl at its creation.
    size_t nbMainThreads = _runParams->getAttributeValue<size_t>("PSD_MADS_NB_SUBPROBLEM");
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    for (size_t mainThreadNum = 1; mainThreadNum < nbMainThreads; mainThreadNum++)
    {
        auto evalStopReason = std::make_shared<NOMAD::StopReason<NOMAD::EvalMainThreadStopType>>();
        auto subProblemEvalContParams = std::unique_ptr<NOMAD::EvaluatorControlParameters>(new NOMAD::EvaluatorControlParameters(*evalContParams));
        subProblemEvalContParams->checkAndComply();
        evc->addMainThread(mainThreadNum, evalStopReason, evaluator, std::move(subProblemEvalContParams));
    }

    _randomPickup.reset();

#ifdef _OPENMP
    omp_init_lock(&_psdMadsLock);
#endif
}


void NOMAD::PSDMads::destroy()
{
#ifdef _OPENMP
    omp_destroy_lock(&_psdMadsLock);
#endif
}


void NOMAD::PSDMads::waitForBarrier() const
{
    // Wait until _barrier is initialized, at the beginning of run.
    while (nullptr == _barrier)
    {
        usleep(10);
    }
}


void NOMAD::PSDMads::startImp()
{
    NOMAD::Algorithm::startImp();

    if (0 == NOMAD::getThreadNum())
    {
        // Initializations that are run in thread 0 only
        size_t k = 0;
        _psdMainMesh = dynamic_cast<NOMAD::MadsInitialization*>(_initialization.get())->getMesh();
        _barrier = _initialization->getBarrier();   // Initialization was run in Algorithm::startImp()
        _megaIteration = std::make_shared<NOMAD::MadsMegaIteration>(this, k, _barrier, _psdMainMesh, NOMAD::SuccessType::NOT_EVALUATED);
    }
    else
    {
        // Other threads must wait for thread 0 to initialize _barrier.
        waitForBarrier();
    }
}


bool NOMAD::PSDMads::runImp()
{
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    const NOMAD::EvalType evalType = evc->getEvalType();
    std::string s;

    bool isPollster = (0 == NOMAD::getThreadNum());   // Pollster is thread 0, which is always a main thread.
    size_t k = _megaIteration->getK();

    while (!_termination->terminate(k))
    {
        // Create a PSDMadsMegaIteration to manage the pollster worker and the regular workers.
        auto bestEvalPoint = (_barrier->nbXFeas() > 0) ? _barrier->getFirstXFeas() : _barrier->getFirstXInf();
        NOMAD::Point fixedVariable(_pbParams->getAttributeValue<size_t>("DIMENSION"));
        if (!isPollster)
        {
            fixedVariable = *(bestEvalPoint->getX());
            generateSubproblem(fixedVariable);
        }

        NOMAD::PSDMadsMegaIteration psdMegaIteration(this, k, _barrier,
                                        _psdMainMesh, _megaIteration->getSuccessType(),
                                        *bestEvalPoint, fixedVariable);
        psdMegaIteration.start();
        bool madsSuccessful = psdMegaIteration.run();    // One Mads (on pollster or subproblem) is run per PSDMadsMegaIteration.

        if (madsSuccessful)
        {
            auto madsOnSubPb = psdMegaIteration.getMads();
            OUTPUT_INFO_START
            s = madsOnSubPb->getName() + " was successful. Barrier is updated.";
            AddOutputInfo(s);
            OUTPUT_INFO_END

#ifdef _OPENMP
            omp_set_lock(&_psdMadsLock);
#endif
            _lastMadsSuccessful = !isPollster;  // Ignore pollster for this flag
            auto evalPointList = madsOnSubPb->getMegaIterationBarrier()->getAllPoints();
            NOMAD::convertPointListToFull(evalPointList, fixedVariable);
            _barrier->updateWithPoints(evalPointList,
                                       evalType, _runParams->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"));
#ifdef _OPENMP
            omp_unset_lock(&_psdMadsLock);
#endif
        }

        psdMegaIteration.end();

        if (isPollster)
        {
            if (doUpdateMesh())
            {
                // Reset values
                // Note: In the context of PSD-Mads, always reset RandomPickup variable.
                #pragma omp critical(randomPickup)
                _randomPickup.reset();
                _lastMadsSuccessful = false;

                // MegaIteration manages the mesh and barrier
                // MegaIteration will not be run.
                NOMAD::MadsUpdate update(_megaIteration.get());

                // Update will take care of enlarging or refining the mesh,
                // based on the current _barrier.
                update.start();
                update.run();
                update.end();

                _psdMainMesh->checkMeshForStopping(_stopReasons);
            }

            // Update iteration number
            _megaIteration->setK(k++);
        }

        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }
    }

    return true;
}


void NOMAD::PSDMads::endImp()
{
    // To be run in thread 0 only.
    #pragma omp master
    {
        // PSD-Mads is done. Ensure other threads are not stuck waiting for points to be evaluated.
        NOMAD::CacheBase::getInstance()->setStopWaiting(true);
        _termination->start();
        _termination->run();
        _termination->end();
    }

    NOMAD::Algorithm::endImp();
}


// Update mesh:
// - If PSD_MADS_ORIGINAL is true, like this is done in NOMAD 3.
// - If a percentage of variables been covered
// - If the last mads on a subproblem was successful.
//
// - Not implemented: if a certain number of mads has been executed.
bool NOMAD::PSDMads::doUpdateMesh() const
{
    if (_runParams->getAttributeValue<bool>("PSD_MADS_ORIGINAL"))
    {
        // Behave like the NOMAD 3 version of PSD-Mads, i.e, always update mesh.
        return true;
    }

    bool doUpdate = false;
    auto coverage = _runParams->getAttributeValue<NOMAD::Double>("PSD_MADS_SUBPROBLEM_PCT_COVERAGE");
    coverage /= 100.0;

    int nbRemaining = 0;
    #pragma omp critical(randomPickup)
    nbRemaining = _randomPickup.getN();

    if (_lastMadsSuccessful && _runParams->getAttributeValue<bool>("PSD_MADS_ITER_OPPORTUNISTIC"))
    {
        doUpdate = true;
        OUTPUT_INFO_START
        AddOutputInfo("Last subproblem iteration was successful. Update mesh.");
        OUTPUT_INFO_END
    }
    // Coverage % of variables must have been covered.
    else if (nbRemaining < (1.0-coverage.todouble()) * _pbParams->getAttributeValue<size_t>("DIMENSION"))
    {
        OUTPUT_INFO_START
        NOMAD::Double pctCov = 100.0 - ((100.0 * nbRemaining) / _pbParams->getAttributeValue<size_t>("DIMENSION"));
        std::string s = pctCov.tostring() + "% of variables were covered by subproblems. Update mesh.";
        AddOutputInfo(s);
        OUTPUT_INFO_END
        doUpdate = true;
    }

    return doUpdate;
}


void NOMAD::PSDMads::generateSubproblem(NOMAD::Point &fixedVariable)
{
    // Setup fixedVariable to define subproblem.

    // The fixed variables of the subproblem are set to the value of best point, the remaining variables are undefined.
    const auto nbVariablesInSubproblem = _runParams->getAttributeValue<size_t>("PSD_MADS_NB_VAR_IN_SUBPROBLEM");
    #pragma omp critical(randomPickup)
    {
        for (size_t k = 0; k < nbVariablesInSubproblem; k++)
        {
            fixedVariable[_randomPickup.pickup()] = NOMAD::Double();
        }
    }
}
