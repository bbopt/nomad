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
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Algos/PSDMads/PSDMads.hpp"
#include "../../Algos/PSDMads/PSDMadsMegaIteration.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/MicroSleep.hpp"

#include <thread>

void NOMAD::PSDMads::init(const std::vector<NOMAD::EvaluatorPtr> &evaluators,
                          const std::shared_ptr<NOMAD::EvaluatorControlParameters> &evalContParams)
{
    setStepType(NOMAD::StepType::ALGORITHM_PSD_MADS);
    verifyParentNotNull();

    // This is important for nested parallel regions. A parallel region
    // for evaluation queue (for) and parallel region for sub-algos (see in runImp).
    omp_set_nested(true);
    omp_set_dynamic(false); // Not sure about this one!

    // Get the evaluator control
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    auto blockSize = evc->getEvaluatorControlGlobalParams()->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
    if (blockSize > 1)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "PSD-Mads: eval points blocks are not supported.");
    }
    
    // Check the number of parallel evaluation handled by the evaluator control
    OUTPUT_INFO_START
    if (evc->getNbThreadsForParallelEval() != 1)
    {

        std::string s = "Warning: In addition to threads used by PSD Workers and Master, several threads are used for parallel evaluations of trial points. This option has not been fully tested!";
        AddOutputInfo(s);
    }
    OUTPUT_INFO_END

    // Instantiate MadsInitialization member
    _initialization = std::make_unique<NOMAD::MadsInitialization>(this);

    // Initialize all the main threads we will need.
    // The main threads will be the ones with thread numbers 0-(nbMainThreads-1).
    // Main thread 0 is already added to EvaluatorControl at its creation.
    size_t nbMainThreads = _runParams->getAttributeValue<size_t>("PSD_MADS_NB_SUBPROBLEM");
    for (int mainThreadNum = 1; mainThreadNum < (int)nbMainThreads; mainThreadNum++)
    {
        auto subProblemEvalContParams = std::make_unique<NOMAD::EvaluatorControlParameters>(*evalContParams);
        subProblemEvalContParams->checkAndComply();

        // add a main thread
        evc->addMainThread(mainThreadNum, std::move(subProblemEvalContParams));
        // Add evaluators to the main thread
        for (const auto &ev : evaluators)
        {
            // The same evaluators are used by all mainThreads (no move(ev)).
            evc->addEvaluator(ev, mainThreadNum);
        }
        // Select BB evaluator. Otherwise, the last one added is used.
        evc->setCurrentEvaluatorType(NOMAD::EvalType::BB, mainThreadNum);
    }

    _randomPickup.reset();
}

void NOMAD::PSDMads::startImp()
{
    NOMAD::Algorithm::startImp();

    // Initializations
    _psdMainMesh = dynamic_cast<NOMAD::MadsInitialization *>(_initialization.get())->getMesh();
    _barrier = _initialization->getBarrier(); // Initialization was run in Algorithm::startImp()
    _masterMegaIteration = std::make_shared<NOMAD::MadsMegaIteration>(this, 0 /*iteration count*/, _barrier, _psdMainMesh, NOMAD::SuccessType::UNDEFINED);
}

bool NOMAD::PSDMads::runImp()
{

    bool terminateAll = false;
    size_t k = _masterMegaIteration->getK(); // The iteration counter for the master mega iteration. Should be zero!

    size_t t = _runParams->getAttributeValue<size_t>("PSD_MADS_NB_SUBPROBLEM");

    ///< Lock access to some elements when they are updated.
    omp_lock_t psdMadsLock;
    // Initialize lock.
    omp_init_lock(&psdMadsLock);

    // Parallel section
#pragma omp parallel num_threads(t) default(none) shared(k, psdMadsLock, terminateAll)
    {
        bool isPollster = (0 == NOMAD::getThreadNum()); // Pollster is thread 0.

        while (!terminateAll)
        {

            // Create a PSDMadsMegaIteration to manage the pollster worker and the regular workers.
            // Lock psd mads before accessing barrier.
            omp_set_lock(&psdMadsLock);
            auto bestEvalPoint = _barrier->getAllPoints()[0];
            auto barrier = _barrier->clone();                      // Barrier points are not transferred during clone.
            auto success = _masterMegaIteration->getSuccessType(); 

            NOMAD::Point fixedVariable(_pbParams->getAttributeValue<size_t>("DIMENSION"));
            if (!isPollster)
            {
                fixedVariable = *(bestEvalPoint.getX());
                generateSubproblem(fixedVariable);
            }
            omp_unset_lock(&psdMadsLock);

            NOMAD::PSDMadsMegaIteration psdMegaIteration(this, k, barrier,
                                                         _psdMainMesh, success,
                                                         bestEvalPoint, fixedVariable);

            psdMegaIteration.start();
            bool madsSuccessful = psdMegaIteration.run(); // One Mads (on pollster or subproblem) is run per PSDMadsMegaIteration.

            //  Update the reference barrier
            omp_set_lock(&psdMadsLock);
            if (madsSuccessful)
            {
                auto madsOnSubPb = psdMegaIteration.getMads();
                OUTPUT_INFO_START
                std::string s = madsOnSubPb->getName() + " was successful. Barrier is updated.";
                AddOutputInfo(s);
                OUTPUT_INFO_END

                _lastMadsSuccessful = !isPollster; // False if pollster. True if mads Successful and not pollster. Ignore pollster for this flag
                auto evalPointList = madsOnSubPb->getMegaIterationBarrier()->getAllPoints();
                NOMAD::convertPointListToFull(evalPointList, fixedVariable);

                _barrier->updateWithPoints(evalPointList,
                                           false /* not used by progressive barrier*/,
                                           true /* true: update incumbents and hMax*/);
            }
            omp_unset_lock(&psdMadsLock);

            psdMegaIteration.end();

            if (isPollster)
            {
                omp_set_lock(&psdMadsLock);
                if (doUpdateMesh())
                {

                    // Reset values
                    // Note: In the context of PSD-Mads, always reset RandomPickup variable.
                    _randomPickup.reset();
                    _lastMadsSuccessful = false;

                    // MegaIteration manages the mesh and barrier
                    // MegaIteration will not be run.
                    NOMAD::MadsUpdate update(_masterMegaIteration.get());

                    // Update will take care of enlarging or refining the mesh,
                    // based on the current _barrier.
                    update.start();
                    update.run();
                    update.end();

                    _psdMainMesh->checkMeshForStopping(_stopReasons);
                }
                omp_unset_lock(&psdMadsLock);

                // Update master mega iteration counter. Flush to all threads. Used by psdMadsMegaIteration.
                k++;
#pragma omp flush(k)

                // Test for termination. Flush to all threads to stop looping
                if (_termination->terminate(k))
                {
                    terminateAll = true;

                    // The pollster is done, send the signal to all that the cache must stop waiting.
                    NOMAD::CacheBase::getInstance()->setStopWaiting(true);

#pragma omp flush(terminateAll)
                }

                // Safeguard lock. Only the pollster (unique) accesses the _masterMegaIteration.
                omp_set_lock(&psdMadsLock);
                _masterMegaIteration->setK(k);
                omp_unset_lock(&psdMadsLock);
            }
        }
    }
    omp_destroy_lock(&psdMadsLock);

    return true;
}

void NOMAD::PSDMads::endImp()
{

    // PSD-Mads is done.
    _termination->start();
    _termination->run();
    _termination->end();

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
    auto coverage = _runParams->getAttributeValue<NOMAD::Double>("PSD_MADS_SUBPROBLEM_PERCENT_COVERAGE");
    coverage /= 100.0;

    int nbRemaining = 0;

    nbRemaining = (int)_randomPickup.getN();

    if (_lastMadsSuccessful && _runParams->getAttributeValue<bool>("PSD_MADS_ITER_OPPORTUNISTIC"))
    {
        doUpdate = true;
        OUTPUT_INFO_START
        AddOutputInfo("Last subproblem iteration was successful. Update mesh.");
        OUTPUT_INFO_END
    }
    // Coverage: a % of variables must have been covered.
    else if (nbRemaining < (1.0 - coverage.todouble()) * _pbParams->getAttributeValue<size_t>("DIMENSION"))
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
    if (nbVariablesInSubproblem >= fixedVariable.size())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "PSD-Mads: NB var in subproblem cannot be greater or equal to the overall dimension.");
    }

    for (size_t k = 0; k < nbVariablesInSubproblem; k++)
    {
        fixedVariable[_randomPickup.pickup()] = NOMAD::Double();
    }
}
