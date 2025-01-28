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
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Algos/COOPMads/COOPMads.hpp"
#include "../../Algos/COOPMads/CacheSearchMethod.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/MicroSleep.hpp"


void NOMAD::COOPMads::init(const std::vector<NOMAD::EvaluatorPtr> & evaluators,
                          const std::shared_ptr<NOMAD::EvaluatorControlParameters>& evalContParams)
{

    // This is important for nested parallel regions. A parallel region
    // for evaluation queue (for) and parallel region for algos (see in runImp).
    omp_set_nested(true);
    omp_set_dynamic(false); // Not sure about this one!

    setStepType(NOMAD::StepType::ALGORITHM_COOP_MADS);
    verifyParentNotNull();

    auto blockSize = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams()->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
    if (blockSize > 1)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "COOP-Mads: eval points blocks are not supported.");
    }

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    
    // Check the number of parallel evaluation handled by the evaluator control
    OUTPUT_INFO_START
    if (evc->getNbThreadsForParallelEval() != 1)
    {

        std::string s = "Warning: In addition to threads used by COOP Mads, several threads are used for parallel evaluations of trial points. This option has not been fully tested!";
        AddOutputInfo(s);
    }
    OUTPUT_INFO_END
    
    // Instanciate MadsInitialization member (to comment)
    _initialization = std::make_unique<NOMAD::MadsInitialization>(this);

    // Initialize all the main threads (algo main threads) we will need.
    // Add an evaluator for each thread and a default evaluator type.
    // The main threads will be the ones with thread numbers 0-(nbMainThreads-1).
    // Main thread 0 is already added to EvaluatorControl at its creation.
    size_t nbMainThreads = _runParams->getAttributeValue<size_t>("COOP_MADS_NB_PROBLEM");
    for (int mainThreadNum = 1; mainThreadNum < (int)nbMainThreads; mainThreadNum++)
    {
        auto problemEvalContParams = std::make_unique<NOMAD::EvaluatorControlParameters>(*evalContParams);
        problemEvalContParams->checkAndComply();

        // add a main thread
        evc->addMainThread(mainThreadNum, std::move(problemEvalContParams));
        // Add evaluators to the main thread
        for (const auto & ev: evaluators )
        {
            // The same evaluators are used by all mainThreads (no move(ev)).
            evc->addEvaluator(ev,mainThreadNum);
        }
        // Select BB evaluator. Otherwise, the last one added is used.
        evc->setCurrentEvaluatorType(NOMAD::EvalType::BB, mainThreadNum);
    }
}

bool NOMAD::COOPMads::runImp()
{

    // Set parameters for problems
    auto problemPbParams = std::make_shared<NOMAD::PbParameters>(*_pbParams);
    auto problemRunParams = std::make_shared<NOMAD::RunParameters>(*_runParams);
    problemPbParams->checkAndComply();

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    auto evcParams = evc->getEvaluatorControlGlobalParams();
    problemRunParams->checkAndComply(evcParams, problemPbParams);

    size_t t = _runParams->getAttributeValue<size_t>("COOP_MADS_NB_PROBLEM");
#pragma omp parallel num_threads(t) default(none) shared(problemRunParams,problemPbParams,evc)
    {
        auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();
        NOMAD::Mads madsOnPb(this, madsStopReasons, problemRunParams, problemPbParams);


        // Special cache search method for coop-mads. Insert at first position.
        // The cache search is used to synchronize the best incumbents between
        // the mads.
        // During Search initialization, the extra search methods will be transferred from
        // Mads to the Search.
        auto cacheSearch = std::make_shared<NOMAD::CacheSearchMethod>(this);
        madsOnPb.insertSearchMethod(0,cacheSearch);

        int mainThreadNum = NOMAD::getThreadNum();

        OUTPUT_INFO_START
        std::string s = "Running " + madsOnPb.getName();
        s += " on thread " + NOMAD::itos(mainThreadNum);
        AddOutputInfo(s);
        OUTPUT_INFO_END

        madsOnPb.start();
        bool madsSuccessful = madsOnPb.run();
        madsOnPb.end();

        OUTPUT_INFO_START
        std::string s = "Done running " + madsOnPb.getName();
        s += " on thread " + NOMAD::itos(mainThreadNum) + ". ";
        s += "Number of evaluations: " + NOMAD::itos(evc->getBbEvalInSubproblem()) + ". "; // Modify when ready for getBbEvalInProblem()
        if (cacheSearch->isEnabled())
        {
            s += "Number of cache search trial points generated (successes): " +  NOMAD::itos(cacheSearch->getNbTrialPointsGenerated(NOMAD::EvalType::BB));

        }

        AddOutputInfo(s);
        OUTPUT_INFO_END

        // Done stop waiting in cache.
        NOMAD::CacheBase::getInstance()->setStopWaiting(true);
    }

    return true;
}

void NOMAD::COOPMads::endImp()
{
    _termination->start();
    _termination->run();
    _termination->end();
    NOMAD::Algorithm::endImp();
}
