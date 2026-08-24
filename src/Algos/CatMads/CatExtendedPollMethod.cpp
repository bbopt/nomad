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


#include "../../Algos/CatMads/CatMads.hpp"
#include "../../Algos/CatMads/CatExtendedPollMethod.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/SimpleMads/SimpleMads.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::CatExtendedPollMethod::init()
{
    // Query the enabling parameter here
    auto catAlgo = getParentOfType<NOMAD::CatAlgoUtils*>(NOMAD::Step::_parentStep);
    if ( nullptr == catAlgo)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatExtendedPollMethod: must have a CatMads ancestor");
    }

    // For some testing, it is possible that _runParams is null or evaluator control is null
    if ( nullptr != _runParams && nullptr != NOMAD::EvcInterface::getEvaluatorControl() )
    {
        if ( _runParams->getAttributeValue<bool>("MEGA_SEARCH_POLL") )
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"CatExtendedPollMethod: does not work for mega search poll");
        }
    }

    const auto nbObj = NOMAD::Algorithm::getNbObj();
    if (nbObj > 1)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatExtendedPollMethod: does not work for multi objective");
    }

    _extendedPollTrigger = _runParams->getAttributeValue<NOMAD::Double>("CAT_EXTENDED_POLL_TRIGGER");
    if (_extendedPollTrigger <= 0 )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"CatExtendedPoll: extended poll trigger should be positive.");
    }
    
    setEnabled(true);

}


bool NOMAD::CatExtendedPollMethod::runImp()
{
    if (nullptr == _iterAncestor )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,getName() + " must have an Iteration ancestor.");
    }

    auto pollTrialPoints = static_cast<NOMAD::MadsIteration*>(_iterAncestor)->getPollTrialPoints();

    // Retrieve incumbent solutions for opportunistic stop
    auto bestFeasible = _iterAncestor->getMegaIterationBarrier()->getCurrentIncumbentFeas();
    auto bestInfeasible = _iterAncestor->getMegaIterationBarrier()->getCurrentIncumbentInf();
    NOMAD::FHComputeType computeType = NOMAD::defaultFHComputeType;
    auto hMax = _iterAncestor->getMegaIterationBarrier()->getHMax();

    bool success = false;
    for(const auto & pp: pollTrialPoints)
    {
        // Only check poll points from categorical poll
        // If the point is feasible, compare with best feasible value
        // Otherwise the point is infeasible and compare with best infeasible value
        
        if (nullptr != bestFeasible &&
            pp.isFeasible(computeType) && pp.getGenStep() == NOMAD::StepType::POLL_METHOD_CAT &&
            (pp.getF(computeType) <= bestFeasible->getF(computeType) + _extendedPollTrigger*(bestFeasible->getF(computeType)).abs()))
        {
            OUTPUT_INFO_START
            std::string s = "Starting feasible point for extended point";
            s += pp.getX()->display() + " with f=" + pp.getF(computeType).display() + " h=" + pp.getH(computeType).display() + " \n";
            s += "Extended point end. ";
            AddOutputInfo(s);
            OUTPUT_INFO_END

            success = runDouble(pp);
        }
        
        // For infeasible points, they must also be within hMax
        else if (nullptr != bestInfeasible && pp.getGenStep() == NOMAD::StepType::POLL_METHOD_CAT &&
                 !pp.isFeasible(computeType) && (pp.getF(computeType) <= bestInfeasible->getF(computeType) + _extendedPollTrigger*(bestInfeasible->getF(computeType)).abs())
                && (pp.getH(computeType) < hMax ))
        {
            OUTPUT_INFO_START
            std::string s = "Starting infeasible point for extended point";
            s += pp.getX()->display() + " with f=" + pp.getF(computeType).display() + " h=" + pp.getH(computeType).display() + " \n";
            s += "Extended point end. ";
            AddOutputInfo(s);
            OUTPUT_INFO_END
            success = runDouble(pp);
        }


        // Stop-loop if we have true, meaning that we had a full success on optim starting on a poll trial point
        if (success)
        {
            break;
        }
    }

    return success;
}

bool NOMAD::CatExtendedPollMethod::runDouble(const NOMAD::EvalPoint & pp)
{

    // Set specific evaluator control
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (evc->getCurrentEvalType() != NOMAD::EvalType::BB)
    {
        return false;
    }

    auto catMadsBarrier = evc->getBarrier();
    auto catMadshMax = catMadsBarrier->getHMax();

    // Two points produced and evaluated. No need to sort.
    evc->setEvalSortType(NOMAD::EvalSortType::RANDOM);

    // Get iteration current frame size
    // Frame size is a more robust criterion:
    //  1 frame size value -> 1 mesh size value
    // but 1 mesh size value -> several frame size value (G-Mesh)
    auto currentIterationFrameSize = _iterAncestor->getMesh()->getdeltaMeshSize();

    //
    // Set overall optim parameters
    //
    auto optPbParams = std::make_shared<NOMAD::PbParameters>(*_pbParams);
    auto optRunParams = std::make_shared<NOMAD::RunParameters>(*_runParams);

    // Reset to a regular optimization
    optRunParams->setAttributeValue("CATMADS_OPTIMIZATION", false);
    optPbParams->resetToDefaultValue("CAT_GROUP");


    // Set initial frame size.
    // Need to reset initial mesh size.
    // When undefined, initial mesh size is computed from initial frame size
    // in checkAndComply
    optPbParams->setAttributeValue("INITIAL_FRAME_SIZE", currentIterationFrameSize);
    optPbParams->resetToDefaultValue("INITIAL_MESH_SIZE");

    // No need to set min mesh size. Let's do one iteration with Double directions
    optRunParams->setAttributeValue("MAX_ITERATIONS", 1);

    // No need for anisotropic mesh changes.
    // Keep the current anisotropy in the mesh/frame
    optRunParams->setAttributeValue("ANISOTROPIC_MESH", false);

    // Just Double poll
    optRunParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    optRunParams->setAttributeValue("NM_SEARCH", false);
    optRunParams->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::DOUBLE);

    // What variables need to be fixed in optim -> categorical variables are in group 1
    auto lvg = _pbParams->getAttributeValue<NOMAD::ListOfVariableGroup>("VARIABLE_GROUP");
    auto it_front = lvg.begin(); // Cat variables are in first group. Otherwise use std::advance(it_front, xxx);
    std::set<size_t>::iterator it_ind;

    NOMAD::Point fixedVar(pp.size());
    for (it_ind = it_front->begin() ; it_ind != it_front->end(); it_ind++)
    {
        fixedVar[*it_ind] = pp[*it_ind];
    }
    optPbParams->setAttributeValue("FIXED_VARIABLE", fixedVar);

    // The overall evaluation budget stopping condition is checked by the evaluator control
    size_t budgetTotal = evc->getEvaluatorControlGlobalParams()->getAttributeValue<size_t>("MAX_BB_EVAL");
    // Remaining budget of evaluations to avoid extra evaluations by the Extended Poll
    size_t nbEvalsDone = evc->getBbEval();

    // Success of the extended poll to be returned.
    bool overallSuccess = false;

    // Set default success type. Maybe updated if better points are obtained
    setSuccessType(NOMAD::SuccessType::UNSUCCESSFUL);

    //
    // Let's loop on PPs.
    //

    // Local success with respect to pp.
    // Used to update pp when looping.
    bool localSuccess = false;

    bool iterateOnPPs = true;
    NOMAD::ArrayOfPoint x0s{*pp.getX()};

    // Starting point for the extended poll:
    NOMAD::EvalPoint ppCurrent = pp;
    while ( iterateOnPPs && nbEvalsDone < budgetTotal)
    {

        optPbParams->setAttributeValue("X0", x0s);

        // We do not want certain warnings appearing in sub-optimization.
        optPbParams->doNotShowWarnings();
        optPbParams->checkAndComply();

        // Check and comply for run params
        auto evcParams = evc->getEvaluatorControlGlobalParams();
        optRunParams->checkAndComply(evcParams, optPbParams);

        // Stop reason not used for the moment.
        // Could be used to indicate what causes the stop.
        auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

        NOMAD::Mads epMads(this, madsStopReasons, optRunParams, optPbParams, false /* false: barrier not initialized from cache */, false /* false: use pb fixed variables, true: use only local fixed variables */);

        // Run a single iteration of mads on Double direction
        epMads.start();

        // Get the barrier after start.
        auto epMadsBarrier = epMads.getInitializationBarrier();
        // Keep track of points in the barrier
        epMadsBarrier->setKeepInsertedPointsTag(true);

        epMads.run();
        epMads.end();

        // Get the tags of points inserted in barrier.
        auto allEvaluatedTrialPointsTags = epMadsBarrier->getInsertedPointsTag();

        NOMAD::FHComputeType computeType = NOMAD::defaultFHComputeType;
        bool fullSuccess = false;
        bool partialSuccess = false;
        iterateOnPPs = false;
        localSuccess = false;
        for ( auto & evTag: allEvaluatedTrialPointsTags)
        {
            // Get the eval point from cache
            auto ev = NOMAD::CacheBase::getInstance()->find(evTag);

            if (!ev.NOMAD::ArrayOfDouble::isDefined())
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"CatExtendedPollMethod: cannot get cache point with provided tag");
            }


            if (ev.dominates(ppCurrent, computeType))
            {
                ppCurrent = ev;
                localSuccess = true;

                OUTPUT_INFO_START
                std::string s = "EXTENDED POLL, current pp point: ";
                s += ppCurrent.display() ;
                AddOutputInfo(s);
                OUTPUT_INFO_END

                x0s.clear();
                x0s.push_back(NOMAD::Point(ppCurrent));
            }

            // Check success wrt best points in cat mads barrier.
            auto ep_temp = std::make_shared<NOMAD::EvalPoint>(ev);
            NOMAD::SuccessType successTypeTrialPoint = NOMAD::SuccessType::UNDEFINED;
            if (ep_temp->isFeasible(computeType))
            {
                successTypeTrialPoint = catMadsBarrier->getSuccessTypeOfPoints(ep_temp, nullptr);
            }
            else
            {
                successTypeTrialPoint = catMadsBarrier->getSuccessTypeOfPoints(nullptr, ep_temp);
            }

            if (successTypeTrialPoint == NOMAD::SuccessType::FULL_SUCCESS)
            {
                // Improvement of an incumbent
                fullSuccess = true;
            }
            else if (successTypeTrialPoint == NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                partialSuccess = true;
            }

            // All evaluated points should be stored (no re-evaluation). The cat mads barrier update will be performed later.
            insertTrialPoint(ev);

        }

        if (fullSuccess)
        {
            // Full success wrt to cat mads. Let's stop.
            setSuccessType(NOMAD::SuccessType::FULL_SUCCESS);
            // Stop extended poll
            overallSuccess = true;
        }
        else if(partialSuccess)
        {
            setSuccessType(NOMAD::SuccessType::PARTIAL_SUCCESS);
            // Continue extended poll iterations only if a local success.
            if (localSuccess)
            {
                iterateOnPPs = true;
            }
        }
        else
        {
            // Unsuccessful wrt to cat mads.
            // Continue extended poll iterations only if a local success.
            if (localSuccess)
            {
                iterateOnPPs = true;
            }
        }

        nbEvalsDone = evc->getBbEval();
    }

    // Set back to cat_sort
    evc->setEvalSortType(NOMAD::EvalSortType::CAT_SORT);

    return overallSuccess;
}
