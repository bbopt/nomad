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
#include "../../Algos/QuadModel/QuadModelSinglePass.hpp"

#include "../../Cache/CacheBase.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/DMultiMads/DMultiMadsIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMegaIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsQuadDMSSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Type/DMultiMadsSearchStrategyType.hpp"


void NOMAD::DMultiMadsQuadDMSSearchMethod::init()
{
    if (nullptr == getParentOfType<DMultiMadsIteration*>())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMadsQuadDMSSearch only works for DMultiMads");
    }

    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    _ref_compute_type = evc->getComputeType();

    if (NOMAD::ComputeType::STANDARD != _ref_compute_type)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do Quad DMS search for DMultiMads on ComputeType other than STANDARD.");
    }

    // Get barrier from upper MegaIteration, if available.
    const auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    if (megaIter == nullptr)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMadsQuadDMSSearch should have a DMultiMads mega iteration parent.");
    }
    const std::shared_ptr<NOMAD::BarrierBase> barrier = megaIter->getBarrier();
    
    if (nullptr == barrier || nullptr == std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using Quad DMS search, we need a DMultiMadsBarrier.");
    }
    
    // Number of objectives. This is used to select the objective combination to consider.
    _m = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier)->getNbObj();
    
    // Temp for testing. Cheap unit test!
    // testObjCombinations();
    
    if (evc->getCurrentEvalType() == NOMAD::EvalType::MODEL)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do Quad DMS search for DMultiMads on EvalType::MODEL.");
    }
    if (!evc->getUseCache())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using QUAD DMS search, we need a cache.");
    }
    
    setStepType(NOMAD::StepType::SEARCH_METHOD_DMULTIMADS_QUAD_DMS);

    const bool runQuadSearch = getRunParams()->getAttributeValue<bool>("QUAD_MODEL_SEARCH");
    const auto quadStrategy = getRunParams()->getAttributeValue<NOMAD::DMultiMadsQuadSearchType>("DMULTIMADS_QUAD_MODEL_STRATEGY");
    const bool isEnabled = runQuadSearch && quadStrategy == NOMAD::DMultiMadsQuadSearchType::DMS;
    setEnabled(isEnabled);
    
    if (isEnabled && std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier)->getNbObj() >= 5)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMads Quad DMS search cannot be performed for 5 or more objectives.");
    }
    
#ifndef USE_SGTELIB
    if (isEnabled())
    {
        OUTPUT_INFO_START
        AddOutputInfo(getName() + " cannot be performed because NOMAD is compiled without sgtelib library");
        OUTPUT_INFO_END
        setEnabled(false);
    }
#endif
}

void NOMAD::DMultiMadsQuadDMSSearchMethod::preRun()
{
    // Initialize these parameters for objective combinations
    _l = 0;
    _activeObjsIndex.clear();
    _indices.clear();
    
    std::shared_ptr<NOMAD::BarrierBase> barrier = nullptr;

    // Keep a reference to compute type. It will be modified in evc during quad model optim.
    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    _ref_compute_type = evc->getComputeType();

    if (NOMAD::ComputeType::STANDARD != _ref_compute_type)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do Quad DMS search for DMultiMads on ComputeType other than STANDARD.");
    }

    // Get barrier from upper MegaIteration, if available.
    const auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    if (megaIter != nullptr)
    {
        barrier = megaIter->getBarrier();
    }

    if (nullptr == barrier || nullptr == std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using Quad DMS search, we need a DMultiMadsBarrier.");
    }

    if (evc->getCurrentEvalType() == NOMAD::EvalType::MODEL)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do Quad DMS search for DMultiMads on EvalType::MODEL.");
    }
    if (!evc->getUseCache())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using Quad DMS search, we need a cache.");
    }

    // Keep a reference to the DMultiMads barrier.
    _ref_dmads_barrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier);

    // Keep Pareto elements for the detection of the Quad DMS model search success
    auto getCurrentParetoFrontApproximation = [](const DMultiMadsBarrier& barrier) -> std::vector<NOMAD::EvalPointPtr>
    {
        if (barrier.getCurrentIncumbentFeas() != nullptr)
        {
            return barrier.getAllXFeas();
        }
        else
        {
            return barrier.getAllXInf();
        }
    };
    _paretoElements.clear();
    _paretoElements = getCurrentParetoFrontApproximation(*_ref_dmads_barrier);
}


bool NOMAD::DMultiMadsQuadDMSSearchMethod::runImp()
{
    // Basic checks before running search.
    preRun();

    // Solve single-objective optimization subproblems by increasing levels
    // of objective combinations.
    size_t previousLevel = 0;
    bool successLevel = false;
    while (selectObjCombination() && !_stopReasons->checkTerminate())
    {
        // We move to the next level of combinations
        // Check if the search has reached a success for the previous level.
        if (previousLevel < _l)
        {
            // No need to continue
            if (successLevel)
                break;

            // Move to the next level
            previousLevel = _l;
        }

        generateTrialPointOnSingleObjCombination();
        
        evalTrialPoints(this);

        // From IterationUtils. Update megaIteration barrier.
        const bool successCombination = postRunUpdates();
        successLevel = std::max(successCombination, successLevel);
    }

    const bool success = (NOMAD::SuccessType::FULL_SUCCESS == getSuccessType());
    return success;
}

void NOMAD::DMultiMadsQuadDMSSearchMethod::prepareSingleObjectiveRun()
{
    // Define the single-objective function for Quad DMS run.
    NOMAD::singleOutputComputeFType singleObjCompute = [&](const BBOutputTypeList& bbOutputTypeList,
                                                        const BBOutput& bbOutput) -> NOMAD::Double
    {
        if (!bbOutput.getEvalOk() || bbOutputTypeList.empty())
        {
            return NOMAD::INF;
        }

        if (!bbOutput.checkSizeMatch(bbOutputTypeList))
        {
            return NOMAD::INF;
        }

        const auto& fvalues = bbOutput.getObjectives(bbOutputTypeList);
        NOMAD::Double maxObjFun = NOMAD::M_INF;
        for (size_t i = 0; i < fvalues.size(); i++)
        {
            if (std::any_of(_activeObjsIndex.begin(), _activeObjsIndex.end(),
                            [i](size_t j)
                            {
                                // Objective number starts at 1 in _activeObjsIndex
                                return j == i + 1;
                            }))
            {
                maxObjFun = max(maxObjFun, fvalues[i]);
            }
        }

        return maxObjFun;
    };

    // Make the objective function available to the evaluator control.
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    evc->setComputeType(NOMAD::ComputeType::DMULTI_COMBINE_F, singleObjCompute);
}


bool NOMAD::DMultiMadsQuadDMSSearchMethod::postRunUpdates()
{
    bool success = false;

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    // Set back the STANDARD compute type for DMultiMads
    evc->setComputeType(_ref_compute_type);

    const NOMAD::FHComputeType completeComputeType = {evc->getCurrentEvalType(),
                                                      evc->getFHComputeTypeS()};

    _ref_dmads_barrier->checkForFHComputeType(completeComputeType);
    NOMAD::Double hMax = _ref_dmads_barrier->getHMax();

    // Determine if the resolution of the subproblem has been a success.
    for (auto& ep: _trialPoints)
    {
        // 0- The search has generated a point that is not eval ok
        if (!ep.isEvalOk(completeComputeType.evalType))
        {
            continue;
        }
        
        
        // 1- The search has generated an infeasible point and there exists a feasible set
        // of current solutions.
        if (!ep.isFeasible(completeComputeType) && _paretoElements[0]->isFeasible(completeComputeType))
        {
            continue;
        }

        // 2- The search has generated a feasible point, and there was no feasible point
        // generated before: FULL_SUCCESS
        if ((_ref_dmads_barrier->getCurrentIncumbentFeas() == nullptr) &&
            ep.isFeasible(completeComputeType))
        {
            success = true;
            break;
        }

        // 3- Feasible case: check that ep is non-dominated.
        if (ep.isFeasible(completeComputeType))
        {
            bool isDominated = false;
            for (const auto& pElt: _paretoElements)
            {
                const auto compFlag = pElt->compMO(ep, completeComputeType);
                isDominated = (compFlag == NOMAD::CompareType::DOMINATING ||
                               compFlag == NOMAD::CompareType::EQUAL);
                if (isDominated)
                {
                    break;
                }
            }
           if (isDominated)
           {
               continue;
           }
           success = true;
           break;
        }

        // 3- Infeasible case: check that ep is non-dominated and below the current threshold.
        if (ep.getH(completeComputeType) > hMax)
        {
            continue;
        }

        bool isDominated = false;
        for (const auto& pElt: _paretoElements)
        {
            const auto compFlag = pElt->compMO(ep, completeComputeType, true /*onlyfvalues*/);
            isDominated = (compFlag == NOMAD::CompareType::DOMINATING) ||
                          (compFlag == NOMAD::CompareType::EQUAL);
            if (isDominated)
            {
                break;
            }
        }
        if (!isDominated)
        {
            success = true;
            break;
        }
    }

    // Determine success type according to DMultiMADS criterion
    auto megaIterSuccess = getSuccessType();
    for (auto& ep: _trialPoints)
    {
        
        // No need to continue; this step has already been marked as a success.
        if (megaIterSuccess == NOMAD::SuccessType::FULL_SUCCESS)
        {
            break;
        }
        
        // 0- The search has generated a point that is not eval ok
        if (!ep.isEvalOk(completeComputeType.evalType))
        {
            continue;
        }

        const auto epPtr = std::make_shared<EvalPoint>(ep);
        megaIterSuccess = ep.isFeasible(completeComputeType) ? _ref_dmads_barrier->getSuccessTypeOfPoints(epPtr, nullptr)
                                                             : _ref_dmads_barrier->getSuccessTypeOfPoints(nullptr, epPtr);
    }

    // Transfer all points from _trialPoints list in evalPointList
    std::vector<NOMAD::EvalPoint> evalPointList;
    for (const auto& ep: _trialPoints)
    {
        evalPointList.push_back(ep);
    }
    _trialPoints.clear();

    // Update the DMultiMads barrier with Quad DMS search points.
    _ref_dmads_barrier->updateWithPoints(evalPointList);

    // Set the success type of the current iteration
    setSuccessType(std::max(megaIterSuccess, getSuccessType()));

    return success;
}

void NOMAD::DMultiMadsQuadDMSSearchMethod::generateTrialPointOnSingleObjCombination()
{
    OUTPUT_DEBUG_START
    std::string s = "DMulti-MADS Quad DMS search: single-objective quadratic optimization launched ";
    s += "for objectives combination ( ";
    for (const auto ind: _activeObjsIndex)
    {
        s += "f" + std::to_string(ind) + " ";
    }
    s += ")";
    AddOutputDebug(s);
    OUTPUT_DEBUG_END
    
    // Save the current frame incumbent.
    auto dMadsBarrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(_ref_dmads_barrier);
    const NOMAD::EvalPointPtr currentFrameCenter = dMadsBarrier->getCurrentIncumbentFeas() != nullptr
                                                   ? dMadsBarrier->getCurrentIncumbentFeas()
                                                   : dMadsBarrier->getCurrentIncumbentInf();
    
    const auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    
    // Define the objective function to use for quad model optimization.
    prepareSingleObjectiveRun();
    
    NOMAD::QuadModelSinglePass singlePass(this, currentFrameCenter, currentFrameCenter->getMesh(),
                                          {} /* no scaling direction */,
                                          true /* true: prior combine objs */);
    
    // Generate the trial points
    singlePass.generateTrialPoints();
    
    // Pass the generated trial pts to this
    const auto& trialPtsSinglePass = singlePass.getTrialPoints();
    for (auto evalPoint : trialPtsSinglePass)
    {
        evalPoint.setPointFrom(currentFrameCenter, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        evalPoint.addGenStep(getStepType(), false /*do not inherit -> just qms step*/);
        
        evalPoint.setMesh(megaIter->getMesh());
        
        if (snapPointToBoundsAndProjectOnMesh(evalPoint, _lb, _ub))
        {
            evalPoint.addGenStep(getStepType());
            bool inserted = insertTrialPoint(evalPoint);

            OUTPUT_INFO_START
            std::string s = "xt:";
            s += (inserted) ? " " : " not inserted: ";
            s += evalPoint.display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
    }
}

void NOMAD::DMultiMadsQuadDMSSearchMethod::generateTrialPointsFinal()
{
    throw NOMAD::Exception(__FILE__, __LINE__, "Not yet implemented");
    
}

bool NOMAD::DMultiMadsQuadDMSSearchMethod::changeLevelAndUpdateIndex()
{
    // Should we increment _l
    // _index[0] cannot exceed _m - _l
    if (_indices[0] + _l == _m)
    {
        // Combinations of max _m objectives
        if (_l + 1 < _m)
        {
            // Start a new level of combinations
            // Example _l = 1 (combination of two), _m = 4
            // Row 1 [1 2] [1 3] [1 4]
            // Row 2 [2 3] [2 4]
            // Row 3 [3 4]
            // Row 3, i0 +_l = 3 + 1 = 4
            // Increment level: _l = 2 (combination of three objectives).
            // Row 1, first combination [ 1 2 3]
            _l++;
            _indices.clear();
            for (size_t l = 1 ; l <= _l+1 ; l++)
            {
                _indices.push_back(l);
            }
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        // Do not increment l, switch row
        // Example _l = 1 (combination of two), _m = 4
        // [1 2] [1 3] [1 4] --> row one is done
        // We need to do [2 3] [2 4]
        size_t I0 = _indices[0];
        _indices.clear();
        for (size_t l = 1 ; l <= _l+1 ; l++)
        {
            _indices.push_back(l + I0);
        }
        return true;
    }
}

bool NOMAD::DMultiMadsQuadDMSSearchMethod::selectObjCombination()
{
    // Perform initialization.
    if ( _indices.empty() )
    {
        _indices.push_back(1);
        _activeObjsIndex.push_back(1);
        return true;
    }
    
    // Next combination
    bool selectWithSuccess = true;
    if (_indices[_l] < _m)
    {
        // Example: _l = 1 (combinations of two), _m=4
        // Possible combinations of row 1: [1 2] [1 3] [1 4]
        // If we have done [1 3] we can increment the last index: [1 4]
        _indices[_l]++;
    }
    else
    {
        // An index has reached m. Switch row or done.
        selectWithSuccess = changeLevelAndUpdateIndex();
    }
    
    _activeObjsIndex.clear();
    if (!selectWithSuccess)
    {
        return false;
    }
    
    // Construct obj index
    for (const auto i: _indices)
    {
        _activeObjsIndex.push_back(i);
    }
    return true;
}


void NOMAD::DMultiMadsQuadDMSSearchMethod::testObjCombinations()
{
    std::cout << "Objective combinations: ";
    while (selectObjCombination())
    {
        std::cout <<" [ ";
        for(const auto i: _activeObjsIndex)
        {
            std::cout << i << " " ;
        }
        std::cout << " ] ";
    }
}
