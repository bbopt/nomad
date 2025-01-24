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
#include "../../Cache/CacheBase.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/DMultiMads/DMultiMadsIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMegaIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsNMSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/NelderMead/NMAllReflective.hpp"
#include "../../Type/DMultiMadsSearchStrategyType.hpp"

void NOMAD::DMultiMadsNMSearchMethod::init()
{
    if (nullptr == getParentOfType<DMultiMadsIteration*>())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMadsNMSearch only works for DMultiMads");
    }

    const bool isEnabled = getRunParams()->getAttributeValue<bool>("NM_SEARCH");
    setEnabled(isEnabled);

    if (isEnabled)
    {
        // Set the lap counter
        const auto nmFactor = _runParams->getAttributeValue<size_t>("NM_SEARCH_MAX_TRIAL_PTS_NFACTOR");
        const auto dim = _pbParams->getAttributeValue<size_t>("DIMENSION");
        if (nmFactor < NOMAD::INF_SIZE_T)
        {
            NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( dim*nmFactor );
        }

        // NM is an algorithm with its own stop reasons.
        _nmStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::NMStopType>>();

        // Create the NM algorithm with its own stop reason
        _nm = std::make_unique<NOMAD::NM>(this,
                                          _nmStopReasons ,
                                          _runParams,
                                          _pbParams);

    }

    setStepType(NOMAD::StepType::SEARCH_METHOD_NM);

    const auto nmStrategy = getRunParams()->getAttributeValue<NOMAD::DMultiMadsNMSearchType>("DMULTIMADS_NM_STRATEGY");

    _use_dom_strategy = nmStrategy == NOMAD::DMultiMadsNMSearchType::DOM;
}

void NOMAD::DMultiMadsNMSearchMethod::preRunValidations()
{
    std::shared_ptr<NOMAD::BarrierBase> barrier = nullptr;

    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    _ref_compute_type = evc->getComputeType();

    if (NOMAD::ComputeType::STANDARD != _ref_compute_type)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do NM search for DMultiMads on ComputeType other than STANDARD.");
    }

    // Get barrier from upper MegaIteration, if available.
    const auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    if (megaIter != nullptr)
    {
        barrier = megaIter->getBarrier();
    }

    if (nullptr == barrier || nullptr == std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using NM search, we need a DMultiMadsBarrier.");
    }

    if (evc->getCurrentEvalType() == NOMAD::EvalType::MODEL)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do NM search for DMultiMads on EvalType::MODEL.");
    }
    if (!evc->getUseCache())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using NM search, we need a cache.");
    }

    // Keep a reference to the DMultiMads barrier.
    _ref_dmads_barrier = barrier;

    // And to the compute type.
    _ref_compute_type = evc->getComputeType();
}

bool NOMAD::DMultiMadsNMSearchMethod::runImp()
{
    // Basic checks before running Nelder-Mead search.
    preRunValidations();

    if (_use_dom_strategy)
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("Start DoM strategy");
        OUTPUT_DEBUG_END

        bool success = runDoMStrategy();
        return success;
    }

    OUTPUT_DEBUG_START
    AddOutputDebug("Start MultiMads strategy");
    OUTPUT_DEBUG_END

    bool success = runMultiMadsStrategy();
    return success;
}

bool NOMAD::DMultiMadsNMSearchMethod::runDoMStrategy()
{
    auto dMadsBarrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(_ref_dmads_barrier);

    // Get current frame center
    const auto currentFrameCenter = dMadsBarrier->getCurrentIncumbentFeas() != nullptr
            ? dMadsBarrier->getCurrentIncumbentFeas() : dMadsBarrier->getCurrentIncumbentInf();

    // Get all the points available in the Pareto front approximation.
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
    std::vector<NOMAD::EvalPointPtr> paretoElements = getCurrentParetoFrontApproximation(*dMadsBarrier);

    // Compute reference vector.
    // It is a vector of size:
    // (lk - 1) * nobj, if lk > 1;
    // nobj             otherwise
    // where lk is the number of elements in the Pareto front approximation and nobj the number of objectives.
    const size_t nbParetoElements = paretoElements.size();
    const size_t nbObj = dMadsBarrier->getNbObj();

    OUTPUT_DEBUG_START
    AddOutputDebug("Nb current Pareto elements: " + std::to_string(nbParetoElements));
    AddOutputDebug("Number of objectives: " + std::to_string(nbObj));
    AddOutputDebug("Current frame center is: " + currentFrameCenter->display());
    OUTPUT_DEBUG_END

    // NB: As this criterion is computationally intensive, we limit the use of NOMAD::Double
    // and work directly with std::double.
    std::vector<double> ref;
    ref.reserve(std::max((int)nbParetoElements - 1, 1) * nbObj);

    // It contains all elements of the Pareto front approximation minus the objective vector of the current frame
    // center if the Pareto front approximation is not composed of 1 element
    const FHComputeType& initFHComputeType = dMadsBarrier->getFHComputeType();
    if (nbParetoElements == 1)
    {
        const NOMAD::ArrayOfDouble& fvalues = currentFrameCenter->getFs(initFHComputeType);
        for (size_t i = 0; i < nbObj; ++i)
        {
            const NOMAD::Double fi = fvalues[i];
            ref.push_back(fi.todouble());
        }
    }
    else
    {
       for (size_t j = 0; j < nbParetoElements; ++j)
       {
           const auto& paretoElt = paretoElements[j];
           if (currentFrameCenter == paretoElt)
           {
               continue;
           }
           const auto& fvalues = paretoElt->getFs(initFHComputeType);
           for (size_t i = 0; i < nbObj; ++i)
           {
               ref.push_back(fvalues[i].todouble());
           }
       }
    }

    OUTPUT_DEBUG_START
    // NB: The reference vector is costly to print
    AddOutputDebug("Reference vector of dimensions " + std::to_string(ref.size() / nbObj) + " x " + std::to_string(nbObj));
    for (size_t j = 0; j < ref.size() / nbObj; ++j)
    {
        std::string s;
        for (size_t i = 0; i < nbObj; ++i)
        {
            s += " " + std::to_string(ref[j * nbObj + i]);
        }
        AddOutputDebug(s);
    }
    OUTPUT_DEBUG_END


    // We also allocate a vector that will contain the values of objective function.
    // These values can be loaded in processor cache, allowing for faster calculus
    std::vector<double> fvalues(nbObj);

    // Compute single-objective function for NM.
    NOMAD::singleOutputComputeFType singleObjCompute = [&ref, &fvalues, nbObj](const BBOutputTypeList& bbOutputTypeList,
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

        const size_t nbParetoElements = ref.size() / nbObj;

        const auto& objValues = bbOutput.getObjectives(bbOutputTypeList);
        for (int i = 0; i < nbObj; ++i)
        {
            fvalues[i] = objValues[i].todouble();
        }

        // The single-objective function is defined by:
        // psi(x) = - min_{y in R} sum_{i = 1}^nbObj max(0, yi - fi(x)) if fi(x) is not dominated by any element
        //                                                               of R
        //            min_{y in R} sum_{i = 1}^nbObj max(0, fi(x) - yi) otherwise.

        // First pass: detect if the point is dominated.
        // Compute doMValue = - min_{y in R} sum_{i = 1}^nbObj max(0, yi - fi(x))
        double doMValue = NOMAD::INF;
        for (size_t j = 0; j < nbParetoElements; ++j)
        {
            double curDoMValue = 0;
            for (size_t i = 0; i < fvalues.size(); i++)
            {
                curDoMValue += std::max(0.0, ref[j * nbObj + i] - fvalues[i]);
            }
            doMValue = std::min(doMValue, curDoMValue);

            // In this case, the point is dominated by the element of R or equal. No need to continue.
            if (doMValue == 0)
            {
                break;
            }
        }

        // The point is not dominated.
        // NB: The use of NOMAD::Double here enables to have a higher threshold for the positiveness
        // of this function than using std::double
        if (NOMAD::Double(doMValue) > 0)
        {
            return -doMValue;
        }

        // Second pass: return the minimum move to reach one of the elements of the Pareto front.
        // Compute doMValue = min_{y in R} sum_{i = 1}^nbObj max(0, fi(x) - yi).
        doMValue = NOMAD::INF;
        for (size_t j = 0; j < nbParetoElements; ++j)
        {
            double curDoMValue = 0;
            for (size_t i = 0; i < fvalues.size(); i++)
            {
                curDoMValue += std::max(0.0, fvalues[i] - ref[j * nbObj + i]);
            }
            doMValue = std::min(doMValue, curDoMValue);
            // The point is equal to one of the points in the Pareto front. No need to continue
            if (doMValue == 0)
            {
                break;
            }
        }

        return doMValue;
    };

    // Make the objective function available to the evaluator control.
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    evc->setComputeType(NOMAD::ComputeType::DMULTI_COMBINE_F, singleObjCompute);

    // Get all the points of the DMultiMads barrier
    const std::vector<NOMAD::EvalPoint> evalPointList = _ref_dmads_barrier->getAllPoints();

    // Create a new barrier for this NM launch: HMax is set to initial value of the DMultiMads barrier.
    const auto hMax = dMadsBarrier->getHMax();

    // Update with all the DMultiMads barrier points.
    auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    // Update with all the DMultiMads barrier points.
    megaIter->setBarrier(std::make_shared<NOMAD::ProgressiveBarrier>(hMax,
                                                                     NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                                                     evc->getCurrentEvalType(),
                                                                     evc->getFHComputeTypeS(),
                                                                     evalPointList));

    _tagBefore = NOMAD::EvalPoint::getCurrentTag();

    // Run Nelder-Mead search
    NOMAD::NMSearchMethod::runImp();
    const bool successPostProcessing = postRunUpdates();

    return successPostProcessing;
}

bool NOMAD::DMultiMadsNMSearchMethod::runMultiMadsStrategy()
{
    auto dMadsBarrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(_ref_dmads_barrier);

    // Compute the reference vector in the objective space.
    const auto ref = computeReferencePoint(*dMadsBarrier);
    const size_t nbObj = dMadsBarrier->getNbObj();
    std::vector<bool> launchSingleObjectiveRuns(nbObj, false);
    for (size_t obj = 0; obj < nbObj; ++obj)
    {
        launchSingleObjectiveRuns[obj] = !ref[obj].isDefined();
    }
    const auto launchSingleObjRuns = std::any_of(launchSingleObjectiveRuns.cbegin(),
                                                 launchSingleObjectiveRuns.cend(),
                                                 [](bool b) {return b;});

    bool success = false;

    // Two cases
    // 1- The current incumbent is an extreme solution. Launch a single-objective Nelder-Mead step for each objective
    //    for which the current incumbent is an extreme solution.
    if (launchSingleObjRuns)
    {
        for (size_t obj = 0; obj < nbObj; ++ obj)
        {
            if (launchSingleObjectiveRuns[obj])
            {
                OUTPUT_DEBUG_START
                    AddOutputDebug("DMulti-MADS NM search: single-objective optimization launch for objective f" + std::to_string(obj + 1));
                OUTPUT_DEBUG_END

                const NOMAD::ArrayOfDouble refSingleObj(nbObj);
                refSingleObj[obj] = 0; // Only the defined objective will be considered in the NM search
                prepareSingleObjectiveRun(refSingleObj);

                // Launch NM search
                NOMAD::NMSearchMethod::runImp();
                const bool successPostProcessing = postRunUpdates();
                success = successPostProcessing || success;
            }
        }
        return success;
    }

    // 2- The current incumbent is not an extreme solution. Launch a ``Multi-MADS'' search starting from the current
    //    incumbent using the reference vector.
    OUTPUT_DEBUG_START
        AddOutputDebug("DMulti-MADS NM search: reference point " + ref.display());
    OUTPUT_DEBUG_END

    prepareMultiMadsRun(ref);
    NOMAD::NMSearchMethod::runImp();
    success = postRunUpdates();

    return success;
}

void NOMAD::DMultiMadsNMSearchMethod::prepareSingleObjectiveRun(const NOMAD::ArrayOfDouble& ref)
{
    // Get all the points of the DMultiMads barrier
    const std::vector<NOMAD::EvalPoint> evalPointList = _ref_dmads_barrier->getAllPoints();

    // Define the single-objective function for NM run.
    NOMAD::singleOutputComputeFType singleObjCompute = [&ref](const BBOutputTypeList& bbOutputTypeList,
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

        const auto& bbo = bbOutput.getBBOAsArrayOfDouble();

        size_t j = 0;
        for (size_t i = 0; i < bbo.size(); i++)
        {
            if (bbOutputTypeList[i].isObjective() && ref[j].isDefined())
            {
                return bbo[i];
            }
            if (bbOutputTypeList[i].isObjective())
            {
                j++;
            }
        }

        // Normally, we should not go here.
        return NOMAD::INF;
    };

    // Make the objective function available to the evaluator control.
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    evc->setComputeType(NOMAD::ComputeType::DMULTI_COMBINE_F, singleObjCompute);

    // Create a new barrier for this NM launch: HMax is set to initial value of the DMultiMads barrier.
    const auto hMax = _ref_dmads_barrier->getHMax();

    // Set the mega iter barrier to a ProgressiveBarrier with all the DMultiMads barrier points.
    auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    megaIter->setBarrier(std::make_shared<NOMAD::ProgressiveBarrier>(hMax,
                                                                     NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                                                     evc->getCurrentEvalType(),
                                                                     evc->getFHComputeTypeS(),
                                                                     evalPointList));

    _tagBefore = NOMAD::EvalPoint::getCurrentTag();
}

void NOMAD::DMultiMadsNMSearchMethod::prepareMultiMadsRun(const NOMAD::ArrayOfDouble &ref)
{
    // Get all the points of the DMultiMads barrier
    const std::vector<NOMAD::EvalPoint> evalPointList = _ref_dmads_barrier->getAllPoints();
    const auto dMadsBarrier = std::dynamic_pointer_cast<DMultiMadsBarrier>(_ref_dmads_barrier);

    NOMAD::singleOutputComputeFType singleObjCompute = [&ref](const BBOutputTypeList& bbOutputTypeList,
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

        const auto& bbo = bbOutput.getBBOAsArrayOfDouble();
        const auto nbDominatedObj = [&]() -> size_t
        {
            size_t j = 0;
            size_t nbDomObj = 0;
            for (size_t i = 0; i < bbo.size(); ++i)
            {
                if (bbOutputTypeList[i].isObjective() && ref[j].isDefined())
                {
                    const NOMAD::Double& refValue = ref[j];
                    const NOMAD::Double& fValue = bbo[i];
                    // The objective component of f is dominated by the same objective component of ref
                    if (fValue >= refValue)
                    {
                        nbDomObj += 1;
                    }
                    j++;
                }
            }
            return nbDomObj;
        }(); // IIFE

        const size_t nbObj = getNbObj(bbOutputTypeList);
        NOMAD::Double minDistDominated = NOMAD::INF;
        NOMAD::Double minDistDominating = NOMAD::INF;
        NOMAD::Double distRefF = 0;
        bool isDominating = true;
        size_t j = 0;
        for (size_t i = 0; i < bbo.size(); i++)
        {
            if (bbOutputTypeList[i].isObjective() && ref[j].isDefined())
            {
                const NOMAD::Double& refValue = ref[j];
                const NOMAD::Double& fValue = bbo[i];
                const NOMAD::Double refmf = (refValue - fValue) * (refValue - fValue);
                if (fValue > refValue)
                {
                    isDominating = false;
                    minDistDominated = std::min(minDistDominated, refmf);
                }
                else
                {
                    minDistDominating = std::min(minDistDominating, refmf);
                }
                distRefF += refmf;
                j++;
            }
        }
        if (nbDominatedObj == nbObj)
        {
            // The objective vector is dominated by the reference point:
            // return || ref - f(x) ||^2
            return distRefF;
        }

        // The objective vector is not dominated by the reference point. Two cases can occur:
        // * f(x) and ref are indifferent: return the distance between the dominance zone and f(x),
        //   i.e., min_{i in I} || ref - f(x) ||^2, where I = {i : fi(x) > ref[i]}.
        // * f(x) is in the dominance zone: return the negative distance between the dominance zone and f(x),
        //   i.e., min || ref - f(x) ||^2.
        const NOMAD::Double objValue = isDominating ? -minDistDominating : minDistDominated;
        return objValue;
    };

    // Make the objective function available to the evaluator control.
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    evc->setComputeType(NOMAD::ComputeType::DMULTI_COMBINE_F, singleObjCompute);

    // Create a new barrier for this NM launch: HMax is set to initial value of the DMultiMads barrier.
    const auto hMax = _ref_dmads_barrier->getHMax();

    // Update with all the DMultiMads barrier points.
    auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    megaIter->setBarrier(std::make_shared<NOMAD::ProgressiveBarrier>(hMax,
                                                                     NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this),
                                                                     evc->getCurrentEvalType(),
                                                                     evc->getFHComputeTypeS(),
                                                                     evalPointList));

    _tagBefore = NOMAD::EvalPoint::getCurrentTag();
}

bool NOMAD::DMultiMadsNMSearchMethod::postRunUpdates()
{
    bool success = false;

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();

    // Get the progressive barrier temporarily stored into MegaIteration.
    auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    const std::shared_ptr<NOMAD::BarrierBase> barrier = megaIter->getBarrier();
    if (nullptr == barrier || nullptr == std::dynamic_pointer_cast<NOMAD::ProgressiveBarrier>(barrier))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads, after NM search, we should have a ProgressiveBarrier.");
    }

    // Reset all DMULTI_COMBINE_F_VALUES in cache
    if (NOMAD::CacheBase::getInstance() != nullptr)
    {
        NOMAD::CacheBase::getInstance()->processOnAllPoints(
                [](EvalPoint& evp)
                {
                    evp.resetDMultiCombineFValue();
                }
        );
        evc->getBestIncumbent(-1)->resetDMultiCombineFValue();
    }

    // Collect all eval points generated by the NM search.
    std::vector<NOMAD::EvalPoint> evalPointList;
    if (_nm->getSuccessType() > NOMAD::SuccessType::NO_TRIALS)
    {
        const auto evalType = evc->getCurrentEvalType();
        std::function<bool(const EvalPoint&)> func = [&, evalType](const EvalPoint& evalPoint) -> bool
        {
            if(evalPoint.getTag() > (int)_tagBefore && evalPoint.isEvalOk(evalType))
            {
                return true;
            }
            return false;
        };

        NOMAD::CacheInterface cacheInterface(this);
        cacheInterface.find(func, evalPointList);
    }

    // Set back the mega iteration barrier
    megaIter->setBarrier(_ref_dmads_barrier);

    // Set back the STANDARD compute type for DMultiMads
    evc->setComputeType(_ref_compute_type);

    const NOMAD::FHComputeType completeComputeType = {evc->getCurrentEvalType(),
                                                      evc->getFHComputeTypeS()};

    auto megaIterSuccess = NOMAD::SuccessType::UNSUCCESSFUL;
    // Set the mesh and determine successType
    for (NOMAD::EvalPoint & ep: evalPointList)
    {
        ep.setMesh(megaIter->getMesh());

        const auto epPtr = std::make_shared<EvalPoint>(ep);
        if (megaIterSuccess != NOMAD::SuccessType::FULL_SUCCESS)
        {
            _ref_dmads_barrier->checkForFHComputeType(completeComputeType);

            if (ep.isFeasible(completeComputeType))
            {
                megaIterSuccess = _ref_dmads_barrier->getSuccessTypeOfPoints(epPtr, nullptr);
            }
            else
            {
                megaIterSuccess = _ref_dmads_barrier->getSuccessTypeOfPoints(nullptr, epPtr);
            }
        }
    }

    // Update the DMultiMads barrier with NM search points.
    success = _ref_dmads_barrier->updateWithPoints(evalPointList);
    setSuccessType(megaIterSuccess);

    return success;
}

void NOMAD::DMultiMadsNMSearchMethod::generateTrialPointsFinal()
{
    throw NOMAD::Exception(__FILE__, __LINE__, "Not yet implemented");
}

NOMAD::ArrayOfDouble NOMAD::DMultiMadsNMSearchMethod::computeReferencePoint(const NOMAD::DMultiMadsBarrier& barrier) const
{
    // Get all the points available in the Pareto front approximation.
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
    std::vector<NOMAD::EvalPointPtr> evalPointList = getCurrentParetoFrontApproximation(barrier);

    const NOMAD::EvalPointPtr currentFrameCenter = barrier.getCurrentIncumbentFeas() != nullptr ? barrier.getCurrentIncumbentFeas()
                                                                                                : barrier.getCurrentIncumbentInf();

    const size_t nbObj = barrier.getNbObj();
    const FHComputeType& initFHComputeType = barrier.getFHComputeType();
    if (evalPointList.size() == 1)
    {
        return NOMAD::ArrayOfDouble(nbObj);
    }

    // When the ie coordinate of the reference point is not defined, it means the current frame incumbent
    // is an extreme solution of the Pareto front for the current ie objective.
    NOMAD::ArrayOfDouble ref(nbObj);
    for (size_t obj = 0; obj < nbObj; ++obj)
    {
        std::sort(evalPointList.begin(), evalPointList.end(),
                  [&initFHComputeType, obj](const NOMAD::EvalPointPtr& ev1, const EvalPointPtr& ev2)->bool
                  {
                      return ev1->getFs(initFHComputeType)[obj] < ev2->getFs(initFHComputeType)[obj];
                  });

        // Get extreme values value according to one objective
        NOMAD::Double fmin = evalPointList[0]->getFs(initFHComputeType)[obj];
        NOMAD::Double fmax = evalPointList[evalPointList.size()-1]->getFs(initFHComputeType)[obj];

        // The current frame incumbent is an extreme solution for obj.
        if (fmin == fmax)
        {
            continue;
        }

        // Find the current frame center
        size_t frameCenterId = 0;
        for (size_t i = 0; i < evalPointList.size(); ++i)
        {
            if (currentFrameCenter == evalPointList[i])
            {
                frameCenterId = i;
                break;
            }
        }

        // The current frame incumbent is an extreme solution for obj.
        if (frameCenterId == 0)
        {
            continue;
        }

        // Compute the reference objective vector. When frameCenterId is the last element of the list, the
        // reference vector takes the corresponding objective value.
        const size_t nextFId = std::min(frameCenterId + 1, evalPointList.size() - 1);
        const auto& nextF = evalPointList[nextFId]->getFs(initFHComputeType);
        ref[obj] = nextF[obj];
    }

    return ref;
}
