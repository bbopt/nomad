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
#include <algorithm>

#include "../../Cache/CacheBase.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/DMultiMads/DMultiMadsIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMegaIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsExpansionIntLineSearchMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"

void NOMAD::DMultiMadsExpansionIntLineSearchMethod::init()
{
    if (nullptr == getParentOfType<DMultiMadsIteration*>())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMadsExpansionIntLineSearchMethod only works for DMultiMads");
    }

    setStepType(NOMAD::StepType::SEARCH_METHOD_DMULTIMADS_EXPANSIONINT_LINESEARCH);

    _bbInputTypes = _pbParams->getAttributeValue<NOMAD::BBInputTypeList>("BB_INPUT_TYPE");
    const bool hasIntegerVariables = !_bbInputTypes.empty() && std::any_of(_bbInputTypes.cbegin(), _bbInputTypes.cend(),
                                                                           [](const NOMAD::BBInputType bbi)
                                                                           {
                                                                               return bbi == NOMAD::BBInputType::INTEGER;
                                                                           });

    const bool runExpansionLinesearch = _runParams->getAttributeValue<bool>("DMULTIMADS_EXPANSIONINT_LINESEARCH");
    const bool isEnabled = runExpansionLinesearch && hasIntegerVariables;
    setEnabled(isEnabled);

    // Save lower and upper bounds
    _lb = getPbParams()->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    _ub = getPbParams()->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");

    // Deactivate the projection on the mesh for this search, as we only change integer variables.
    _projectOnMesh = false;
}

void NOMAD::DMultiMadsExpansionIntLineSearchMethod::preRunValidations()
{
    std::shared_ptr<NOMAD::BarrierBase> barrier = nullptr;

    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    const auto refComputeType = evc->getComputeType();

    if (NOMAD::ComputeType::STANDARD != refComputeType)
    {
        NOMAD::Exception(__FILE__,__LINE__,"Cannot do expansion integer linesearch for DMultiMads on ComputeType other than STANDARD.");
    }

    // Get barrier from upper MegaIteration, if available.
    const auto megaIter = getParentOfType<NOMAD::DMultiMadsMegaIteration*>(false);
    if (megaIter != nullptr)
    {
        barrier = megaIter->getBarrier();
    }

    if (nullptr == barrier || nullptr == std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(barrier))
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"For DMultiMads using expansion integer linesearch, we need a DMultiMadsBarrier.");
    }

    if (evc->getCurrentEvalType() == NOMAD::EvalType::MODEL)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot do expansion integer linesearch for DMultiMads on EvalType::MODEL.");
    }

    // Keep a reference to the DMultiMads barrier.
    _ref_dmads_barrier = barrier;
}

NOMAD::Direction NOMAD::DMultiMadsExpansionIntLineSearchMethod::computePrimitiveDirection(const NOMAD::Point& frameCenter,
                                                                                          const NOMAD::Point& pointFrom,
                                                                                          int& initStepSize) const
{
    // Compute the generative direction for integer variables
    bool hasRealDirCoordinates = false;
    auto dir = NOMAD::Point::vectorize(pointFrom, frameCenter);
    for (size_t i = 0; i < dir.size(); ++i)
    {
        if (_bbInputTypes[i] != NOMAD::BBInputType::INTEGER)
        {
            if (dir[i] != 0)
            {
                hasRealDirCoordinates = true; // The point was generated from non-integer coordinates
            }
            dir[i] = 0; // We will not move along these coordinates
        }
        else
        {
            dir[i] = dir[i].round(); // Be sure to round to the nearest integer
        }
    }
    
    auto gcd = [](int a, int b)
    {
        a = std::abs(a);
        b = std::abs(b);
        if (a == 0 || b == 0)
        {
            return std::max(a, b);
        }

        while (b != 0)
        {
            int t = b;
            b = a % b;
            a = t;
        }
        return a;
    };

    // Extract the gcd of all dir coordinates
    int gDivisor = dir[0].round();
    for (size_t i = 1; i < dir.size(); ++i)
    {
        gDivisor = gcd(gDivisor, dir[i].round());
    }

    // The frame center has been generated from a direction composed only of non-integer directions
    if (gDivisor == 0)
    {
       dir *= gDivisor;
       return dir;
    }

    // Save the initStepSize: it is equal to gcd if gcd is a multiple of 2 and the frame center has been generated
    // from integer directions; otherwise -1.
    // NB: an explanation for this black-magic code can be found here:
    // https://stackoverflow.com/questions/57025836/how-to-check-if-a-given-number-is-a-power-of-two
    const bool isPowerOf2 = (gDivisor > 0) && ((gDivisor & (gDivisor-1)) == 0);
    initStepSize = !hasRealDirCoordinates && isPowerOf2 ? gDivisor : -1;

    // Compute the primitive direction
    for (size_t i = 0; i < dir.size(); ++i)
    {
        dir[i] = dir[i].round() / gDivisor;
    }

    return dir;
}

int NOMAD::DMultiMadsExpansionIntLineSearchMethod::computeMaxStepSize(const NOMAD::Point &frameCenter,
                                                                      const NOMAD::Direction& dir,
                                                                      const NOMAD::ArrayOfDouble &lb,
                                                                      const NOMAD::ArrayOfDouble &ub) const
{
    // Compute maximal step size to reach bounds
    int stepSize = std::numeric_limits<int>::max();
    for (size_t i = 0; i < dir.size(); ++i)
    {
        if (dir[i] == 0)
        {
            continue;
        }

        // We will never reach bounds according to these coordinates
        if (dir[i] > 0 && (!ub[i].isDefined() || ub[i] == NOMAD::INF))
        {
            continue;
        }
        if (dir[i] < 0 && (!lb[i].isDefined() || lb[i] == NOMAD::M_INF))
        {
            continue;
        }

        const NOMAD::Double xi = frameCenter[i];
        if (dir[i] > 0)
        {
            const NOMAD::Double ubmx = ub[i] - xi;
            stepSize = std::min(stepSize, ubmx.round() / std::abs(dir[i].round()));
            continue;
        }

        if (dir[i] < 0)
        {
            const NOMAD::Double lbmx =  xi - lb[i];
            stepSize = std::min(stepSize, lbmx.round() / std::abs(dir[i].round()));
        }
    }

    // No bound is defined. We affect a stepsize value equal to 10
    if (stepSize == std::numeric_limits<int>::max())
    {
        stepSize = 10;
    }
    return stepSize;
}

bool NOMAD::DMultiMadsExpansionIntLineSearchMethod::isInBarrier(const NOMAD::Point& x) const
{
    const auto dMadsBarrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(_ref_dmads_barrier);
    for (size_t i = 0; i < dMadsBarrier->nbXFeas(); ++i)
    {
        if (*dMadsBarrier->getXFeas(i).getX() == x)
        {
            return true;
        }
    }

    for (size_t i = 0; i < dMadsBarrier->nbXFilterInf(); ++i)
    {
        if (*dMadsBarrier->getXFilterInf(i).getX() == x)
        {
            return true;
        }
    }
    return false;
}

void NOMAD::DMultiMadsExpansionIntLineSearchMethod::generateTrialPointsFinal()
{
    preRunValidations();

    const auto dMadsBarrier = std::dynamic_pointer_cast<NOMAD::DMultiMadsBarrier>(_ref_dmads_barrier);

    const auto currentFrameCenter = dMadsBarrier->getCurrentIncumbentFeas() != nullptr ?
            dMadsBarrier->getCurrentIncumbentFeas() : dMadsBarrier->getCurrentIncumbentInf();
    if (currentFrameCenter == nullptr)
    {
        return;
    }

    const auto pointFrom = currentFrameCenter->getPointFrom(NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
    if (nullptr == pointFrom || *pointFrom == *currentFrameCenter)
    {
        OUTPUT_INFO_START
        AddOutputInfo("No available direction: stop");
        OUTPUT_INFO_END
        return;
    }

    // Compute the primitive direction
    int initStepSize = -1;
    auto dir = computePrimitiveDirection(*currentFrameCenter, *pointFrom, initStepSize);
    if (dir.norm() == 0)
    {
        OUTPUT_INFO_START
        AddOutputInfo("No available primitive direction: stop");
        OUTPUT_INFO_END
        return;
    }

    OUTPUT_INFO_START
    const auto initDir = NOMAD::Point::vectorize(*pointFrom, *currentFrameCenter);
    AddOutputInfo("Frame center: " + currentFrameCenter->display());
    AddOutputInfo("Initial direction: " + initDir.display());
    AddOutputInfo("Primitive direction: " + dir.display());
    OUTPUT_INFO_END

    // Compute the maximum step size parameter to reach problem integer bounds.
    const int maxStepSize = computeMaxStepSize(*currentFrameCenter, dir, _lb, _ub);

    // We have already reach one of the problem bounds
    if (maxStepSize == 0)
    {
        OUTPUT_INFO_START
        AddOutputInfo("Max step size parameter value is 0: stop");
        OUTPUT_INFO_END
        return;
    }

    // Obtain cache.
    const auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    const bool useCache = evc->getUseCache();

    // We try to guess if the frame center has already been generated from a primitive direction multiplied by
    // some scalar, i.e. initStepSize, which must be superior to 1 if it is the case; and is equal to -1 otherwise.
    // If the frame center has been generated from a primitive direction, we start a search with a step size parameter
    // equal to 2 * initStepSize; but we must stay in bounds.
    int stepSize = std::min(std::max(1, 2 * initStepSize), maxStepSize);
    auto scaledDir = dir;
    scaledDir *= stepSize;

    OUTPUT_INFO_START
    AddOutputInfo("Step size parameter: " + std::to_string(stepSize));
    AddOutputInfo("Scaled direction: " + scaledDir.display());
    OUTPUT_INFO_END

    // First candidate tentative
    auto ev1 = NOMAD::EvalPoint(*currentFrameCenter->getX() + scaledDir);
    ev1.setPointFrom(std::make_shared<NOMAD::EvalPoint>(*currentFrameCenter),
                           NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
    ev1.addGenStep(getStepType());
    snapPointToBoundsAndProjectOnMesh(ev1, _lb, _ub);

    const bool inBarrier1 = isInBarrier(*ev1.getX());

    // In this case, the candidate has never been evaluated.
    if (!useCache && !inBarrier1)
    {
        insertTrialPoint(ev1);
        return;
    }

    if (useCache)
    {
        NOMAD::CacheInterface cacheInterface(this);
        NOMAD::EvalPoint foundEvalPoint;

        // The point has not been evaluated. Generate it.
        if (!cacheInterface.find(*ev1.getX(), foundEvalPoint))
        {
            insertTrialPoint(ev1);
            return;
        }
    }

    // At this point, the candidate is already in the cache.
    OUTPUT_INFO_START
    std::string s = "Search candidate ";
    s += ev1.getX()->display();
    s += " already in cache";
    AddOutputInfo(s);
    OUTPUT_INFO_END

    // We cannot further decrease the step size parameter.
    if (stepSize == 1)
    {
        OUTPUT_INFO_START
        AddOutputInfo("Step size parameter has reached 0: stop");
        OUTPUT_INFO_END
        return;
    }

    // Divide the stepsize by 2 and retry.
    stepSize = std::max(1, stepSize / 2);
    scaledDir = dir;
    scaledDir *= stepSize;

    OUTPUT_INFO_START
    AddOutputInfo("New tentative");
    AddOutputInfo("Step size parameter: " + std::to_string(stepSize));
    AddOutputInfo("Scaled direction: " + scaledDir.display());
    OUTPUT_INFO_END

    auto ev2 = NOMAD::EvalPoint(*currentFrameCenter->getX() + scaledDir);

    ev2.setPointFrom(std::make_shared<NOMAD::EvalPoint>(*currentFrameCenter),
                           NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
    ev2.addGenStep(getStepType());
    snapPointToBoundsAndProjectOnMesh(ev2, _lb, _ub);

    const bool inBarrier2 = isInBarrier(*ev2.getX());

    if (!useCache && !inBarrier2)
    {
        insertTrialPoint(ev2);
        return;
    }

    if (useCache)
    {
        NOMAD::CacheInterface cacheInterface(this);
        NOMAD::EvalPoint foundEvalPoint;

        // The point has not been evaluated. Generate it.
        if (!cacheInterface.find(*ev2.getX(), foundEvalPoint))
        {
            insertTrialPoint(ev2);
            return;
        }
    }

    OUTPUT_INFO_START
    std::string s = "Search candidate n2 ";
    s += ev1.getX()->display();
    s += " already in cache: stop.";
    AddOutputInfo(s);
    OUTPUT_INFO_END
}
