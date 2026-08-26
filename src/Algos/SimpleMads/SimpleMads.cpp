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

#include "../../Algos/SimpleMads/SimpleMads.hpp"
#include "../../Algos/Mads/SinglePollMethod.hpp"
#include "../../Algos/SubproblemManager.hpp"


void NOMAD::SimpleMads::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_MADS);

    // We can accept Mads with more than one objective when doing a PhaseOneSearch of DMultiMads optimization.
    if (!_runParams->getAttributeValue<bool>("DMULTIMADS_OPTIMIZATION") && NOMAD::Algorithm::getNbObj() > 1)
    {
        throw NOMAD::InvalidParameter(__FILE__,__LINE__,"Mads solves single objective problems. To handle several objectives please use DMultiMads: DMULTIMADS_OPTIMIZATION yes");
    }

    if (nullptr == _poll.getBarrier())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"A simple mads must have a poll barrier.");
    }

    _vnsEnabled = _runParams->getAttributeValue<bool>("SIMPLE_MADS_WITH_VNS");
    if (_vnsEnabled)
    {

        if (_model !=nullptr)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Simple Mads with VNS: for now, cannot run VNS with model.");
        }

        if (_runParams->getAttributeValue<bool>("ANISOTROPIC_MESH"))
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Simple Mads with VNS: disable anisotropic mesh.");
        }

        _rng = std::make_shared<NOMAD::SimpleRNG>();

        _lb = _pbParams->getAttributeValue<ArrayOfDouble>("LOWER_BOUND");
        _ub = _pbParams->getAttributeValue<ArrayOfDouble>("UPPER_BOUND");

        // Strategy for _neighParameter update depends on the bounds
        // The VNS is conducted within a box or not
        if(_lb.isComplete() && _ub.isComplete())
        {
            // We have a box
            _neighParameter = 0.1;
            _boxedVNS = true;
        }
        else
        {
            _neighParameter = 1.0;
            _boxedVNS = false;
        }

    }

}

bool NOMAD::SimpleMads::runImp()
{
    size_t k = 0;   // Iteration number (incremented at start)

    bool runOk = true;

    // X0 eval in poll should be counted because cache is not used
    _nbEval++;

    // Termination: 1- Reach max model eval, 2- fail to create trial points on mesh
    while ( _nbEval < _maxEval && runOk)
    {

        _poll.start();
        runOk = _poll.run();
        _poll.end();
        _nbEval += _poll.getSingleNbEval();

        // Unlike in regular Mads, the VNS Search is done only when when poll is unsuccessful
        // If run not ok we stop looping
        if (runOk && _poll.getSuccessType() == NOMAD::SuccessType::UNSUCCESSFUL &&  _vnsEnabled)
        {
            bool vnsRunOk = runVNSSearch();
            if (vnsRunOk)
            {
                updatePollBarrierFromVNS();
            }
        }

        k++;
    }
    if (k == 1 && !runOk)
    {
        return false;
    }

    return true;
}

const NOMAD::SimpleEvalPoint & NOMAD::SimpleMads::getBestSimpleSolution(bool bestFeas) const
{

    const auto & barrier = _poll.getBarrier();
    bool runInPhaseOneSearch = _poll.getPhaseOneSearch();

    if (bestFeas && !runInPhaseOneSearch)
    {
        // If poll is in PhaseOneSearch the incumbentFeas is
        // PhaseOne feasible but not Standard feasible
        return barrier->getCurrentIncumbentFeas();
    }
    else
    {
        // If poll is in PhaseOneSearch the incumbentFeas is
        // Standard infeasible
        if (!runInPhaseOneSearch)
        {
            return barrier->getCurrentIncumbentInf();
        }
        else
        {
            return barrier->getCurrentIncumbentFeas();
        }
    }

}


NOMAD::EvalPoint NOMAD::SimpleMads::getBestSolution(bool bestFeas) const
{
    NOMAD::EvalPoint bestSol;

    const auto & barrier = _poll.getBarrier();
    bool runInPhaseOneSearch = _poll.getPhaseOneSearch();
    if (nullptr != barrier)
    {
        if (bestFeas)
        {
            // If poll is in PhaseOneSearch the incumbentFeas is
            // PhaseOne feasible but not Standard feasible
            if (!runInPhaseOneSearch)
            {
                bestSol = static_cast<NOMAD::EvalPoint>(barrier->getCurrentIncumbentFeas());
            }
        }
        else
        {
            // If poll is in PhaseOneSearch the incumbentFeas is
            // Standard infeasible
            if (!runInPhaseOneSearch)
            {
                bestSol = static_cast<NOMAD::EvalPoint>(barrier->getCurrentIncumbentInf());
            }
            else
            {
                bestSol = static_cast<NOMAD::EvalPoint>(barrier->getCurrentIncumbentFeas());
            }
        }

        if (bestSol.isComplete())
        {
            auto fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this);
            bestSol = bestSol.makeFullSpacePointFromFixed(fixedVariable);
        }
    }
    return bestSol;
}


void NOMAD::SimpleMads::endImp()
{

    if ( _endDisplay )
    {
        endDisplay();
    }
}


void NOMAD::SimpleMads::endDisplay() const
{
    // Display the best feasible solutions.
    std::string sFeas;

    // Output level is info if this algorithm is a sub part of another algorithm.
    NOMAD::OutputLevel outputLevel = NOMAD::OutputLevel::LEVEL_INFO;

    if ( ! NOMAD::OutputQueue::GoodLevel(outputLevel) )
    {
        return;
    }

    NOMAD::FHComputeTypeS computeType; /*default initializer is used*/

    auto solFormat = NOMAD::OutputQueue::getInstance()->getSolFormat();

    NOMAD::OutputInfo displaySolFeas(getName(), sFeas, outputLevel);

    sFeas = "Best feasible solution";

    auto bestFeas = getBestSolution(true);

    if (!bestFeas.isComplete())
    {
        sFeas += ":     Undefined.";
        displaySolFeas.addMsg(sFeas);
    }
    else
    {
        sFeas += ":     ";
        displaySolFeas.addMsg(sFeas + bestFeas.display(computeType,
                                                       solFormat,
                                                       NOMAD::DISPLAY_PRECISION_FULL,
                                                       false));
    }

    NOMAD::OutputQueue::Add(std::move(displaySolFeas));

    // Display the best infeasible solutions.
    std::string sInf;
    NOMAD::OutputInfo displaySolInf(getName(), sInf, outputLevel);
    sInf = "Best infeasible solution";

    auto bestInf = getBestSolution(false);

    if ( !bestInf.isComplete())
    {
        sInf += ":   Undefined.";
        displaySolInf.addMsg(sInf);
    }
    else
    {
        sInf += ":   ";
        displaySolInf.addMsg(sInf + bestInf.display(computeType,
                                                    solFormat,
                                                    NOMAD::DISPLAY_PRECISION_FULL,
                                                    false));
    }
    NOMAD::OutputQueue::Add(std::move(displaySolInf));

    std::string sNbEval = "Function evaluations: " + NOMAD::itos(_poll.getNbEval());
    NOMAD::OutputQueue::Add(sNbEval);
}


bool NOMAD::SimpleMads::runVNSSearch()
{
    if (!_vnsEnabled)
        return false;

    bool success = false;

    // Max eval for this single VNS search
    const size_t maxSingleVNSEval = _maxEval - _nbEval;
    if (maxSingleVNSEval > 0 && static_cast<double>(_nbEvalByVNS)/static_cast<double>(_nbEval) < _vnsTrigger)
    {

        // Manage the display with a block
        OUTPUT_INFO_START
        AddOutputInfo("Start Simple Mads VNS Search",true,false);
        OUTPUT_INFO_END

        auto frameSize = _poll.getMesh()->getDeltaFrameSize();
        auto initialFrameSize = _poll.getMesh()->getInitialFrameSize();

        // Continue only if current frame size is smaller than initial frame size
        if ( initialFrameSize < frameSize )
        {
            return false;
        }

        // MegaIteration's barrier member is already in sub dimension.
        auto bestXFeas = _poll.getBarrier()->getCurrentIncumbentFeas();
        auto bestXInf  = _poll.getBarrier()->getCurrentIncumbentInf();

        NOMAD::Point frameCenter;
        if (bestXFeas.isDefined())
        {
            frameCenter = bestXFeas;
        }
        else if (bestXInf.isDefined())
        {
            frameCenter = bestXInf;
        }
        else
            return false;

        size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
        NOMAD::Direction dirUnit(n, 0.0);
        NOMAD::Direction scaledDir(n, 0.0);
        NOMAD::Direction::computeDirOnUnitSphere(dirUnit, _rng);
        if (!dirUnit.isDefined())
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"VNS: single scaled direction not defined");
        }
        // Scaling direction with poll mesh
        auto pollMesh = _poll.getMesh();

        // Compute infinite norm for direction pointed by itDir.
        NOMAD::Double infiniteNorm = dirUnit.infiniteNorm();
        if (0 == infiniteNorm)
        {
            std::string err("Cannot handle an infinite norm of zero");
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        if (_boxedVNS)
        {

            // NOTE: Strategy suggested by S. Alarie

            scaledDir = dirUnit.forceExtendToBounds(frameCenter, _lb, _ub);
        }
        else
        {
            for (size_t i = 0; i < n; ++i)
            {
                // Scaling and projection on the mesh
                scaledDir[i] = pollMesh->scaleAndProjectOnMesh(i, dirUnit[i] / infiniteNorm);
            }
            scaledDir = dirUnit;
        }

        // Multiply shake direction by VNS neighborhood parameter
        scaledDir *= _neighParameter;

        OUTPUT_INFO_START
        AddOutputInfo("Shaking direction scaled: " + scaledDir.display());
        OUTPUT_INFO_END

        // Create shake point. Project on mesh + snap to bounds.
        NOMAD::Point shakePoint(frameCenter+scaledDir);
        shakePoint = pollMesh->projectOnMesh(shakePoint, frameCenter);
        shakePoint.snapToBounds(_lb, _ub);

        OUTPUT_INFO_START
        AddOutputInfo("Shake point after scaling, projecting on mesh and snapping to bounds: " + shakePoint.display());
        OUTPUT_INFO_END

        if ( shakePoint == frameCenter)
        {
            // Update the neighParameter
            if (_boxedVNS)
            {
                // For a boxed vns: loop 0.1, ..., 1.0, 0.1, 0.2, ....
                _neighParameter = (_neighParameter >= 1.0) ? 0.1 : (_neighParameter + 0.1);
            }
            else
            {
                // Simple update.
                // IMPORTANT: the increase maybe insufficient with respect to the mesh decrease -> shakePoint could be very close to frameCenter.
                _neighParameter ++;
            }
            OUTPUT_INFO_START
            AddOutputInfo("Shaking failed. ShakePoint=FrameCenter.");
            AddOutputInfo("End Simple Mads VNS Search",false,true);
            OUTPUT_INFO_END
            return false;
        }
        _vnsPoll = std::make_unique<NOMAD::SimplePoll>(this, _eval_x, shakePoint);

        bool vnsPollOk = true;
        bool frameSizeTooSmall = false;
        size_t k = 0;
        // Termination: 1- Reach max model eval, 2- fail to create trial points on mesh, 3- mesh size too small
        while ( _vnsPoll->getNbEval() < maxSingleVNSEval && vnsPollOk && !frameSizeTooSmall)
        {
            _vnsPoll->start();
            vnsPollOk = _vnsPoll->run();
            _vnsPoll->end();

            frameSizeTooSmall = (_vnsPoll->getMesh()->getDeltaFrameSize() < frameSize  );
            k++;
        }
        // Update the eval counters
        _nbEvalByVNS += _vnsPoll->getNbEval();
        _nbEval += _vnsPoll->getNbEval();

        if (k>0)
        {
            success = true;
        }

        // Update the neighParameter
        if (_boxedVNS)
        {
            // For a boxed vns: loop 0.1, ..., 1.0, 0.1, 0.2, ....
            _neighParameter = (_neighParameter >= 1.0) ? 0.1 : (_neighParameter + 0.1);
        }
        else
        {
            // Simple update.
            // IMPORTANT: the increase maybe insufficient with respect to the mesh decrease -> shakePoint could be very close to frameCenter.
            _neighParameter ++;
        }

        // Manage the display with a block
        OUTPUT_INFO_START
        AddOutputInfo("End Simple Mads VNS Search",false,true);
        OUTPUT_INFO_END

    }

    return success;
}


void NOMAD::SimpleMads::updatePollBarrierFromVNS()
{
    const auto & vnsBarrier = _vnsPoll->getBarrier();

    // Get vns best
    std::vector<SimpleEvalPoint> vnsBest;

    auto bestOne = vnsBarrier->getCurrentIncumbentFeas();
    if (bestOne.isDefined())
    {
        vnsBest.push_back(bestOne);
    }
    auto bestTwo = vnsBarrier->getCurrentIncumbentInf();
    if (bestTwo.isDefined())
    {
        vnsBest.push_back(bestTwo);
    }

    // Update poll barrier with best points
    auto & pollBarrier = _poll.getBarrier();
    bool updated = pollBarrier->updateWithPoints(vnsBest);
}
