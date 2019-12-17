/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#include "../Algos/EvcInterface.hpp"

#include "../Algos/Mads/MadsIteration.hpp"
#include "../Algos/Mads/MegaSearchPoll.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
std::shared_ptr<NOMAD::EvaluatorControl> NOMAD::EvcInterface::_evaluatorControl = nullptr;


void NOMAD::EvcInterface::init()
{
    verifyStepNotNull();
    verifyEvaluatorControlNotNull();
}


void NOMAD::EvcInterface::verifyStepNotNull()
{
    if (nullptr == _step)
    {
        std::string err = "Step for EvcInterface should not be NULL";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::EvcInterface::verifyEvaluatorControlNotNull()
{
    if (nullptr == _evaluatorControl)
    {
        std::string err = "EvaluatorControl for EvcInterface should not be NULL";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::EvcInterface::setEvaluatorControl(const std::shared_ptr<NOMAD::EvaluatorControl>& evaluatorControl)
{
    _evaluatorControl = evaluatorControl;
    verifyEvaluatorControlNotNull();
}


// For each point, look if it is in the cache.
// If it is, remove it from the EvalPointSet.
// If not, add it to EvaluatorControl's Queue.
void NOMAD::EvcInterface::keepPointsThatNeedEval(const NOMAD::EvalPointSet &trialPoints, bool useMesh)
{
    // Create EvalPoints and send them to EvaluatorControl
    if (nullptr == _evaluatorControl)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, _step->getName() + ": EvaluatorControl not found");
    }

    // Currently, this method may be used inside an Iteration (Search or Poll, NM, ...),
    // or inside a MegaSearchPoll.
    auto iteration = dynamic_cast<const NOMAD::Iteration*>(_step->getParentOfType<NOMAD::Iteration*>());
    auto megaSearchPoll = dynamic_cast<const NOMAD::MegaSearchPoll*>(_step);
    
    if (useMesh && nullptr == iteration && nullptr == megaSearchPoll)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, _step->getName() + ": In keepPointsThatNeedEval: need a parent of type Iteration or MegaSearchPoll");
    }

    if (trialPoints.size() > 0)
    {
        _step->AddOutputInfo("Add points to eval queue for step " + _step->getName(), true, false);
        _step->AddOutputDebug(NOMAD::itos(trialPoints.size()) + " points to add to eval queue");
    }

    for (auto trialPoint : trialPoints)
    {
        // Try to insert the trialPoint in the cache.
        // If it is not in the cache, it will be added.
        // If the point is already in the cache, depending on its
        // EvalStatus we might want to evaluate it.

        // First, convert trial point to full dimension, since we are
        // now only working with the cache and the EvaluatorControl.
        auto trialPointSub = trialPoint;    // Used to get iteration
        trialPoint = trialPoint.makeFullSpacePointFromFixed(_step->getSubFixedVariable());

        // maxNumberEval is used to compute if we should re-evaluate this point.
        // Default value is 1.
        // This will be a parameter in the future. Currently not implemented.
        const int maxNumberEval = 1;
        bool doEval = NOMAD::CacheBase::getInstance()->smartInsert(trialPoint, maxNumberEval, _step->getEvalType());

        if (doEval)
        {
            NOMAD::EvalQueuePointPtr evalQueuePoint(new NOMAD::EvalQueuePoint(trialPoint));
            if (useMesh && nullptr == iteration)
            {
                iteration = megaSearchPoll->getIterForPoint(trialPointSub).get();
                if (nullptr == iteration)
                {
                    std::string s = _step->getName();
                    s += ": In keepPointsThatNeedEval: Could not determine iteration for point ";
                    s += trialPoint.display();
                    throw NOMAD::Exception(__FILE__,__LINE__, s);
                }
            }
            if ( useMesh )
            {
                auto mesh = iteration->getMesh();
                if ( mesh != nullptr )
                {
                    evalQueuePoint->setMeshSize(mesh->getdeltaMeshSize());
                    evalQueuePoint->setFrameSize(mesh->getDeltaFrameSize());
                    evalQueuePoint->setK(iteration->getK());
                }
            }
                

            evalQueuePoint->setComment(NOMAD::MainStep::getAlgoComment());
            evalQueuePoint->setGenStep(_step->getName());

            _evaluatorControl->addToQueue(evalQueuePoint);
            _step->AddOutputDebug("New point added to eval queue: " + trialPoint.display());
        }
        else
        {
            // Cache hit
            _step->AddOutputDebug("Point already found in cache: " + trialPoint.display());
        }
    }

    size_t evcNbPoints = _evaluatorControl->getQueueSize();
    if (evcNbPoints > 0)
    {
        // Not conditional to trialPoints size:
        // There could be leftover points to evaluate from previous evaluations.
        _step->AddOutputDebug("Eval queue now has " + NOMAD::itos(evcNbPoints) + " points.");
    }

    if (trialPoints.size() > 0)
    {
        _step->AddOutputInfo("Add points to eval queue for step " + _step->getName(), false, true);
    }

    NOMAD::OutputQueue::Flush();
}


void NOMAD::EvcInterface::setBarrier(const std::shared_ptr<NOMAD::Barrier>& subBarrier)
{
    if ( subBarrier == nullptr )
        return;

    // Input is the barrier from MegaIteration, which may belong to a subspace.
    // EvaluatorControl's barrier must be in full dimension.
    auto fixedVariable = _step->getSubFixedVariable();
    auto fullBarrier = std::make_shared<NOMAD::Barrier>(*subBarrier);

    // Clear xFeas and xInf lists and recompute them
    fullBarrier->clearXFeas();
    fullBarrier->clearXInf();
    for (auto xFeas : subBarrier->getAllXFeas())
    {
        auto xFeasFull = std::make_shared<NOMAD::EvalPoint>(xFeas->makeFullSpacePointFromFixed(fixedVariable));
        fullBarrier->addXFeas(xFeasFull, _step->getEvalType());
    }
    for (auto xInf : subBarrier->getAllXInf())
    {
        auto xInfFull = std::make_shared<NOMAD::EvalPoint>(xInf->makeFullSpacePointFromFixed(fixedVariable));
        fullBarrier->addXInf(xInfFull);
    }

    _evaluatorControl->setBarrier(fullBarrier);
}


// When points are generated and added to queue,
// we can start evaluation.
NOMAD::SuccessType NOMAD::EvcInterface::startEvaluation()
{
    _step->AddOutputInfo("Evaluate points for " + _step->getName(), true, false);

    NOMAD::SuccessType success = NOMAD::SuccessType::UNSUCCESSFUL;
    std::shared_ptr<NOMAD::AllStopReasons> stopReasons = _step->getAllStopReasons();

    if ( ! stopReasons->checkTerminate() )
    {
        // Evaluate points
        success = _evaluatorControl->run();
    }

    std::string s = _step->getName() + ": " + NOMAD::enumStr(success);
    s += ". Stop reasons: " + stopReasons->getStopReasonAsString() ;
    _step->AddOutputDebug(s);

    NOMAD::OutputQueue::Flush();

    _step->AddOutputInfo("Evaluate points for " + _step->getName(), false, true);

    return success;
}


bool NOMAD::EvcInterface::evalSinglePoint(NOMAD::EvalPoint &evalPoint, const NOMAD::Double &hMax)
{
    // Convert to full dimension before calling EvaluatorControl
    evalPoint = evalPoint.makeFullSpacePointFromFixed(_step->getSubFixedVariable());
    bool ret = _evaluatorControl->evalSinglePoint(evalPoint, hMax);
    // Convert back to subspace dimension
    evalPoint = evalPoint.makeSubSpacePointFromFixed(_step->getSubFixedVariable());

    return ret;
}
