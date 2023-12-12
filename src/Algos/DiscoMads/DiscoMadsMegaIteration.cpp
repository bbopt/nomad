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
#include "../../Algos/SubproblemManager.hpp"
#include "../../Algos/DiscoMads/DiscoMads.hpp"
#include "../../Algos/DiscoMads/DiscoMadsMegaIteration.hpp"
#include "../../Algos/DiscoMads/DiscoMadsUpdate.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"

#include "../../Cache/CacheBase.hpp"
#include "../../Cache/CacheSet.hpp"


bool NOMAD::DiscoMadsMegaIteration::discontinuityTest(const NOMAD::EvalPoint & x1, const NOMAD::EvalPoint & x2){
        // Retur True if (x1,x2) is a weak discontinuity is detected between x1 and x2 (called in evaluator callback)
        bool critvalue=false;

        // Do not compute criteria if the two points are equal as distance will be zero
        const NOMAD::Point * x1Point = x1.getX();
        const NOMAD::Point * x2Point = x2.getX();
        if(*x1Point==*x2Point)
            {
                return critvalue;
            }

        // If both points are correctly evaluated...
        if(x1.getEvalStatus(NOMAD::EvalType::BB)== NOMAD::EvalStatusType::EVAL_OK && x2.getEvalStatus(NOMAD::EvalType::BB)== NOMAD::EvalStatusType::EVAL_OK)
        {
            //As this function is used in eval callback, we use "Add" instead of "AddOutputInfo" to avoid segmentation faults as callback may be called deep in code
            OUTPUT_DEBUGDEBUG_START
                std::string s = "Test if revelation between points "+  std::to_string(x1.getTag()) + " and "+ std::to_string(x2.getTag());
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
                NOMAD::OutputQueue::Flush();
            OUTPUT_DEBUGDEBUG_END

            NOMAD::Double d=NOMAD::Point::dist(x1,x2);  // distance between 2 points

            // Check for numerical problem
            if(d==0)
            {
                OUTPUT_DEBUG_START
                std::string s = "Warning: DiscoMadsMegaIteration:: Revelation:: distance between tested points is null.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                NOMAD::OutputQueue::Flush();
                OUTPUT_DEBUG_END
                throw NOMAD::Exception(__FILE__,__LINE__,"Numerical precision problem");
            }

            // Revelation if distance between x1 and x2 < detectionRadius and rate of change of at least one revealing output exceeds the limit rate
            if (d < _detectionRadius)
            {
                auto arrayOutputx1 = x1.getEval(NOMAD::EvalType::BB)->getBBOutput().getBBOAsArrayOfDouble();  // BB output X1
                auto arrayOutputx2 = x2.getEval(NOMAD::EvalType::BB)->getBBOutput().getBBOAsArrayOfDouble();  // BB output x1

                // Loop on revealing outputs
                for(const int idxOutput: _idxRevealingOutput){
                    NOMAD::Double outputDiff = (arrayOutputx1[idxOutput]-arrayOutputx2[idxOutput]).abs();
                    if (outputDiff>_limitRate*d)
                    {
                        critvalue=true;
                        OUTPUT_DEBUGDEBUG_START
                        std::string s = "Revelation between points "+  std::to_string(x1.getTag()) + " and "+ std::to_string(x2.getTag());
                        s += " distance:  "+ d.display() + ", output variation: "+ outputDiff.display()+ ", output idx: "+std::to_string(idxOutput);
                        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                        NOMAD::OutputQueue::Flush();
                        OUTPUT_DEBUGDEBUG_END
                    }
                }
            }
        }
        return critvalue;
    }


bool NOMAD::DiscoMadsMegaIteration::proximityTestOnRevealingPoint(const NOMAD::Point & x1, const NOMAD::EvalPoint & x2){
    // Return true if dist(x1,x2) < exclusionRadius and x2 is a revealing point
    bool critvalue=false;

    // check first if point x2 is revealing
     if (x2.getEvalStatus(NOMAD::EvalType::BB)==NOMAD::EvalStatusType::EVAL_OK && x2.getRevealingStatus()>0)
     {
        // then compute distance
        NOMAD::Double d=NOMAD::Point::dist(x1,x2);  // distance between 2 points
        if(d< _exclusionRadius)
        {
            critvalue=true;
        }
     }

     //TODO: could be optimized by keeping track of a distance matrix for the whole run

    return critvalue;
}


void NOMAD::DiscoMadsMegaIteration::callbackCheckIfRevealingAndUpdate(NOMAD::EvalQueuePointPtr & evalQueuePoint)
{

    std::string s;     //for display

    // If evalQueuePoint was a MODEL or SURROGATE eval, it is useless to check for revealation as no new information on BB output is known
    if(evalQueuePoint->getEvalType()==NOMAD::EvalType::BB)
    {

        // Callback is called even after failed eval, discard treatment in this case //TODO: voir avec C.T. si on ne ferait pas plutôt le RunEvalCallback dans evalblock seulement si eval correcte
        if(!evalQueuePoint->isEvalOk(NOMAD::EvalType::BB))
        {
            return;
        }

        // Point violating EB constraints cannot be revealing
        if(!evalQueuePoint->isEBOk(NOMAD::EvalType::BB))
        {
            OUTPUT_DEBUG_START
                s = std::to_string(evalQueuePoint->getTag())+" is not a revealing point (violates at least one EB constraint).";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                NOMAD::OutputQueue::Flush();
            OUTPUT_DEBUG_END
            return;
        }

        auto cache = CacheBase::getInstance().get();
        std::vector<NOMAD::EvalPoint> revealingPointList;

        // 1) Try to reveal either hidden constraints or discontinuity
        // Reveal hidden constraints
        if(_detectHiddConst)
        {
            if(evalQueuePoint->getRevealingStatus()==2){
                evalQueuePoint->setRevealedConstraint(1.0);
                cache->update(*evalQueuePoint,NOMAD::EvalType::BB);  // update revealing status and revealed constraint
                revealingPointList.push_back(*evalQueuePoint);
            }
        }
        // or reveal discontinuities
        else{
            // Revelation test between recently evaluated point evalQueuePoint and cache points
            auto crittest = [&](const EvalPoint& x2){return this->discontinuityTest(*evalQueuePoint, x2);};
            cache->find(crittest,revealingPointList);   // NB: revealingPointList does not contain evalQueuePoint

            // If we have detected at least one revealing point thanks to evalQueuePoint...
            if(revealingPointList.size()>0)
            {
                // ...then evalQueuePoint is a new revealing point and shoul be updated
                OUTPUT_DEBUG_START    // NOTE: safer to use "Add" instead of "AddOutputInfo" in callback to avoid segmentation faults when accessing step name
                    s = std::to_string(evalQueuePoint->getTag())+" is a revealing point.";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    NOMAD::OutputQueue::Flush();
                OUTPUT_DEBUG_END

                evalQueuePoint->setRevealingStatus(2);           // set to 2 to indicate that the point was found revealing at this iteration
                evalQueuePoint->setRevealedConstraint(1.0);
                cache->update(*evalQueuePoint,NOMAD::EvalType::BB);   // update revealing status and revealed constraint

                // update revealing status of other new revealing points (their revealed constraints are set in DiscoMadsBarrier)
                for(auto evalPoint : revealingPointList)
                {
                    if(evalPoint.getRevealingStatus()==0)         // discard revealing points already known (with status ==1)
                    {
                        evalPoint.setRevealingStatus(2);
                        cache->update(evalPoint,NOMAD::EvalType::BB);   // update only revealing status
                        // DEV : va poser problème avec le parallelisme
                    }
                }

                // add evalQueuePoint to the list of new revealing points for display
                revealingPointList.push_back(*evalQueuePoint);
            }

        }


        // 2) If evalQueuePoint is not a revealing point, we check if it is close to revealing points
        if(revealingPointList.size()==0){
            OUTPUT_DEBUG_START
                std::string s = std::to_string(evalQueuePoint->getTag())+" is not a revealing point.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                NOMAD::OutputQueue::Flush();
            OUTPUT_DEBUG_END

            // locate revealing neighbours points
            std::vector<NOMAD::EvalPoint> revealingNeighbours;
            auto crittest = [&](const EvalPoint& x2){return this->proximityTestOnRevealingPoint(*evalQueuePoint, x2);};
            cache->find(crittest,revealingNeighbours);
    
            if(revealingNeighbours.size()>0)
            {
                NOMAD::Double revealedConstraint = evalQueuePoint->getRevealedConstraint();

                // compute arfificial constraint of evalQueuePoint
                for(auto revealingPoint : revealingNeighbours)
                {
                    OUTPUT_DEBUG_START
                    std::string s = std::to_string(revealingPoint.getTag())+" is a revealing point close to "+std::to_string(evalQueuePoint->getTag());
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    NOMAD::OutputQueue::Flush();
                    OUTPUT_DEBUG_END

                    NOMAD::Double d=NOMAD::Point::dist(*evalQueuePoint,revealingPoint); // Distance //TODO: should be optimized with a distance matrix (computed here for the 3rd time)
                    NOMAD::Double tmpConstraint =1-d/_exclusionRadius;
                    if(tmpConstraint>revealedConstraint)
                    {
                        revealedConstraint=tmpConstraint;
                    }

                }

                // Update value of revealed constraint
                if (revealedConstraint> evalQueuePoint->getRevealedConstraint())
                {
                    evalQueuePoint->setRevealedConstraint(revealedConstraint);
                    cache->update(*evalQueuePoint,NOMAD::EvalType::BB);
                }
            }
        }

        // Debug output to check that all new revealing points have been detected
        OUTPUT_DEBUG_START
        std::vector<NOMAD::EvalPoint> newRevealingPointsBis;
        auto crittestNewRevealing = [&](const EvalPoint& x){return x.getRevealingStatus()==2;};
        cache->find(crittestNewRevealing,newRevealingPointsBis);
        if (newRevealingPointsBis.size()>0){
            std::string s;
            s = "List of new revealing points (" +std::to_string(newRevealingPointsBis.size())+"):";
            for(auto point: newRevealingPointsBis)
            {
                s += " "+std::to_string(point.getTag());
            }
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            NOMAD::OutputQueue::Flush();
        }
        OUTPUT_DEBUG_END

        // Export cache at the end of megaiteration (for debug)
        bool cacheExport = false;
        if(cacheExport)
        {
            // Create special cache file (tmp folder should exist)
            std::string specialCacheFile ="tmp/cache_dynamique_"+std::to_string(evalQueuePoint->getTag())+".txt";
            exportCache(specialCacheFile);
        }
    }
}

// Only used for advanced debug
void NOMAD::DiscoMadsMegaIteration::exportCache(std::string cacheFile)
{

    NOMAD::EvalPointPtr refBestInf = nullptr, refBestFeas = nullptr;
    // Best points at beginning of iteration
    if(_barrier!=nullptr)
    {
        refBestInf=_barrier->getRefBestInf();      // NB: creates segmentation fault if called in init because there is no megaIteration ancestor
        refBestFeas = _barrier->getRefBestFeas();
    }

    // get all cache points
    auto cache = CacheBase::getInstance().get();
    std::vector<NOMAD::EvalPoint> cachePoints;
    cache->getAllPoints(cachePoints);

    // write cache file
    ofstream myCacheFile;
    myCacheFile.open (cacheFile);
    for (const auto & evalPoint : cachePoints)
    {
        if (nullptr != evalPoint.getEval(NOMAD::EvalType::BB) && evalPoint.getEval(NOMAD::EvalType::BB)->goodForCacheFile())
        {
            // Tag
            myCacheFile<<evalPoint.getTag()<<" ";

            // BB output
            myCacheFile<<evalPoint.getBBO(NOMAD::EvalType::BB)<<" ";

            // f value
            myCacheFile<<evalPoint.getF()<<" ";

            // h value
            myCacheFile<<evalPoint.getH()<<" ";
            
            // Revealing status
            myCacheFile<<evalPoint.getRevealingStatus()<<" ";

            // Is this a best feasible?
            int solutionStatus = 0;
            if (refBestFeas!=nullptr && *refBestFeas==evalPoint) {
                solutionStatus=2;
            }
            // Is this a best infeasible?
            else if (refBestInf!=nullptr && *refBestInf==evalPoint){
                solutionStatus=1;
            }
            myCacheFile<<solutionStatus;

            myCacheFile<<std::endl;
        }
    }

}




void NOMAD::DiscoMadsMegaIteration::callbackEvalOpportStop(bool &opportunisticIterStop, NOMAD::EvalQueuePointPtr & evalQueuePoint)
{
    
    // If evalQueuePoint was a MODEL or SURROGATE eval, it is useless to check for revealation as no new information on BB output is known
    if(evalQueuePoint->getEvalType()==NOMAD::EvalType::BB)
    {
        // Check if the point is revealing
        if(evalQueuePoint->getRevealingStatus()==2)
        {
            _isRevealing = true;           //then iteration is revealing
            opportunisticIterStop = true;      //and evaluations at this iteration should be stopped.
        }
    }
}



void NOMAD::DiscoMadsMegaIteration::callbackFailedEval(EvalQueuePointPtr & evalQueuePoint)
{
    // NB: Only used when discomads is used to reveal hidden constraints regions
    if (nullptr != evalQueuePoint && evalQueuePoint->getEvalType()==NOMAD::EvalType::BB)
    {
        auto eval = evalQueuePoint->getEval(NOMAD::EvalType::BB);
        if ( nullptr != eval )
        {
            NOMAD::Double highValue = _hiddConstOutputValue;     // value used in paper: 1e+20

            // Access to list of BB output types
            auto bbOutputTypeList = NOMAD::EvcInterface::getEvaluatorControl()->getCurrentBBOutputTypeList();

            // Create new values of bb Output
            auto bboList = evalQueuePoint->getEval(NOMAD::EvalType::BB)->getBBOutput().getBBOAsArrayOfDouble();
            string bboutput;
            for(size_t i = 0 ; i < bbOutputTypeList.size(); ++i)
            {
                NOMAD::BBOutputType bbot = bbOutputTypeList[i];
                // Set FOBJ and PB constraints to high value
                if (bbot.isObjective() || bbot==NOMAD::BBOutputType::PB)
                {
                    bboutput+=highValue.tostring()+" ";
                }
                else if(bbot==NOMAD::BBOutputType::EB){
                // Set EB constraints to 0 so that they are satisfied
                    bboutput+=std::to_string(0.0)+" ";
                }
                else if(bbot==NOMAD::BBOutputType::RPB){
                // Revealed constraint is always the last, it may be skipped (set automatically in setBBO)
                    continue;
                }
                else{
                    // Check by security : other types of constraints are not treated
                    throw NOMAD::Exception(__FILE__,__LINE__,"Discomads for hidden constraints: callback for failed eval only treat OBJ/PB/EB/RPB constraints.");
                }
            }

            // Assign new bb output to point (change eval status too)
            bool eval_ok = true;
            eval->setBBO(bboutput,bbOutputTypeList,eval_ok);
            
            //  Update revealing status of point
            evalQueuePoint->setRevealingStatus(2);
            auto cache = CacheBase::getInstance().get();
            cache->update(*evalQueuePoint,NOMAD::EvalType::BB);
        }
    }
}



void NOMAD::DiscoMadsMegaIteration::callbackPostProcessing(const NOMAD::Step & step, bool &stop)
{
    // Treat special case of revealation during a search algo: in this case remaining evaluations for this search algo should be stopped
    // as well as parent search evaluations for this DiscoMads iteration
    // NB: this situation is not taken into account in IterationUtils::updateStopReasonForIterStop)

    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    stop = false;  // reset as we don't control a global stop in this callback
    
    // This is postprocessing for BB only
    if (NOMAD::EvalType::BB != evc->getCurrentEvalType())
    {
        return;
    }
    auto evcStopReason = evc->getStopReason(-1);
    
    // If there was a revelation, stop type of evaluator was changed to opportunistic
    if (evcStopReason.checkStopType(NOMAD::EvalMainThreadStopType::CUSTOM_OPPORTUNISTIC_ITER_STOP))
    {

         // Is this step done during a search ?
        NOMAD::Search * searchStep = step.getParentOfType<NOMAD::Search*>(false);
        if(nullptr!= searchStep)
        {
            // Is this done during a search algo ? Look for an algorithm among the parents. It should be a search algorithm, that is, not the root Algorithm. Stop it completely if found.
            // If we simply search for an algorithm among the parents we may endup with the root Mads algorithm.
            auto algoSM = step.getFirstAlgorithm();
            if (algoSM != step.getRootAlgorithm())
            {
                OUTPUT_DEBUG_START
                // NOTE: it is safer to use "Add" instead of "AddOutputInfo" in this callback to avoid segmentation faults as they may be called deep in code
                NOMAD::OutputQueue::Add("User stop of the search algo "+algoSM->getName(), NOMAD::OutputLevel::LEVEL_DEBUG);
                NOMAD::OutputQueue::Flush();
                OUTPUT_DEBUG_END

                // stop the parent search step iteration (may contains severeal searches)
                searchStep->getAllStopReasons()->set(NOMAD::IterStopType::USER_ITER_STOP);
                // stop the search algo used in search step
                algoSM->getAllStopReasons()->set(NOMAD::IterStopType::USER_ALGO_STOP);

                // required to not erase USER_ALGO_STOP by USER_ITER_STOP in IterationUtils::updateStopReasonForIterStop
                evc->setStopReason(-1, NOMAD::EvalMainThreadStopType::STARTED);
            }
        }

        // Is this step done during a revealing poll or a poll
        // => stop is managed by IterationUtils::updateStopReasonForIterStop and passed to checkTerminate
    }
}



void NOMAD::DiscoMadsMegaIteration::init()
{
    setStepType(NOMAD::StepType::MEGA_ITERATION);

    // Check computeType (DiscoMadsMegaIteration should not be done within PhaseOne)
    if(NOMAD::EvcInterface::getEvaluatorControl()->getComputeType()!=NOMAD::ComputeType::STANDARD)
    {
        string s = "DiscoMadsMegaIteration: Only STANDARD compute type is handled";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }
 
    // Get values of discomads parameters


    _detectionRadius = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_DETECTION_RADIUS");   // only for discontinuity revealation
    _limitRate = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_LIMIT_RATE");               // only for discontinuity revealation
    _exclusionRadius  = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_EXCLUSION_RADIUS");


    _detectHiddConst = _runParams->getAttributeValue<bool>("DISCO_MADS_HID_CONST");                   // only for hidden constraints revealation
    _hiddConstOutputValue = _runParams->getAttributeValue<NOMAD::Double>("DISCO_MADS_HID_CONST_OUTPUT_VALUE");               // only for hidden constraints revealation




    // Build vector of indices of revealing output
    const auto bbotList = NOMAD::Algorithm::getBbOutputType();
    std::vector<int> revealingOutputIdx;
        for(unsigned int idxOutput=0; idxOutput<bbotList.size();idxOutput++)
        {
            if(bbotList[idxOutput].isRevealing())
            {
                revealingOutputIdx.push_back(idxOutput);
            }
        }
    _idxRevealingOutput = revealingOutputIdx;
    
    // Set evaluator control callbacks (revelation of discontinuities OR hidden constraints and exclusion)

        // Callback to check if the point that has just been evaluated is revealing and update its revealed constraint
    NOMAD::EvalCallbackFunc<NOMAD::CallbackType::EVAL_UPDATE> cbInterEvalUpdate = [&](EvalQueuePointPtr & evalQueuePoint){return this->callbackCheckIfRevealingAndUpdate(evalQueuePoint);};
    NOMAD::EvcInterface::getEvaluatorControl()->addEvalCallback<NOMAD::CallbackType::EVAL_UPDATE>(cbInterEvalUpdate);

        // Callback during eval to trigger an iteration opportunistic stop (that will cancel all remaining evaluations for this iteration) if a revealing point has been evaluated
    NOMAD::EvalCallbackFunc<NOMAD::CallbackType::EVAL_OPPORTUNISTIC_CHECK> cbInter = [&](EvalQueuePointPtr & evalQueuePoint, bool &opportunisticEvalStop, bool &opportunisticIterStop){opportunisticEvalStop = false; return this->callbackEvalOpportStop(opportunisticIterStop, evalQueuePoint);};
    NOMAD::EvcInterface::getEvaluatorControl()->addEvalCallback<NOMAD::CallbackType::EVAL_OPPORTUNISTIC_CHECK>(cbInter);
    
        // Callback called at postprocessing to specifically stop a search algo (e.g Nelder-Mead) if a revealation occured during this search algo
    auto cbInterPostProcessing= [&](const NOMAD::Step & step, bool &stop){return this->callbackPostProcessing(step,stop);};
    this->addCallback(NOMAD::CallbackType::POSTPROCESSING_CHECK, cbInterPostProcessing);

    if(_detectHiddConst)
    {
        // Callback specific to hidden constraints detection to put high value of f to points for which evaluation failed
        NOMAD::EvalCallbackFunc<NOMAD::CallbackType::EVAL_FAIL_CHECK> cbInterFailedEval = [&](EvalQueuePointPtr & evalQueuePoint){return this->callbackFailedEval(evalQueuePoint);};
        NOMAD::EvcInterface::getEvaluatorControl()->addEvalCallback<NOMAD::CallbackType::EVAL_FAIL_CHECK>(cbInterFailedEval);
    }

}

void NOMAD::DiscoMadsMegaIteration::startImp()
{
    // Update main mesh and barrier.
    NOMAD::DiscoMadsUpdate update( this );
    update.start();
    update.run();
    update.end();

    // Verify mesh stop conditions.
    _mainMesh->checkMeshForStopping(_stopReasons);

    OUTPUT_DEBUG_START
    AddOutputDebug("Mesh Stop Reason: " + _stopReasons->getStopReasonAsString());
    OUTPUT_DEBUG_END
}



bool NOMAD::DiscoMadsMegaIteration::runImp()
{
    std::string s;
    if ( _stopReasons->checkTerminate() )
    {
        OUTPUT_DEBUG_START
        s = "MegaIteration: stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
        return false;
    }

    // Get DiscoMads ancestor to call terminate(k)
    NOMAD::DiscoMads* discomads = getParentOfType<NOMAD::DiscoMads*>();
    if (nullptr == discomads)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "DiscoMads Iteration without DiscoMads ancestor");
    }

    // Until termination criterion is reached
    if (!discomads->terminate(_k))
    {
        // Reset _isRevealing for this iteration
        _isRevealing = false;
        
        OUTPUT_DEBUG_START
        AddOutputDebug("Iteration generated:");
        AddOutputDebug(_madsIteration->getName());
        NOMAD::ArrayOfDouble meshSize  = _madsIteration->getMesh()->getdeltaMeshSize();
        NOMAD::ArrayOfDouble frameSize = _madsIteration->getMesh()->getDeltaFrameSize();
        AddOutputDebug("Mesh size:  " + meshSize.display());
        AddOutputDebug("Frame size: " + frameSize.display());
        OUTPUT_DEBUG_END

        _madsIteration->start();
        bool iterSuccessful = _madsIteration->run();
        _madsIteration->end();

        if (iterSuccessful)
        {
            OUTPUT_DEBUG_START
            s = getName() + ": new success " + NOMAD::enumStr(_success);
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }

        // Update MegaIteration's stop reason
        if (_stopReasons->checkTerminate())
        {
            OUTPUT_DEBUG_START
            s = getName() + " stop reason set to: " + _stopReasons->getStopReasonAsString();
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }

        // Note: Delta (frame size) will be updated in the Update step next time it is called.
        if (getUserInterrupt())
        {
            hotRestartOnUserInterrupt();
        }
    }

    // MegaIteration is a success if either a better xFeas or
    // a dominating or partial success for xInf was found.
    // See Algorithm 12.2 from DFBO.

    // return true if we have a partial or full success.
    return (_success >= NOMAD::SuccessType::PARTIAL_SUCCESS);

}


void NOMAD::DiscoMadsMegaIteration::endImp()
{
    // Increment number of revealing iteration
    if(_isRevealing)
    {
        NOMAD::EvcInterface::getEvaluatorControl()->incrementNbRevealingIter();
    }

    // Call default endImp
    NOMAD::MegaIteration::endImp();

}
