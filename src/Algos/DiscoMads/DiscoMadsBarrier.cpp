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
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Algos/DiscoMads/DiscoMadsBarrier.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Output/OutputDirectToFile.hpp"
#include "../../Algos/EvcInterface.hpp"


bool NOMAD::DiscoMadsBarrier::proximityTest(const NOMAD::Point & x1, const NOMAD::EvalPoint & x2){
    // Return true if dist(x1,x2) < exclusionRadius and x2 has been correctly evaluated
    bool critvalue=false;
    NOMAD::Double d=NOMAD::Point::dist(x1,x2);  // distance between 2 points
     if (x2.getEvalStatus(NOMAD::EvalType::BB)==NOMAD::EvalStatusType::EVAL_OK && d< _exclusionRadius)
     {
        critvalue=true;
     } 

    return critvalue;
}



const NOMAD::Double NOMAD::DiscoMadsBarrier::getKiemeHvalue(const std::vector<EvalPointPtr>& evalPointList, const size_t k, NOMAD::EvalType evalType) const
{
    NOMAD::Double kiemeHValue = NaN;
    std::string s; // for info output

    if (evalPointList.size() == 0)
    {
        OUTPUT_INFO_START
        s = "Warning: DiscoMadsBarrier::getKiemeHvalue called on an empty evalPoints list";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
        NOMAD::OutputQueue::Flush();
        OUTPUT_INFO_END
        return kiemeHValue;
    }

    // Get hvalues of all points
    std::vector<NOMAD::Double> hValues;
    for(auto it:evalPointList)
    {   
        hValues.push_back(it->getEval(evalType)->getH());
    }
    if(evalPointList.size()!= hValues.size())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Problem to compute k-ieme h value of the evalPoint list.");
    }

    // Order points with increasing h values
    std::vector<EvalPointPtr> sortedEvalPointList = evalPointList;
    std::sort(sortedEvalPointList.begin(), sortedEvalPointList.end(), [](const NOMAD::EvalPointPtr& x1,const NOMAD::EvalPointPtr& x2) {
      return x1->getEval(NOMAD::EvalType::BB)->getH() < x2->getEval(NOMAD::EvalType::BB)->getH() ;
   });

    // Get k-ieme value
    kiemeHValue = sortedEvalPointList.at(k-1)->getEval(NOMAD::EvalType::BB)->getH();

    OUTPUT_DEBUG_START
    s = "List of evalPoints sorted with increasing h values: \n ";
    for(auto point: sortedEvalPointList)
    {   
        s +=point->display()+"\n";
    }
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    s = "h-value of k-ieme point (k= "+std::to_string(k) +"): "+kiemeHValue.tostring();
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    NOMAD::OutputQueue::Flush();
    OUTPUT_DEBUG_END

    return kiemeHValue;
}


const size_t NOMAD::DiscoMadsBarrier::getNonDominatedInfPoints(std::vector<NOMAD::EvalPointPtr> &evalPointList,const NOMAD::EvalType evalType,const NOMAD::ComputeType computeType) 
{
    evalPointList.clear();

    if (_xInf.size()> 0)
    {
        std::vector<NOMAD::EvalPointPtr>::const_iterator it1= _xInf.begin(), it2;
        while ( it1 !=_xInf.end() )
            {
                // Do not consider points with h>hMax
                if ((*it1)->getEval(evalType)->getH(computeType) > _hMax)
                {
                    continue;
                }

                // Find all infeasible non-dominated points
                bool isDominated = false;
                for ( it2 = _xInf.begin(); it2 < _xInf.end(); it2++)
                {
                    if (it1 == it2)
                    {
                        continue;
                    }
                // Case ep2 dominates ep1 --> ep1 is dominated
                    if ((*it2)->dominates(*(*it1), evalType, computeType))
                    {
                        isDominated = true;
                        break;
                    }
                }

                // No point dominates ep1
                if (! isDominated)
                {   
                    evalPointList.push_back(*it1);
                }
                it1++;    
            }
    }
    return evalPointList.size();
};



// Points from evalPointList are already in subproblem dimension.
bool NOMAD::DiscoMadsBarrier::updateWithPoints(const std::vector<NOMAD::EvalPoint>& evalPointList,
                                      NOMAD::EvalType evalType,
                                      NOMAD::ComputeType computeType,
                                      const bool keepAllPoints /* Not used here */,
                                      const bool updateIncumbentsAndHmax )
{

    // Just in case
    if(evalType==NOMAD::EvalType::MODEL)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "DiscoMAdsBarrier:: shoudl not be used on quadratic model optimization because it may be very slow.");
    }

    bool updatedInc = false;
    bool revealingIteration = false;
    string s; // for output
    
    // Detect if there has been a revelation (so the iteration is revealing) 
    std::vector<NOMAD::EvalPoint> newRevealingPoints; 
    auto cache = CacheBase::getInstance().get();
    auto crittestNewRevealing = [&](const EvalPoint& x){ return x.getRevealingStatus()==2;};
    cache->find(crittestNewRevealing,newRevealingPoints);

    if (newRevealingPoints.size()>0){
        revealingIteration=true;

        OUTPUT_DEBUG_START
        s = "Iteration is revealing; tags of "+ std::to_string(newRevealingPoints.size()) +" new revealing points seen by the DiscoMadsBarrier: ";
        for(auto point: newRevealingPoints)
        {   
            s+= std::to_string(point.getTag())+" ";
        }
        s+= "\n";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        NOMAD::OutputQueue::Flush();
        OUTPUT_DEBUG_END
        }

    // If the iteration is revealing, special update of the barrier 
    if (revealingIteration){
        // In case of revelation, evaluations at this iteration are stopped. 
        // Then do the update of incumbents and hmax reggardless of updateIncumbentsAndHmax, which may be false if revelation occurs during a search for example
        //(_updateIncumbentsAndHmax is set to True in search only if full success)

        // --- Step 1 : store information of current barrier
        Double hMaxBeforeUpdate = _hMax;

        // Non dominated infeasible points with h<hmax
        size_t nbNonDomInfPtsBelowHmaxBeforeUpdate = 0;
        std::vector<EvalPointPtr> nonDomInfPtsBelowHmaxBeforeUpdate;     
        nbNonDomInfPtsBelowHmaxBeforeUpdate = getNonDominatedInfPoints(nonDomInfPtsBelowHmaxBeforeUpdate, evalType, computeType);

        // --- Step 2 : update RPB constraints considering revealed information
        for (auto & revealingPoint : newRevealingPoints)
        {// For each new revealing point...
            // set revealed constraint to 1.0 (only the constraint of the current evaluated point was updated in callbackCheckIfRevealingAndUpdate)
            const NOMAD::Double valeur=revealingPoint.getRevealedConstraint();
            revealingPoint.setRevealedConstraint(1.0);
        
            // reset revealing status to 1
            revealingPoint.setRevealingStatus(1);
            cache->update(revealingPoint,NOMAD::EvalType::BB);

            // locate neighbours points around the new revealing point
            std::vector<NOMAD::EvalPoint> evalPointToUpdate; 
            auto crittest = [&](const EvalPoint& x2){return this->proximityTest(revealingPoint, x2);};
            cache->find(crittest,evalPointToUpdate);

            OUTPUT_DEBUG_START
                s = "Points close to revealing point "+ std::to_string(revealingPoint.getTag())+": ";
                for(auto point: evalPointToUpdate)
                {   
                    s+= std::to_string(point.getTag())+" ";
                }
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END

            // loop on neighbours to update their RPB constraint
            for (auto & ep : evalPointToUpdate)
            {
                //do not consider update of points with respect to itself, useless as revealing point has RPB = 1.0
                if(revealingPoint==ep){
                    continue;
                }

                //do not update revealed constraint if point violate EB constraints
                if(!ep.isEBOk(NOMAD::EvalType::BB))
                {
                    OUTPUT_DEBUG_START
                    s= "\t RPB constraint of point "+std::to_string(ep.getTag())+" not updated (point violates at least one EB constraint).";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                    continue;
                }

                // compute revealed constraint
                NOMAD::Double d=NOMAD::Point::dist(revealingPoint,ep); 
                NOMAD::Double value_constraint =1-d/_exclusionRadius;

                // update revealed constraint only if increasing
                if(ep.getRevealedConstraint()<value_constraint)
                {
                    OUTPUT_DEBUG_START
                    s= "\t RPB constraint of point " +std::to_string(ep.getTag())+ " updated from "+ep.getRevealedConstraint().tostring()+"to"+value_constraint.tostring();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END

                    ep.setRevealedConstraint(value_constraint);
                    cache->update(ep,NOMAD::EvalType::BB);
                }
                else
                {
                    OUTPUT_DEBUG_START
                    s= "\t RPB constraint of point " +std::to_string(ep.getTag())+ " not updated";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                }
            }
        }

        // --- Step 3 : update points in the barrier
        // access all infeasible points in cache
        std::vector<EvalPoint> cachePoints;
        cache->getAllPoints(cachePoints);

        // reset xFeas, xInf and hmax of the barrier
        clearXFeas();
        clearXInf();
        setHMax(INF);

        // reinsert all points (no points discarded as hmax = Inf)
        ProgressiveBarrier::updateWithPoints(cachePoints, evalType, computeType, true, true /*true update infeasible incumbent and hmax*/ );  
            // Important : hmax stays at INF cause barrier was reset just before

        // get number of infeasible non dominated points(these points may have changed because of increase in RPB constraint)
        size_t nbNonDomInfPts =0;
        std::vector<EvalPointPtr> nonDomInfPoints;    // U^{k+1} in paper  
        nbNonDomInfPts = getNonDominatedInfPoints(nonDomInfPoints, evalType, computeType);   // no points discarded as hmax = Inf

        // determine number of infeasible points to keep in the barrier (p. 1850 of paper) // H1 : faire test unitaire sur ça
        size_t nbInfPointsToKeep;                        // N(k+,h(bar(x))) in paper
        if (nbNonDomInfPtsBelowHmaxBeforeUpdate==0)    // if N(k,h^k_max)==0
        {   // N(k+1,h(bar(x)))= U^{k+1} 
            nbInfPointsToKeep = nbNonDomInfPts;        
        }
        else
        {   //N(k+1,h(bar(x)))=       min(N(k,h^k_max),U^{k+1}) 
            nbInfPointsToKeep = std::min(nbNonDomInfPtsBelowHmaxBeforeUpdate,nbNonDomInfPts); 
        }
            
        OUTPUT_DEBUG_START
        s= "Number of non dominated infeasible points below hmax at beginning of iteration: "+std::to_string(nbNonDomInfPtsBelowHmaxBeforeUpdate);
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        s= "Number of non dominated infeasible points regarless hmax: "+std::to_string(nbNonDomInfPts);
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        s= "Number of non dominated infeasible points to keep below hmax for next iteration: " +std::to_string(nbInfPointsToKeep);
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END


        // Set hMax to the hvalue of the nbInfPointsToKeep-eme infeasible undominated point of the barrier (bar(x) in paper)
        NOMAD::Double hValueToSet = getKiemeHvalue(nonDomInfPoints, nbInfPointsToKeep, evalType); // H1 do print for debug
        setHMax(hValueToSet);

        OUTPUT_DEBUG_START
        s= "hmax updated from "+hMaxBeforeUpdate.tostring()+" to "+hValueToSet.tostring();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END

        NOMAD::OutputQueue::Flush();
        // Recompute infeasible incumbents with the new hmax
        setInfeasibleIncumbents(evalType, computeType);

        // --- Step 4: update solution file, which was written before update of revealed constraints of neighbours points 
        // so best may have changed
        OUTPUT_DIRECTTOFILE_START
        // Access best feas
        std::vector<EvalPointPtr> bestFeasPoints;
        bestFeasPoints = getAllXIncFeas();
        NOMAD::EvalPointPtr BestFeas = bestFeasPoints[0];
    
        // Write in solution file
        NOMAD::StatsInfo info;
        if(BestFeas!=nullptr)
        { 
            info.setBBO(BestFeas->getBBO(NOMAD::EvalType::BB));
            info.setSol(*(BestFeas->getX()));   
        }
        else{
            // There may be no best feasible as feasibility of points may change during run
            info.setBBO("No best feasible solution at this iteration.");
        }
        NOMAD::OutputDirectToFile::Write(info, true);
        OUTPUT_DIRECTTOFILE_END
    }

    // If iteration is not revealing, do usual update from progressive barrier
    else
    {
        updatedInc = ProgressiveBarrier::updateWithPoints(evalPointList,
                                                        evalType,
                                                        computeType,
                                                        keepAllPoints /* Not used here */,
                                                        updateIncumbentsAndHmax);
        }
      
    return updatedInc;
}
