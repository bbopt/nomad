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

#include "../../Eval/ComputeSuccessType.hpp"
#include "../SimpleMads/SimpleProgressiveBarrier.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::SimpleProgressiveBarrier::init(const NOMAD::Point& fixedVariable,
                                           const std::vector<NOMAD::SimpleEvalPoint>& evalPointList)
{

    if (fixedVariable.isEmpty())
    {
        std::string s = "Error: Fixed variable of dimension 0";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
    
    // Constructor's call to update should not update ref best points.
    updateWithPoints(evalPointList);

}


void NOMAD::SimpleProgressiveBarrier::updateRefBests()
{
    _refBestFeas = getCurrentIncumbentFeas();
    _refBestInf  = getCurrentIncumbentInf();
}

// This update of the simple progressive barrier keeps only at most two undominated points with h<=hMax
bool NOMAD::SimpleProgressiveBarrier::updateWithPointsKeep2(const std::vector<NOMAD::SimpleEvalPoint>& evalPointList)
{

    bool updated = false;
    bool updatedFeas = false;
    bool updatedInf = false;
    bool updatedInc = false;
    bool updatedIncFeas = false;
    bool updatedIncInf = false;
    bool updatedHMax = false;

    _xTmpInf.clear();

    std::string s;  // for output info

    // Loop to identify the feasible barrier point. Will become the feasible incumbent later.
    for (const auto & evalPoint : evalPointList)
    {
        if (evalPoint.getH()==0)
        {
            const NOMAD::Double& f = evalPoint.getF();
            if ( !f.isDefined() || f == NOMAD::INF)
            {
                continue;
            }
            if ( ! _xFeas.isDefined() || f < _xFeas.getF())
            {
                _xFeas = evalPoint;
                updatedFeas = true;
            }
        }
    }

    // Second loop updates the barrier infeasible points.
    for (const auto & evalPoint : evalPointList)
    {

        const NOMAD::Double& h = evalPoint.getH();
        // Consider infeasible points with h < INF
        // Do not consider points for which h > hmax
        if (!h.isDefined() || (h == NOMAD::INF) ||
            ((_hMax < NOMAD::INF) && (h > _hMax)))
        {
            continue;
        }
        else if (h > 0)
        {
            updatedInf = true;
            _xTmpInf.push_back(evalPoint);
        }
    }

    if (updatedInf || updatedFeas)
    {
        _incumbentsAndHMaxUpToDate = false;
    }

    //
    // Ready to update hMax and the incumbents (if requested)
    NOMAD::Double hMaxPrev = _hMax;
    if (!_incumbentsAndHMaxUpToDate)
    {
        NOMAD::SuccessType feasSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
        
        // Step One: Set xIncFeas (feasible incumbent) from xFeas
        if (_xFeas.isDefined())
        {
            NOMAD::Double f = _xFeas.getF();
            NOMAD::Double fMinRef=_xIncFeas.getF(); // fMinRef not defined if _xIncFeas not defined
            if ( ! _xIncFeas.isDefined() || f < fMinRef)
            {
                _xIncFeas = _xFeas;
                updatedIncFeas = true;
                feasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
            }
        }

        // Step two: Detect the infeasible success type
        NOMAD::SuccessType infeasSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
        if ( !_xIncInf.isDefined() && !_xTmpInf.empty())
        {
            // Case with no previous incumbent (hence, no barrier point).
            infeasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else if (!_xTmpInf.empty())
        {
            // Note: An improving point (with respect to the infeasible incumbent) may have been added at the current iteration  -> loop on _xTmpInf.
            std::vector<NOMAD::SimpleEvalPoint>::iterator it;
            for ( it = _xTmpInf.begin(); it < _xTmpInf.end(); it++ )
            {
                NOMAD::SuccessType successType = computeSuccessType( *it, _xIncInf);
                if (successType > infeasSuccessType)
                {
                    infeasSuccessType = successType;
                }
                if (infeasSuccessType == NOMAD::SuccessType::FULL_SUCCESS)
                {
                    break;
                }
            }
        }
        
        // Step three: remove dominated points
        updateBarrierPointsKeep2();
        
//        // Temp for testing
//        std::cout<<"================================="<<std::endl;
//        for(const NOMAD::SimpleEvalPoint & x: _xInf)
//        {
//            std::cout<<"F="<<x.getF()<<" H="<<x.getH()<<std::endl;
//        }

        // Step four: update hMax
        if (infeasSuccessType == NOMAD::SuccessType::FULL_SUCCESS
            || feasSuccessType == NOMAD::SuccessType::FULL_SUCCESS)
        {
            // We have a new incumbent: infeasSuccessType or feasSuccessType is FULL_SUCCESS
            if ( _xIncInf.isDefined())
            {
                _hMax = _xInf[0].getH();
                _xIncInf = _xInf[0];
            }
            // Else. Case with no prior infeasible incumbent: we do not update hMax. Typically, hMax = INF.
        }
        else
        {
            _hMax = _xInf[0].getH();
            _xIncInf = _xInf[0];
        }

        if (_hMax < hMaxPrev)
        {
            updatedHMax = true;
        }

        _incumbentsAndHMaxUpToDate = true; // For the next update
    }

    updated = updatedFeas || updatedInf;
    updatedInc = updatedIncFeas || updatedIncInf;

    return updatedInc;
}

// Points from evalPointList are already in subproblem dimension.
bool NOMAD::SimpleProgressiveBarrier::updateWithPoints(const std::vector<NOMAD::SimpleEvalPoint>& evalPointList)
{

    bool updated = false;
    bool updatedFeas = false;
    bool updatedInf = false;
    bool updatedInc = false;
    bool updatedIncFeas = false;
    bool updatedIncInf = false;
    bool updatedHMax = false;

    _xTmpInf.clear();

    // The simple progressive barrier keeps only the undominated points with h<=hMax
    
    std::string s;  // for output info

    // Loop to identify the feasible barrier point. Will become the feasible incumbent later.
    for (const auto & evalPoint : evalPointList)
    {

        if (evalPoint.getH()==0)
        {
            const NOMAD::Double& f = evalPoint.getF();
            if ( !f.isDefined() || f == NOMAD::INF)
            {
                continue;
            }
            if ( ! _xFeas.isDefined() || f < _xFeas.getF())
            {
                _xFeas = evalPoint;
                updatedFeas = true;
            }
        }
    }

    // Second loop updates the barrier infeasible points.
    for (const auto & evalPoint : evalPointList)
    {

        const NOMAD::Double& h = evalPoint.getH();
        // Consider infeasible points with h < INF
        // Do not consider points for which h > hmax
        if (!h.isDefined() || (h == NOMAD::INF) ||
            ((_hMax < NOMAD::INF) && (h > _hMax)))
        {
            continue;
        }
        else if (h > 0)
        {
            updatedInf = true;
            _xTmpInf.push_back(evalPoint);
        }
    }

    if (updatedInf || updatedFeas)
    {
        _incumbentsAndHMaxUpToDate = false;
    }

    //
    // Ready to update hMax and the incumbents (if requested)
    NOMAD::Double hMaxPrev = _hMax;
    if (!_incumbentsAndHMaxUpToDate)
    {
        NOMAD::SuccessType feasSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
        // Step One: Set xIncFeas (feasible incumbent) from xFeas

        if (_xFeas.isDefined())
        {
            NOMAD::Double f = _xFeas.getF();
            NOMAD::Double fMinRef=_xIncFeas.getF(); // fMinRef not defined if _xIncFeas not defined
            if ( ! _xIncFeas.isDefined() || f < fMinRef)
            {
                _xIncFeas = _xFeas;
                updatedIncFeas = true;
                feasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
            }
        }

        // Step two: Detect the infeasible success type
        NOMAD::SuccessType infeasSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
        if ( !_xIncInf.isDefined() && !_xTmpInf.empty())
        {
            // Case with no previous incumbent (hence, no barrier point).
            infeasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else if (!_xTmpInf.empty())
        {
            // Note: An improving point (with respect to the infeasible incumbent) may have been added at the current iteration  -> loop on _xTmpInf.
            std::vector<NOMAD::SimpleEvalPoint>::iterator it;
            for ( it = _xTmpInf.begin(); it < _xTmpInf.end(); it++ )
            {
                NOMAD::SuccessType successType = computeSuccessType( *it, _xIncInf);
                if (successType > infeasSuccessType)
                {
                    infeasSuccessType = successType;
                }
                if (infeasSuccessType == NOMAD::SuccessType::FULL_SUCCESS)
                {
                    break;
                }
            }
        }
        
        // Step three: remove dominated points
        bool updatedIncInf = removeInfeasibleDominatedPoints();
        
//        // TEMP for testing
//        auto compFunc = [](const NOMAD::SimpleEvalPoint & a, const NOMAD::SimpleEvalPoint & b) { return a.getH() < b.getH(); };
//        std::sort(_xInf.begin(), _xInf.end(), compFunc);
//        std::cout<<"================================="<<std::endl;
//        for(const NOMAD::SimpleEvalPoint & x: _xInf)
//        {
//            std::cout<<"F="<<x.getF()<<" H="<<x.getH()<<std::endl;
//        }

        // Step four: update hMax
        if (infeasSuccessType == NOMAD::SuccessType::FULL_SUCCESS
            || feasSuccessType == NOMAD::SuccessType::FULL_SUCCESS)
        {
            // We have a new incumbent: infeasSuccessType or feasSuccessType is FULL_SUCCESS
            if ( _xIncInf.isDefined())
            {
                // Case with a prior infeasible incumbent: it is used to set hMax.
                _hMax = _xIncInf.getH();

            }
            // Else. Case with no prior infeasible incumbent: we do not update hMax. Typically, hMax = INF.
        }
        else if (infeasSuccessType == NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            _hMax = getWorstHInBarrier();
        }
        else if (infeasSuccessType == NOMAD::SuccessType::UNSUCCESSFUL && _xIncInf.isDefined())
        {
            _hMax = _xIncInf.getH();
        }
        if (_hMax < hMaxPrev)
        {
            updatedHMax = true;
        }

        // Step five: remove barrier points with h > hmax and update infeasible incumbent;
        updatedIncInf = setInfeasibleIncumbents();
        
        // Reduction of the number of barrier points using slices between hmin and hmax of infeasible points
        // reduceBarrierInfeasible();

        _incumbentsAndHMaxUpToDate = true; // For the next update
    }

    updated = updatedFeas || updatedInf;
    updatedInc = updatedIncFeas || updatedIncInf;

    return updatedInc;
}



// Remove dominated points from barrier. Keep only one of the points if two points have the same f and h.
// xIncFeas and xIncInf are singletons
bool NOMAD::SimpleProgressiveBarrier::removeInfeasibleDominatedPoints()
{
    bool updatedXInf = false;
    
    auto it1= _xTmpInf.begin();
    std::vector<NOMAD::SimpleEvalPoint>::iterator it2 ;
    for ( it1 = _xTmpInf.begin(); it1 < _xTmpInf.end(); it1++)
    {
        // Find all infeasible non-dominated points
        bool keepXTmp = true;
        for ( it2 = _xInf.begin(); it2 < _xInf.end(); )
        {
            
            // Case 2 dominates 1 --> 1 is dominated. Do not consider.
            // Test and remove doublons
            if ( dominates(*it2, *it1) )
            {
                keepXTmp = false;
                // std::cout << " ---> dominated" ;
                break;
            }
            // Case 1 dominates 2 --> 2 must be removed.
            if (dominates( *it1 , *it2))
            {
                it2 = _xInf.erase(it2);

            }
            else
            {
                if( (it2)->getF() == (it1)->getF() && (it2)->getH() == (it1)->getH())
                {
                    it2 = _xInf.erase(it2);
                }
                else
                {
                    it2++;
                }
            }
        }
        if (keepXTmp)
        {
            _xInf.push_back(*it1);
            updatedXInf = true;
        }
        // std::cout << std::endl;
    }
    
    

    // Done with xTmpInf.
    _xTmpInf.clear();
 
    return updatedXInf;
}

// Keep one infeasible point with smallest f and one infeasible point with smallest h.
bool NOMAD::SimpleProgressiveBarrier::updateBarrierPointsKeep2()
{
    // Case where no points are added. _xInf[0] is removed
    if (_xTmpInf.empty())
    {
        // Remove _xInf[0]
        if (_xInf.size() == 2)
        {
            _xInf.erase(_xInf.begin());
        }
        return true;
    }
    
    // The simple barrier 2 must always have at least one point and two points at most (P0 and P1).
    // Point 0 to replace _xInf[0], the best F, the current incumbent, must always exist
    // Point 1 to replace _xInf[1], the best H, may not exist.
    auto itForP0 = _xTmpInf.end(), itForP1 = _xTmpInf.end();
    
    NOMAD::Double FForP0=NOMAD::INF, HForP1=NOMAD::INF, tmpF, tmpH;
    if (_xInf.size() == 2)
    {
        // Add _xInf[1]. It could be Point 1 again
        _xTmpInf.push_back(_xInf[1]);
    }
    
    for (auto it = _xTmpInf.begin(); it< _xTmpInf.end(); it++)
    {
        // Search new infeasible point with the smallest f
        tmpF = it->getF();
        if (tmpF < FForP0 )
        {
            itForP0 = it;
            FForP0 = tmpF;
        }
        
        // Search new infeasible point with smallest h
        tmpH = it->getH();
        if ( tmpH < HForP1 )
        {
            itForP1 = it;
            HForP1 = tmpH;
        }
    }
    if ( itForP0 == _xTmpInf.end())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,"We must have a P0");
    }
        
    // Update _xInf[0] with P0
    if (!_xInf.empty())
    {
        _xInf[0]=*itForP0;
    }
    else
    {
        _xInf.push_back(*itForP0);
    }

    // If available,update _xinf[1] with P1
    if ( itForP1 != _xTmpInf.end())
    {
        if (_xInf.size() > 1)
        {
            _xInf[1]=*itForP1;
        }
        else
        {
            _xInf.push_back(*itForP1);
        }
    }

    // Done with xTmpInf.
    _xTmpInf.clear();
 
    return true;
}



// Remove the barrier points with h > hmax;
// Set the infeasible incumbents from the barrier:
// keeps only the undominated infeasible point (no duplicates) satisfying h <= hmax with the lowest f.
bool NOMAD::SimpleProgressiveBarrier::setInfeasibleIncumbents()
{
    if (_xInf.empty())
    {
        return false;
    }
    

    NOMAD::Double prevXIncInfH = NOMAD::INF;
    NOMAD::SimpleEvalPoint prevXIncInf;
    if ( _xIncInf.isDefined() )
    {
        prevXIncInfH = _xIncInf.getH(); // NOT necessarily hMax.
        prevXIncInf = _xIncInf;
    }
    
    if (_xInf.size() == 1)
    {
        _xIncInf = _xInf[0];
        if (prevXIncInfH == _xInf[0].getH())
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    
    bool updatedInfInc = false;
    NOMAD::Double maxH = 0.0;
    auto it1= _xInf.begin();
    NOMAD::SimpleEvalPoint ptMaxH = *it1; // Base value for case with a single _xInf.
    it1++;
    while ( it1 !=_xInf.end() )
    {
        // std::cout << "f=" << (*it1_fs)[0] << " h="<< (*it1_h) ;
        // Remove barrier points with h>hMax
        NOMAD::Double h = (it1)->getH();
        if (h > _hMax)
        {
            it1 = _xInf.erase(it1);
            // std::cout << " ---> erased (h>hMax)" <<std::endl;
            continue;
        }
        
        if (h > maxH)
        {
            ptMaxH = *it1;
            maxH = h;
        }
        
        it1++;
    }

    // std::cout << "--------------------------------" << std::endl;
    
    if (_xInf.empty() || !ptMaxH.isDefined())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"All barrier points have been removed.");
    }
    
    _xIncInf= ptMaxH;
    
    // Infeasible incumbent(s) have been updated.
    if (prevXIncInf != _xIncInf )
    {
        updatedInfInc = true;
    }
    
    return updatedInfInc;
}

void NOMAD::SimpleProgressiveBarrier::reduceBarrierInfeasible()
{
    if (_xInf.size() < 5)
    {
        return;
    }
    
    std::vector<NOMAD::SimpleEvalPoint>::iterator it1;
    
    
    auto compFunc = [](const NOMAD::SimpleEvalPoint & a, const NOMAD::SimpleEvalPoint & b) { return a.getH() < b.getH(); };
    
    std::sort(_xInf.begin(), _xInf.end(), compFunc);
    
    
    size_t nbDH = 2;
    NOMAD::Double Hbase = _xInf[0].getH();
    NOMAD::Double DH = (_hMax - Hbase)/ (double)nbDH;
    
    
    it1 = _xInf.begin();
    // Slice DH in nbDH
    for (size_t i= 1 ; i < nbDH+1 ; i++)
    {
        NOMAD::Double h1 = Hbase + (double)i*DH;
        NOMAD::Double h2 = Hbase + (double)(i+1)*DH;
        
        // Remove points in slice, except the first one
        it1++;
        while (it1 < _xInf.end())
        {
           if (it1->getH() < h1)
           {
               it1 = _xInf.erase(it1);
           }
           else
           {
               break;
           }
        }
        // Find the next slice i that contains the current element
        while (it1 < _xInf.end() && it1->getH() > h2 && h2 < _hMax)
        {
            i++;
            h2 = Hbase + (double)(i+1)*DH;
        }
    }
}


NOMAD::Double NOMAD::SimpleProgressiveBarrier::getWorstHInBarrier() const
{
    if (_xInf.empty())
    {
        return _hMax;
    }
    // Loop on all points to find the one with max h and h < hMax (previous hMax)
    NOMAD::Double maxH = 0.0;
    NOMAD::Double hIncInf = _xIncInf.getH();
    NOMAD::Double H;
    for (const auto & xInf : _xInf)
    {
        // Detect max h infeasible barrier (non dominated) point such that h<hMax (the point with h just below hmax).
        H =  xInf.getH();
        if (H > maxH && H < hIncInf)
        {
            maxH = xInf.getH();
        }
    }
    return maxH;
}


bool NOMAD::SimpleProgressiveBarrier::dominates(const NOMAD::SimpleEvalPoint & p1, const NOMAD::SimpleEvalPoint & p2) const
{



    // Comparing objective vectors of different size is undefined
    if ( !p1.isDefined() || !p2.isDefined())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"F1, H1, F2 and H2 must be defined.");
    }
    
    if (p1.getH() <=0 && p2.getH() <=0) // Both feasible
    {
        // 1 and 2 are both feasible, and 1 dominates 2.
        if (p1.getF() < p2.getF())
        {
            return true;
        }
    }
    else if (p1.getH() > 0 && p2.getH() > 0) // Both infeasible
    {
        if (p1.getH() != NOMAD::INF && p1.getH() <= _hMax)
        {
            if ((p1.getF() <= p2.getF()) && (p1.getH() <= p2.getH()) && ((p1.getF() < p2.getF()) || (p1.getH() < p2.getH())))
            {
                return true;
            }
        }
    }
    return false;
}


NOMAD::SuccessType NOMAD::SimpleProgressiveBarrier::computeSuccessType(const NOMAD::SimpleEvalPoint & p1, const NOMAD::SimpleEvalPoint & p2) const
{
    // UNSUCCESSFUL,       // Failure
    // PARTIAL_SUCCESS,    // Partial success (improving). Found an infeasible
    //                        solution with a better h. f is worse.
    // FULL_SUCCESS        // Full success (dominating)
    NOMAD::SuccessType success = NOMAD::SuccessType::UNDEFINED; /// Will trigger exception if not changed
    
    if (! p1.isDefined())
    {
        return NOMAD::SuccessType::UNSUCCESSFUL;
    }
    
    if (! p2.isDefined())
    {
        if (p1.getH() > _hMax || p1.getH() == NOMAD::INF)
        {
            // Even if f2 is not defined, this case is not successful.
            success = NOMAD::SuccessType::UNSUCCESSFUL;
        }
        else
        {
            // A new infeasible point, without prior infeasible point, is partial success,
            // not a full success.
            if ( p1.getH() <= 0)
            {
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
            else
            {
                success = NOMAD::SuccessType::PARTIAL_SUCCESS;
            }
        }
    }
    else
    {

        
        success = NOMAD::SuccessType::UNSUCCESSFUL;
        if (p1.getH() <=0 && p2.getH() <=0) // Both feasible
        {
            // 1 and 2 are both feasible, and 1 dominates 2.
            if (p1.getF() < p2.getF())
            {
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
        }
        else if (p1.getH() > 0 && p2.getH() > 0) // Both infeasible
        {
            if (p1.getH() != NOMAD::INF && p1.getH() <= _hMax)
            {
                if ((p1.getF() <= p2.getF()) && (p1.getH() <= p2.getH()) && ((p1.getF() < p2.getF()) || (p1.getH() < p2.getH())))
                {
                    success = NOMAD::SuccessType::FULL_SUCCESS;
                }
                else if (p1.getH() < p2.getH() && p1.getF() > p2.getF())
                {
                    // Partial success (improving). Found an infeasible
                    // solution with a better h. f is worse.
                    success = NOMAD::SuccessType::PARTIAL_SUCCESS;
                }
            }
        }
        // Comparing feasible and infeasible point => unsuccessful.
        else
        {
            success = NOMAD::SuccessType::UNSUCCESSFUL;
        }
    }
    
    return success;
}
