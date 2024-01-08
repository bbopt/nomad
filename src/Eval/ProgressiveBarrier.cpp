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

#include "../Cache/CacheBase.hpp"
#include "../Eval/ComputeSuccessType.hpp"
#include "../Eval/ProgressiveBarrier.hpp"
#include "../Output/OutputQueue.hpp"

void NOMAD::ProgressiveBarrier::init(const NOMAD::Point& fixedVariable,
                          NOMAD::EvalType  evalType,
                          NOMAD::ComputeType computeType,
                          bool barrierInitializedFromCache)
{
    // NOTE: the _xFeas[0] and _xInf[0] are not necessarily the incumbents
    // after this init.
    // The second init function calls updateWithPoints and makes sure the incumbents
    // and hMax are set.


    if (fixedVariable.isEmpty())
    {
        std::string s = "Error: Fixed variable of dimension 0";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    if (barrierInitializedFromCache)
    {
        checkCache();

        std::vector<NOMAD::EvalPoint> cachePoints;

        // Get best feasible and infeasible solutions from cache.
        // Points from cache are in full dimension. Convert them
        // to subproblem dimension.
        // Update the barrier with points.
        // Important: do not update infeasible incumbent and hmax. It is is done in the second init. We don't want to update infeasible incumbent and hmax two times!
        if (NOMAD::CacheBase::getInstance()->findBestFeas(cachePoints, fixedVariable, evalType, computeType) > 0)
        {
            for (const auto &evalPoint : cachePoints)
            {
                NOMAD::EvalPointPtr evalPointSub = std::make_shared<NOMAD::EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariable));
                _xFeas.push_back(evalPointSub);
            }
            _incumbentsAndHMaxUpToDate = false;
        }
        if (NOMAD::CacheBase::getInstance()->findFilterInf(cachePoints, _hMax, fixedVariable, evalType, computeType) > 0)
        {
            for (const auto &evalPoint : cachePoints)
            {
                // Points in progressive barrier must have h < INF.
                if (evalPoint.getH(evalType,computeType) < NOMAD::INF)
                {
                    NOMAD::EvalPointPtr evalPointSub = std::make_shared<NOMAD::EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariable));
                    _xInf.push_back(evalPointSub);
                }
            }
            _incumbentsAndHMaxUpToDate = false;
        }
    }
}

void NOMAD::ProgressiveBarrier::init(const NOMAD::Point& fixedVariable,
                          NOMAD::EvalType  evalType,
                          const std::vector<NOMAD::EvalPoint>& evalPointList,
                          NOMAD::ComputeType computeType)
{

    // Constructor's call to update should not update ref best points.
    updateWithPoints(evalPointList, evalType, computeType, true, true /*true update infeasible incumbent and hmax*/ );

    // Check: xIncFeas or xIncInf could be non-evaluated, but not both.
    auto xIncFeas = getCurrentIncumbentFeas();
    auto xIncInf = getCurrentIncumbentInf();
    if (   (nullptr == xIncFeas || nullptr == xIncFeas->getEval(evalType))
        && (nullptr == xIncInf  || nullptr == xIncInf->getEval(evalType)))
    {
        std::string s = "Barrier constructor: no xIncFeas and xIncInf  properly defined. This may cause problems. \n";
        if (nullptr != xIncFeas)
        {
            s += "There are " + std::to_string(_xIncFeas.size()) + " feasible incumbents, the first one is:\n";
            s += xIncFeas->displayAll();
        }
        if (nullptr != xIncInf)
        {
            s += "There are " + std::to_string(_xInf.size()) + " infeasible incumbents, the first one is:\n";
            s += xIncInf->displayAll();
        }
    }

    checkHMax();
}



void NOMAD::ProgressiveBarrier::setHMax(const NOMAD::Double &hMax)
{
    _hMax = hMax;

    checkHMax();
}



void NOMAD::ProgressiveBarrier::updateRefBests()
{
    _refBestFeas = getCurrentIncumbentFeas();;
    _refBestInf  = getCurrentIncumbentInf();;
}


NOMAD::SuccessType NOMAD::ProgressiveBarrier::getSuccessTypeOfPoints(const EvalPointPtr xFeas,
                                                          const EvalPointPtr xInf,
                                                          NOMAD::EvalType evalType,
                                                          NOMAD::ComputeType computeType)
{
    NOMAD::SuccessType successType = SuccessType::UNSUCCESSFUL;
    NOMAD::SuccessType successType2 = SuccessType::UNSUCCESSFUL;

    NOMAD::EvalPointPtr newBestFeas,newBestInf;

    // Get the reference best points (should work for opportunistic or not)
    auto refBestFeas = getCurrentIncumbentFeas();
    auto refBestInf = getCurrentIncumbentInf();

    // Set the iter success based on best feasible and best infeasible points found compared to initial point.
    if (nullptr != refBestFeas || nullptr != refBestInf)
    {
        // Compute success
        // Get which of newBestFeas and newBestInf is improving
        // the solution. Check newBestFeas first.
        NOMAD::ComputeSuccessType computeSuccess(evalType, computeType);

        if (nullptr != refBestFeas)
        {
            successType = computeSuccess(xFeas, refBestFeas,_hMax);
        }
        if (nullptr != refBestInf)
        {
            successType2 = computeSuccess(xInf, refBestInf,_hMax);
        }
        if (successType2 > successType)
        {
            successType = successType2;
        }
    }
    return successType;
}


// Points from evalPointList are already in subproblem dimension.
bool NOMAD::ProgressiveBarrier::updateWithPoints(const std::vector<NOMAD::EvalPoint>& evalPointList,
                                      NOMAD::EvalType evalType,
                                      NOMAD::ComputeType computeType,
                                      const bool keepAllPoints /* Not used here */,
                                      const bool updateIncumbentsAndHmax )
{

    bool updated = false;
    bool updatedFeas = false;
    bool updatedInf = false;
    bool updatedInc = false;
    bool updatedIncFeas = false;
    bool updatedIncInf = false;
    bool updatedHMax = false;


    std::string s;  // for output info

    OUTPUT_DEBUG_START
    s = "Updating barrier with " + std::to_string(evalPointList.size()) + " suggested points";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    s = "Current barrier: ";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    std::vector<std::string> vs = display(4);
    for (const auto & s: vs)
    {
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    OUTPUT_DEBUG_END

    // Loop to add the feasible points into the barrier.
    // The feasible incumbents are identified later.
    for (const auto & evalPoint : evalPointList)
    {

        auto eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            OUTPUT_DEBUG_START
            if (nullptr == eval)
            {
                s = "Eval is NULL, continue";
            }
            else if (NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
            {
                s = "Eval status is " + NOMAD::enumStr(eval->getEvalStatus());
                s += ", continue";
            }
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            NOMAD::OutputQueue::Flush();
            OUTPUT_DEBUG_END

            // Suggested point is not good.
            continue;
        }


        if (eval->isFeasible(computeType))
        {
            OUTPUT_DEBUG_START
            s = "Point suggested to update barrier (feasible): " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END

            // If eval point is already in barrier, do not add.
            // Tag is the unique identifier of EvalPoint
            size_t tag = evalPoint.getTag();
            std::vector<NOMAD::EvalPointPtr>::iterator it = std::find_if(_xFeas.begin(), _xFeas.end(),[tag](const NOMAD::EvalPointPtr& pt){ return pt->getTag() == tag ; });
            if (it !=_xFeas.end())
            {
                OUTPUT_DEBUG_START
                s = "Point already in barrier.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
                continue;
            }

            OUTPUT_DEBUG_START
            s = "Point inserted in feasible barrier: " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            _xFeas.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
            updatedFeas = true;

        }
    }

    // Second loop updates the barrier infeasible points.
    for (const auto & evalPoint : evalPointList)
    {
        auto eval = evalPoint.getEval(evalType);

        // Exclude points that are not eval ok
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus() )

        {
            // Suggested point is not good.
            continue;
        }
        // Consider infeasible points with h < INF
        if (!eval->isFeasible(computeType) && eval->getH(computeType) < NOMAD::INF)
        {
            OUTPUT_DEBUG_START
            s = "Point suggested to update barrier (infeasible): " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END

            // If eval point is already in barrier, do not add
            // Tag is the unique identifier of EvalPoint
            size_t tag = evalPoint.getTag();
            std::vector<NOMAD::EvalPointPtr>::iterator it = std::find_if(_xInf.begin(), _xInf.end(),[tag](const NOMAD::EvalPointPtr& pt){ return pt->getTag() == tag ; });
            if (it !=_xInf.end())
            {
                OUTPUT_DEBUG_START
                s = "Point already in barrier.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
                continue;
            }

            // Do not consider points for which h > hmax
            NOMAD::Double h = eval->getH(computeType);
            if (!h.isDefined() || (h == NOMAD::INF) ||
                ((_hMax < NOMAD::INF) && (h > _hMax)))
            {
                OUTPUT_DEBUG_START
                if (h.isDefined())
                {
                    s = "H is too large: ";
                    s += h.display(NOMAD::DISPLAY_PRECISION_FULL) + " > " + _hMax.display(NOMAD::DISPLAY_PRECISION_FULL) + ", continue.";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                }
                else
                {
                    s = "H is undefined";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                }
                OUTPUT_DEBUG_END
                continue;
            }

            updatedInf = true;
            _xInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
            OUTPUT_DEBUG_START
            s = "Add in barrier xInf: " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
        }
    }

    if (updatedInf || updatedFeas)
    {
        _incumbentsAndHMaxUpToDate = false;
    }

    //
    // Ready to update hMax and the incumbents (if requested)
    NOMAD::Double hMaxPrev = _hMax;
    if (updateIncumbentsAndHmax && !_incumbentsAndHMaxUpToDate)
    {
        NOMAD::SuccessType feasSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
        // Step One: Select the best xFeas (feasible incumbent) and remove points with f>f(xFeas)
        // After that step xFeas and xIncFeas contain the same points
        if (_xFeas.size() > 0)
        {
            // Detect fMin
            NOMAD::Double fMin = NOMAD::INF;
            NOMAD::Double fMinRef=fMin;
            if (_xIncFeas.size() > 0)
            {
                fMinRef = _xIncFeas[0]->getEval(evalType)->getF(computeType);
            }
            std::vector<NOMAD::EvalPointPtr>::iterator it;
            for ( it = _xFeas.begin(); it < _xFeas.end(); it++ )
            {
                fMin = min(fMin, (*it)->getEval(evalType)->getF(computeType));
            }

            // Remove all points above fMin. There could be more than one point remaining.
            _xIncFeas.clear();
            for ( it = _xFeas.begin(); it < _xFeas.end(); )
            {
                if ((*it)->getEval(evalType)->getF(computeType) > fMin)
                {
                    it = _xFeas.erase(it);
                }
                else
                {
                    ++it;
                }
            }
            // Copy points to _xIncFeas
            copy(_xFeas.begin(), _xFeas.end(), back_inserter(_xIncFeas));

            // Test for FULL_SUCCESS
            if (fMin<fMinRef)
            {
                updatedIncFeas = true;
                feasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
                OUTPUT_DEBUG_START
                s = "New dominating xFeas: " + _xFeas[0]->display();
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
            }
        }

        // Step two: Detect the infeasible success type
        NOMAD::SuccessType infeasSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
        if (_xIncInf.size() == 0 && _xInf.size() > 0)
        {
            // Case with no previous incumbent (hence, no barrier point).
            infeasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
        }
        else if (_xInf.size() > 0)
        {
            // Note: An improving point (with respect to the infeasible incumbent) may have been added at the current iteration or was present from a previous iteration -> loop on _xInf
            std::vector<NOMAD::EvalPointPtr>::iterator it;
            for ( it = _xInf.begin()+1; it < _xInf.end(); it++ )
            {
                auto eval = (*it)->getEval(evalType);
                NOMAD::SuccessType successType = NOMAD::Eval::computeSuccessType(eval, _xIncInf[0]->getEval(evalType), computeType);
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

        // Step three: update hMax
        if (infeasSuccessType == NOMAD::SuccessType::FULL_SUCCESS
            || feasSuccessType == NOMAD::SuccessType::FULL_SUCCESS)
        {
            // We have a new incumbent: infeasSuccessType or feasSuccessType is FULL_SUCCESS
            if (_xIncInf.size() != 0)
            {
                // Case with a prior infeasible incumbent: it is used to set hMax.
                _hMax = _xIncInf[0]->getEval(evalType)->getH(computeType);

            }
            // Else. Case with no prior infeasible incumbent: we do not update hMax. Typically hMax = INF.
        }
        else if (infeasSuccessType == NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            _hMax = getWorstHInBarrier(evalType, computeType);
        }
        else if (infeasSuccessType == NOMAD::SuccessType::UNSUCCESSFUL && _xIncInf.size() > 0)
        {
            _hMax = _xIncInf[0]->getEval(evalType)->getH(computeType);
        }
        if (_hMax < hMaxPrev)
        {
            updatedHMax = true;
        }

        // Step four:
        // - remove points from xInf that have h>hMax.
        // - transfer points with h<hMax into xIncInf
        updatedIncInf = setInfeasibleIncumbents(evalType, computeType);

        _incumbentsAndHMaxUpToDate = true;
    }

    updated = updatedFeas || updatedInf;
    updatedInc = updatedIncFeas || updatedIncInf;

    // Set n and check that all points in barrier have the same dimension
    if (updated || updatedInc)
    {
        setN();
    }

    OUTPUT_DEBUG_START
    if (updatedHMax)
    {
        s = "Updated barrier hMax: from " + hMaxPrev.display() + " to " + _hMax.display() ;
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    if (updatedInc)
    {
        s = "Update incumbent(s) in barrier. ";
        if (updatedIncFeas)
        {
            s += "New feasible incumbent: " + _xIncFeas[0]->display();
        }
        else if (updatedIncInf && _xInf.size() > 0)
        {
            s += "New infeasible incumbent: " + _xIncInf[0]->display();
        }
        else if (updatedIncInf && _xInf.size() == 0)
        {
            s += "Previous infeasible incumbent has been removed. No more infeasible incumbent";
        }
    }
    else if(updated)
    {
        s = "Barrier incumbents not updated but points added into the barrier";
    }
    else
    {
        s = "No points added to the barrier";

    }
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    return updatedInc;
}



// Nice formatted display.
std::vector<std::string> NOMAD::ProgressiveBarrier::display(const size_t max) const
{
    std::vector<std::string> vs;

    auto allXFeas = getAllXFeas();
    auto allXInf  = getAllXInf();
    size_t nbXFeas = 0;
    size_t nbXInf = 0;

    for (auto xFeas : allXFeas)
    {
        vs.push_back("X_FEAS " + xFeas->displayAll());
        nbXFeas++;
        if (nbXFeas >= max && allXFeas.size() > max)
        {
            vs.push_back("... (total " + std::to_string(allXFeas.size()) + ")");
            break;
        }
    }
    for (auto xInf : allXInf)
    {
        vs.push_back("X_INF " + xInf->displayAll());
        nbXInf++;

        if (nbXInf >= max && allXInf.size() > max)
        {
            vs.push_back("... (total " + std::to_string(allXInf.size()) + ")");
            break;
        }
    }
    vs.push_back("H_MAX " + getHMax().display(NOMAD::DISPLAY_PRECISION_FULL));
    vs.push_back("Ref Best Feasible:   " + (_refBestFeas ? _refBestFeas->displayAll() : "NULL"));
    vs.push_back("Ref Best Infeasible: " + (_refBestInf ? _refBestInf->displayAll() : "NULL"));

    return vs;
}

// Set the infeasible incumbents.
// This function also remove points with h > hMax.
// The infeasible incumbents are the undominated infeasible points satisfying h < hMax with the lowest f.
bool NOMAD::ProgressiveBarrier::setInfeasibleIncumbents(NOMAD::EvalType evalType, NOMAD::ComputeType computeType)
{
    if (_xInf.size() == 0)
    {
        return false;
    }


    std::vector<NOMAD::EvalPointPtr>::iterator it1= _xInf.begin();

    std::vector<NOMAD::ArrayOfDouble> f_xInf, f_xIncInf;
    std::vector<NOMAD::Double> h_xInf, h_xIncInf;
    // Vectors for fs and hs
    while ( it1 !=_xInf.end() )
    {
        f_xInf.push_back((*it1)->getEval(evalType)->getFs(computeType));
        h_xInf.push_back((*it1)->getEval(evalType)->getH(computeType));
        it1++;
    }


    NOMAD::Double prevXIncInfH = NOMAD::INF;
    NOMAD::EvalPointPtr prevXIncInf = nullptr;
    if (_xIncInf.size() > 0 )
    {
        prevXIncInfH = _xIncInf[0]->getEval(evalType)->getH(computeType); // Not necessarily hMax!
        prevXIncInf = _xIncInf[0];
    }

    _xIncInf.clear();

    bool updatedInfInc = false;
    NOMAD::Double maxH = 0.0;


    std::vector<NOMAD::EvalPointPtr>::iterator it2;
    std::vector<NOMAD::ArrayOfDouble>::iterator it1_fs = f_xInf.begin() ;
    std::vector<NOMAD::Double>::iterator it1_h = h_xInf.begin();

    it1 = _xInf.begin();

    while ( it1 !=_xInf.end() )
    {
        // Remove points with h>hMax
        if ((*it1_h) > _hMax)
        {
            it1 = _xInf.erase(it1);
            it1_h = h_xInf.erase(it1_h);
            it1_fs = f_xInf.erase(it1_fs);
            continue;
        }

        // Find all infeasible non-dominated points
        bool isDominated = false;
        std::vector<NOMAD::ArrayOfDouble>::iterator it2_fs = f_xInf.begin() ;
        std::vector<NOMAD::Double>::iterator it2_h = h_xInf.begin();
        for ( it2 = _xInf.begin(); it2 < _xInf.end(); it2++, it2_fs++, it2_h++)
        {
            if (it1 == it2)
            {
                continue;
            }
            // Case ep2 dominates ep1 --> ep1 is dominated and cannot be the infeasible incumbent.
            if ( dominates(*it2_fs,*it2_h, *it1_fs, *it1_h))
            {
                isDominated = true;
                break;
            }
        }

        // No point dominates ep1 --> it is a POTENTIAL infeasible incumbent
        // Also determine what is the max h value of all non-dominated points
        if (! isDominated)
        {
            _xIncInf.push_back(*it1);
            f_xIncInf.push_back(*it1_fs);
            h_xIncInf.push_back(*it1_h);
            NOMAD::Double H = *it1_h;
            maxH = max(maxH,H);
        }
        it1++;
        it1_h++;
        it1_fs++;
        // We may find other non dominated points. We only want those with h = max H (closer to hMax).
    }

    // Find infeasible incumbents among non-dominated points (=those with lowest f-value aka highest h value)
    // Do another loop to remove points that do not have h=maxH

    if (_xIncInf.size() > 1)
    {
        it1 = _xIncInf.begin();
        it1_fs = f_xIncInf.begin();
        it1_h = h_xIncInf.begin();
        while( it1 != _xIncInf.end())
        {
            if (*it1_h < maxH)
            {
                it1 = _xIncInf.erase(it1);
                it1_fs = f_xIncInf.erase(it1_fs);
                it1_h = h_xIncInf.erase(it1_h);
            }
            else
            {
                it1++;
                it1_h++;
                it1_fs++;
            }
        }
        if (_xIncInf.size() == 0)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"All infeasible incumbents have been removed.");
        }
    }

    // H of infeasible incumbent(s) has changed. First infeasible incumbents has changed. Infeasible incumbent(s) have been updated.
    if (maxH < prevXIncInfH || prevXIncInf != _xIncInf[0])
    {
        updatedInfInc = true;
    }

    return updatedInfInc;
}

NOMAD::Double NOMAD::ProgressiveBarrier::getWorstHInBarrier(NOMAD::EvalType evalType, NOMAD::ComputeType computeType) const
{
    if (_xInf.size() == 0)
    {
        return _hMax;
    }
    // Loop on all points to find the one with max h and h < hMax (previous hMax)
    NOMAD::Double maxH = 0.0;
    NOMAD::Double hIncInf = _xIncInf[0]->getEval(evalType)->getH(computeType);
    NOMAD::Double H;
    for (size_t i = 0; i < _xInf.size() ; i++ )
    {
        // Detect max h infeasible barrier point such that h<hMax (the point with h just below hmax).
        H =  _xInf[i]->getEval(evalType)->getH(computeType);
        if (H > maxH && H < hIncInf)
        {
            maxH = _xInf[i]->getEval(evalType)->getH(computeType);
        }
    }
    return maxH;
}

NOMAD::EvalPointPtr NOMAD::ProgressiveBarrier::getCurrentIncumbentFeas() const
{
    if (_xIncFeas.size() > 0)
    {
        return _xIncFeas[0];
    }
    else
    {
        return nullptr;
    }

}

NOMAD::EvalPointPtr NOMAD::ProgressiveBarrier::getCurrentIncumbentInf() const
{
    if (_xIncInf.size() > 0)
    {
        return _xIncInf[0];
    }
    else
    {
        return nullptr;
    }

}


bool NOMAD::ProgressiveBarrier::dominates(const NOMAD::ArrayOfDouble & f1, const NOMAD::Double & h1, const NOMAD::ArrayOfDouble & f2, const NOMAD::Double & h2) const
{

    NOMAD::CompareType compareFlag = NOMAD::CompareType::UNDEFINED;

    // Comparing objective vectors of different size is undefined
    if (f1.size() != f2.size())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot compare vectors of different size.");
    }

    // This is the same code than Eval::compMO. It is adapted for given f and h.
    // The comparison code has been adapted from
    // Jaszkiewicz, A., & Lust, T. (2018).
    // ND-tree-based update: a fast algorithm for the dynamic nondominance problem.
    // IEEE Transactions on Evolutionary Computation, 22(5), 778-791.
    if (h1.todouble() < NOMAD::Double::getEpsilon() && h2.todouble() < NOMAD::Double::getEpsilon())
    {
        bool isbetter = false;
        bool isworse = false;
        for (size_t i = 0; i < f1.size(); ++i)
        {
            if (f1[i].todouble() < f2[i].todouble())
            {
                isbetter = true;
            }
            if (f2[i].todouble() < f1[i].todouble())
            {
                isworse = true;
            }
            if (isworse && isbetter)
            {
                break;
            }
        }
        if (isworse)
        {
            compareFlag = isbetter ? NOMAD::CompareType::INDIFFERENT : NOMAD::CompareType::DOMINATED;
        }
        else
        {
            compareFlag = isbetter ? NOMAD::CompareType::DOMINATING : NOMAD::CompareType::EQUAL;
        }
    }
    else if (h1.todouble() >= NOMAD::Double::getEpsilon() && h2.todouble() >= NOMAD::Double::getEpsilon())
    {
        if (h1.todouble() != NOMAD::INF)
        {
            bool isbetter = false;
            bool isworse = false;
            for (size_t i = 0; i < f1.size(); ++i)
            {
                if (f1[i].todouble() < f2[i].todouble())
                {
                    isbetter = true;
                }
                if (f2[i].todouble() < f1[i].todouble())
                {
                    isworse = true;
                }
                if (isworse && isbetter)
                {
                    break;
                }
            }
            if (!(isworse && isbetter))
            {
                if (h1.todouble() < h2.todouble())
                {
                    isbetter = true;
                }
                if (h2.todouble() < h1.todouble())
                {
                    isworse = true;
                }
            }
            if (isworse)
            {
                compareFlag = isbetter ? NOMAD::CompareType::INDIFFERENT : NOMAD::CompareType::DOMINATED;
            }
            else
            {
                compareFlag = isbetter ? NOMAD::CompareType::DOMINATING : NOMAD::CompareType::EQUAL;
            }
        }
    }

    // Comparing an infeasible objective vector with a feasible objective vector is always UNDEFINED.
    return (compareFlag== NOMAD::CompareType::DOMINATING);
}
