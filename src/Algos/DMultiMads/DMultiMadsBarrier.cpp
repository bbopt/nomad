/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
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
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Cache/CacheSet.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Eval/MeshBase.hpp"
#include "../../Output/OutputQueue.hpp"

// Initialization from cache
void NOMAD::DMultiMadsBarrier::init(const NOMAD::Point& fixedVariables,
                            NOMAD::EvalType evalType,
                            NOMAD::ComputeType computeType,
                            bool barrierInitializedFromCache)
{
    std::list<NOMAD::EvalPoint> cachePoints;

    if (fixedVariables.isEmpty())
    {
        std::string s = "Error: Fixed variable of dimension 0";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
    if (barrierInitializedFromCache)
    {
        checkCache();


        // Get best feasible and infeasible solutions from cache.
        // Point from cache are in full dimension.
        // Convert them to subproblem dimension.
        // NB: all solutions (and not only the non dominated ones are considered)
        auto cache = CacheBase::getInstance().get();
        if (cache->findBestFeas(cachePoints, fixedVariables, evalType, computeType) > 0)
        {
            for (const auto & evalPoint : cachePoints)
            {
                NOMAD::EvalPointPtr evalPointSub = std::make_shared<EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariables));
                _xFeas.push_back(evalPointSub);
            }
            cachePoints.clear();
        }
        if (cache->findBestInf(cachePoints, _hMax, fixedVariables, evalType, computeType) > 0)
        {
            for (const auto &evalPoint : cachePoints)
            {
                // Consider points with h < INF. That is, points that are not excluded by extreme barrier constraints.
                if (evalPoint.getH(evalType,computeType) < NOMAD::INF)
                {
                    NOMAD::EvalPointPtr evalPointSub = std::make_shared<EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariables));
                    _xInf.push_back(evalPointSub);
                }
            }
            cachePoints.clear();
        }

        // Get non dominated infeasible points from cache.
        // Points from cache are in full dimension.
        // Convert them to subproblem dimension.
        if (cache->findFilterInf(cachePoints, _hMax, fixedVariables, evalType, computeType) > 0)
        {
            for (const auto & evalPoint : cachePoints)
            {
                NOMAD::EvalPointPtr evalPointSub = std::make_shared<EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariables));
                _xFilterInf.push_back(evalPointSub);
            }
            cachePoints.clear();
        }
    }

    if (_xFeas.size() > 0 || _xInf.size() > 0)
    {
        setN();

        for (const auto& xFeas: _xFeas)
        {
            checkMeshParameters(*xFeas);
        }
        for (const auto& xInf: _xInf)
        {
            checkMeshParameters(*xInf);
        }
        
        // Update the incumbents used by DMultiMads algo as frameCenter
        updateCurrentIncumbents();
        
    }
}


// And from the list of points given
void NOMAD::DMultiMadsBarrier::init(const NOMAD::Point& fixedVariables,
                            NOMAD::EvalType evalType,
                            const std::vector<NOMAD::EvalPoint>& evalPointList,
                            NOMAD::ComputeType computeType)
{

    bool updated = updateWithPoints(evalPointList, evalType, computeType, true); // All points are considered.

    if (! updated)
         return;
    
    // Need to set bb input type
    // Initialize _bbInputsType. By default, it is a vector of real.
    if (_bbInputsType.empty())
    {
        _bbInputsType = BBInputTypeList(_n, NOMAD::BBInputType::CONTINUOUS);
    }
    if (_bbInputsType.size() != _n)
    {
        std::string s = "Error: Inputs dimensions of DMultiMadsBarrier do not match dimensions of provided input types.";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    // Check: xFeas or xInf could be non-evaluated, but not both.
    if (   (_xFeas.empty() || nullptr == _xFeas[0]->getEval(evalType))
        && (_xInf.empty() || nullptr == _xInf[0]->getEval(evalType)))
    {
        std::string s = "Barrier constructor: xFeas or xInf must be in the barrier.\n";
        if (!_xFeas.empty())
        {
            s += "There are " + std::to_string(_xFeas.size()) + " xFeas, the first one is:\n";
            s += _xFeas[0]->displayAll();
        }
        if (!_xInf.empty())
        {
            s += "There are " + std::to_string(_xInf.size()) + " xInf, the first one is:\n";
            s += _xInf[0]->displayAll();
        }
        if (_xFeas.size() == 0 && _xInf.size() == 0)
        {
            s += "There are no xFeas and no xInf defined.";
        }
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    checkHMax();
}


void NOMAD::DMultiMadsBarrier::setHMax(const NOMAD::Double &hMax)
{
    NOMAD::Double oldHMax = _hMax;

    _hMax = hMax;
    checkHMax();

    if (_hMax < oldHMax)
    {
        updateXInfAndFilterInfAfterHMaxSet();
    }
    // Always update the current incumbent infeasible after changing hmax.
    updateCurrentIncumbentInf();
}


void NOMAD::DMultiMadsBarrier::checkMeshParameters(const NOMAD::EvalPoint &x) const
{
    auto mesh = x.getMesh();
    
    // Mesh can be in sub dimension and point can be in full dimension. This can happen for the EvaluatorControl barrier with fixed variables.
    size_t meshSizeCorrection = 0;
    
    if (mesh->getdeltaMeshSize().size()!= x.size())
    {
        meshSizeCorrection = _fixedVariables.nbDefined();
    }
    
    if (mesh->getdeltaMeshSize().size()+meshSizeCorrection != x.size() ||
        mesh->getDeltaFrameSize().size()+meshSizeCorrection != x.size() ||
        mesh->getMeshIndex().size()+meshSizeCorrection != x.size())
    {
        std::string s = "Error: Mesh parameters dimensions are not compatible with EvalPoint dimension.\n";
        s += "EvalPoint dimensions: " + std::to_string(x.size()) + "\n";
        s += "MeshSize dimensions: " + std::to_string(mesh->getdeltaMeshSize().size()) + "\n";
        s += "FrameSize dimensions: " + std::to_string(mesh->getDeltaFrameSize().size()) + "\n";
        s += "MeshIndex dimensions: " + std::to_string(mesh->getMeshIndex().size());
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    if (!mesh->getdeltaMeshSize().isDefined())
    {
        std::string s = "Error: some MeshSize components of EvalPoint passed to MO Barrier ";
        s += "are not defined.";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
    if (!mesh->getDeltaFrameSize().isDefined())
    {
        std::string s = "Error: some FrameSize components of EvalPoint passed to MO Barrier ";
        s += "are not defined.";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
    if (!mesh->getMeshIndex().isDefined())
    {
        std::string s = "Error: some MeshIndex components of EvalPoint passed to MO Barrier ";
        s += "are not defined.";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
}



void NOMAD::DMultiMadsBarrier::clearXFeas()
{
    _xFeas.clear();
    
    // Update the current incumbents. Both the feasible and  infeasible ones.
    updateCurrentIncumbents();
    
}

void NOMAD::DMultiMadsBarrier::checkXFeas(const NOMAD::EvalPoint &xFeas,
                                  NOMAD::EvalType evalType,
                                  NOMAD::ComputeType computeType)
{
    // If evalType is UNDEFINED, skip this check.
    if (NOMAD::EvalType::UNDEFINED != evalType) {
        if (nullptr == xFeas.getEval(evalType)) {
            throw NOMAD::Exception(
                __FILE__, __LINE__,
                "DMultiMadsBarrier: xFeas must be evaluated before being set.");
        }
        checkXFeasIsFeas(xFeas, evalType, computeType);
    }
}


void NOMAD::DMultiMadsBarrier::checkXFeasIsFeas(const NOMAD::EvalPoint &xFeas,
                                        NOMAD::EvalType  evalType,
                                        NOMAD::ComputeType computeType)
{
    // If evalType is UNDEFINED, skip this check.
    if (NOMAD::EvalType::UNDEFINED != evalType)
    {
        auto eval = xFeas.getEval(evalType);
        if (nullptr != eval && NOMAD::EvalStatusType::EVAL_OK == eval->getEvalStatus())
        {
            NOMAD::Double h = eval->getH(computeType);
            if (!h.isDefined() || 0.0 != h)
            {
                std::string err = "Error: DMultiMadsBarrier: xFeas' h value must be 0.0, got: " + h.display();
                throw NOMAD::Exception(__FILE__,__LINE__,err);
            }
            if (computeType == NOMAD::ComputeType::STANDARD && eval->getFs(computeType).size() != _nobj)
            {
                std::string err = "Error: DMultiMadsBarrier: xFeas' F must be of size " + std::to_string(_nobj);
                err += ", got: F.size() = " + std::to_string(eval->getFs(computeType).size());
                err += " with following F values " + eval->getFs(computeType).display();
                throw NOMAD::Exception(__FILE__,__LINE__,err);
            }
        }
    }
}


NOMAD::EvalPointPtr NOMAD::DMultiMadsBarrier::getFirstXInfNoXFeas() const
{
    NOMAD::EvalPointPtr xInf = nullptr;
    if (_xFilterInf.size() == 0)
    {
        return xInf;
    }

    // Select candidates
    NOMAD::Double minFrameSizeInfElts = getMeshMaxFrameSize(getXInfMinH());
    std::vector<bool> canBeFrameCenter(_xInf.size(), false);
    size_t nbSelectedCandidates = 0;
    for (size_t i = 0; i < _xInf.size(); ++i)
    {
        NOMAD::Double maxFrameSizeElt = getMeshMaxFrameSize(_xInf[i]);
        if (minFrameSizeInfElts <= maxFrameSizeElt)
        {
            canBeFrameCenter[i] = true;
            nbSelectedCandidates += 1;
        }
    }

    // The selection must always work
    if (nbSelectedCandidates == 0)
    {
        xInf = _xInf[0];
    }
    else if (nbSelectedCandidates == 1)
    {
        auto it = std::find(canBeFrameCenter.begin(), canBeFrameCenter.end(), true);
        if (it == canBeFrameCenter.end())
        {
            std::string s = "Error: DMultiMadsBarrier, should not reach this condition";
            throw NOMAD::Exception(__FILE__,__LINE__,s);
        }
        else
        {
            size_t selectedInd = std::distance(canBeFrameCenter.begin(), it);
            xInf = _xInf[selectedInd];
        }
    }
    else if ((nbSelectedCandidates == 2) && (_xInf.size() == 2))
    {
        const NOMAD::Eval* eval1 = _xInf[0]->getEval(NOMAD::EvalType::BB);
        const NOMAD::Eval* eval2 = _xInf[1]->getEval(NOMAD::EvalType::BB);
        auto objv1 = eval1->getFs();
        auto objv2 = eval2->getFs();
        if (objv1.abs().max() > objv2.abs().max())
        {
            xInf = _xInf[0];
        }
        else
        {
            xInf = _xInf[1];
        }
    }
    // More than two points in the barrier.
    else
    {
        // First case: biobjective optimization. Points are already ranked by lexicographic order.
        if (_nobj == 2)
        {
            size_t currentBestInd = 0;
            NOMAD::Double maxGap = -1.0;
            NOMAD::Double currentGap;
            for (size_t obj = 0; obj < _nobj; ++obj)
            {
                // Get extreme values value according to one objective
                NOMAD::Double fmin = _xInf[0]->getEval(NOMAD::EvalType::BB)->getFs()[obj];
                NOMAD::Double fmax = _xInf[_xInf.size()-1]->getEval(NOMAD::EvalType::BB)->getFs()[obj];

                // In this case, it means all elements of _xInf are equal (return the first one)
                if (fmin == fmax)
                {
                    break;
                }

                // Intermediate points
                for (size_t i = 1; i < _xInf.size()-1;++i)
                {
                    currentGap = _xInf[i+1]->getEval(NOMAD::EvalType::BB)->getFs()[obj] -
                        _xInf[i-1]->getEval(NOMAD::EvalType::BB)->getFs()[obj];
                    currentGap /= (fmax -fmin);
                    if (canBeFrameCenter[i] && currentGap >= maxGap)
                    {
                        maxGap = currentGap;
                        currentBestInd = i;
                    }
                }

                // Extreme points
                currentGap = 2 * (_xInf[_xInf.size() - 1]->getEval(NOMAD::EvalType::BB) ->getFs()[obj] -
                                  _xInf[_xInf.size() - 2] ->getEval(NOMAD::EvalType::BB) ->getFs()[obj]);
                currentGap /= (fmax - fmin);
                if (canBeFrameCenter[_xInf.size() - 1] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = _xInf.size()-1;
                }

                currentGap = 2 * (_xInf[1]->getEval(NOMAD::EvalType::BB)->getFs()[obj] -
                                  _xInf[0]->getEval(NOMAD::EvalType::BB)->getFs()[obj]);
                currentGap /= (fmax -fmin);
                if (canBeFrameCenter[0] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = 0;
                }
            }
            xInf = _xInf[currentBestInd];
        }
        // More than 2 objectives
        else
        {
            std::vector<std::pair<EvalPoint, size_t> > tmpXInfPInd(_xInf.size());

            // Initialize it.
            for (size_t i = 0; i < tmpXInfPInd.size(); ++i)
            {
                tmpXInfPInd[i] = std::make_pair(*_xInf[i], i);
            }

            size_t currentBestInd = 0;
            NOMAD::Double maxGap = -1.0;
            NOMAD::Double currentGap;

            for (size_t obj = 0; obj < _nobj; ++obj)
            {
                // Sort elements of tmpXInfPInd according to objective obj (in ascending order)
                std::sort(tmpXInfPInd.begin(), tmpXInfPInd.end(),
                          [obj](const std::pair<EvalPoint, size_t>& t1, const std::pair<EvalPoint, size_t> t2)->bool
                          {
                              const NOMAD::Eval* eval1 = t1.first.getEval(NOMAD::EvalType::BB);
                              const NOMAD::Eval* eval2 = t2.first.getEval(NOMAD::EvalType::BB);
                              return eval1->getFs()[obj] < eval2->getFs()[obj];
                          });

                // Get extreme values value according to one objective
                NOMAD::Double fmin = tmpXInfPInd[0].first.getEval(NOMAD::EvalType::BB)->getFs()[obj];
                NOMAD::Double fmax = tmpXInfPInd[tmpXInfPInd.size()-1].first.getEval(NOMAD::EvalType::BB)->getFs()[obj];

                // Can happen for exemple when we have several minima or for more than three objectives
                if (fmin == fmax)
                {
                    fmin = 0.0;
                    fmax = 1.0;
                }

                // Intermediate points
                for (size_t i = 1; i < tmpXInfPInd.size()-1;++i)
                {
                    currentGap = tmpXInfPInd[i+1].first.getEval(NOMAD::EvalType::BB)->getFs()[obj] -
                        tmpXInfPInd[i-1].first.getEval(NOMAD::EvalType::BB)->getFs()[obj];
                    currentGap /= (fmax -fmin);
                    if (canBeFrameCenter[tmpXInfPInd[i].second] && currentGap >= maxGap)
                    {
                        maxGap = currentGap;
                        currentBestInd = tmpXInfPInd[i].second;
                    }
                }

                // Extreme points
                currentGap = 2 * (tmpXInfPInd[tmpXInfPInd.size() - 1].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj] -
                                  tmpXInfPInd[tmpXInfPInd.size() - 2].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj]);
                currentGap /= (fmax - fmin);
                if (canBeFrameCenter[tmpXInfPInd[tmpXInfPInd.size()-1].second] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = tmpXInfPInd[tmpXInfPInd.size()-1].second;
                }

                currentGap = 2 * (tmpXInfPInd[1].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj] -
                                  tmpXInfPInd[0].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj]);
                currentGap /= (fmax -fmin);
                if (canBeFrameCenter[tmpXInfPInd[0].second] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = tmpXInfPInd[0].second;
                }
            }
            xInf = _xInf[currentBestInd];
        }

    }
    return xInf;
}


NOMAD::EvalPointPtr NOMAD::DMultiMadsBarrier::getXInfMinH() const
{
    size_t indXInfMinH = 0;
    NOMAD::Double hMinVal = NOMAD::INF;

    for (size_t i = 0; i < _xInf.size(); ++i)
    {
        const NOMAD::Eval* eval = _xInf[i]->getEval(NOMAD::EvalType::BB);
        NOMAD::Double h = eval->getH();

        // By definition, all elements of _xInf or _xFilterInf have a well-defined
        // h value. So, no need to check.
        if (h < hMinVal)
        {
            hMinVal = h;
            indXInfMinH = i;
        }
    }
    return _xInf[indXInfMinH];
}


void NOMAD::DMultiMadsBarrier::checkXInf(const NOMAD::EvalPoint &xInf, NOMAD::EvalType evalType)
{
    // If evalType is UNDEFINED, skip this check.
    if (NOMAD::EvalType::UNDEFINED != evalType)
    {
        if (nullptr == xInf.getEval(evalType))
        {
            throw NOMAD::Exception(__FILE__, __LINE__,
                                   "Barrier: xInf must be evaluated before being set.");
        }
    }
}



void NOMAD::DMultiMadsBarrier::clearXInf()
{
    _xInf.clear();
    _xFilterInf.clear();
    
    // Update the current incumbent inf. Only the infeasible one depends on XInf (not the case for the feasible one).
    updateCurrentIncumbentInf();
    
}


void NOMAD::DMultiMadsBarrier::updateRefBests()
{
    throw NOMAD::Exception(__FILE__, __LINE__, "Ref Bests are not used with DMultiMads.");
}


// The code is the very similar to what is in Barrier. Here we use current incumbents, but these points are accessible with the getFirstXFeas and getFirstXInf functions.
NOMAD::SuccessType NOMAD::DMultiMadsBarrier::getSuccessTypeOfPoints(const EvalPointPtr xFeas,
                                                            const EvalPointPtr xInf,
                                                            EvalType evalType,
                                                            ComputeType computeType)
{
    NOMAD::SuccessType successType = SuccessType::UNSUCCESSFUL;
    NOMAD::SuccessType successType2 = SuccessType::UNSUCCESSFUL;

    NOMAD::EvalPointPtr newBestFeas,newBestInf;


    // Set the iter success based on current incumbent feasible and infeasible points compared to initial point.
    if (nullptr != _currentIncumbentFeas || nullptr != _currentIncumbentInf)
    {
        // Compute success
        // Get which of newBestFeas and newBestInf is improving
        // the solution. Check newBestFeas first.
        NOMAD::ComputeSuccessType computeSuccess(evalType, computeType);

        if (nullptr != _currentIncumbentFeas)
        {
            successType = computeSuccess(xFeas, _currentIncumbentFeas,_hMax);
        }
        if (nullptr != _currentIncumbentInf)
        {
            successType2 = computeSuccess(xInf, _currentIncumbentInf,_hMax);
        }
        if (successType2 > successType)
        {
            successType = successType2;
        }
    }
    return successType;
}


bool NOMAD::DMultiMadsBarrier::updateFeasWithPoint(const EvalPoint & evalPoint,
                                               EvalType evalType,
                                               ComputeType computeType,
                                               const bool keepAllPoints)
{

    bool updated = false;

    auto eval = evalPoint.getEval(evalType);
    if (eval->isFeasible(computeType))
    {
        std::string s; // for output info
        
        OUTPUT_DEBUG_START
        s = "Point suggested to update DMultiMadsBarrier (feasible): " + evalPoint.display();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        
        if (eval->getFs(computeType).size() != _nobj)
        {
            s = "DMultiMadsBarrier update: number of objectives is equal to " + std::to_string(_nobj);
            s += ". Trying to add this point with number of objectives " + std::to_string(eval->getFs(computeType).size());
            s += ": " + evalPoint.display();
            throw NOMAD::Exception(__FILE__, __LINE__, s);
        }
        
        // Ensure evalPoint is as good as previous points in xFeas
        if (_xFeas.empty())
        {
            // New point is first point
            _xFeas.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
            updated = true;
            OUTPUT_DEBUG_START
            s = "New dominating xFeas: " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            
            // Update new best feasible point.
            _currentIncumbentFeas = _xFeas[0];
        }
        else
        {
            // Dealing with Non Standard computing order
            if (computeType != ComputeType::STANDARD)
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "Update using non standard compute type is not yet implemented");
            }
            else
            {
                bool insert = true;
                std::vector<bool> keepInXFeas(_xFeas.size(), true);
                int currentInd = 0;
                for (const auto& xFeas : _xFeas)
                {
                    auto compFlag = evalPoint.compMO(*xFeas, evalType);
                    if (compFlag == CompareType::DOMINATED)
                    {
                        OUTPUT_DEBUG_START
                        s = "evalPoint is dominated by " + xFeas->display() + "\n";
                        s += "evalPoint is rejected";
                        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                        OUTPUT_DEBUG_END
                        insert = false;
                        break;
                    }
                    else if (compFlag == CompareType::DOMINATING)
                    {
                        OUTPUT_DEBUG_START
                        s = "At least a point in xFeas is dominated by " + evalPoint.display();
                        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                        OUTPUT_DEBUG_END
                        updated = true;
                        keepInXFeas[currentInd] = false;
                    }
                    else if (compFlag == CompareType::EQUAL)
                    {
                        if (!keepAllPoints)
                        {
                            insert = false;
                            break;
                        }
                        
                        // If new point is not already there, add it
                        if (findEvalPoint(_xFeas.begin(), _xFeas.end(), evalPoint) != _xFeas.end())
                        {
                            OUTPUT_DEBUG_START
                            s = "EvalPoint " + evalPoint.display() + "is already in xFeas";
                            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                            OUTPUT_DEBUG_END
                            insert = false;
                        }
                        else
                        {
                            updated = true;
                        }
                        break;
                    }
                    currentInd++;
                }
                if (insert)
                {
                    // Remove all dominated elements.
                    currentInd = 0;
                    _xFeas.erase(std::remove_if(_xFeas.begin(), _xFeas.end(),
                                                [&currentInd, &keepInXFeas](const EvalPointPtr evalPoint)
                        {
                            bool isRemoved = !keepInXFeas[currentInd];
                            ++currentInd;
                            return isRemoved;
                        }), _xFeas.end());
                    
                    // Update mesh and insert the new element.
                    OUTPUT_DEBUG_START
                    s = "Adding new non dominated eval point in xFeas: " + evalPoint.display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    s = "Update the mesh of eval point to be added";
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                    updated = true;
                    
                    // Update the mesh -> success
                    auto dir = evalPoint.getDirection();
                    if (nullptr != dir)
                    {
                        evalPoint.getMesh()->enlargeDeltaFrameSize(*dir);
                    }
                    
                    // insert new element
                    _xFeas.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
                    
                    // Sort according to lexicographic order.
                    std::sort(_xFeas.begin(), _xFeas.end(),
                              [&evalType](const EvalPointPtr evalPoint1, const EvalPointPtr evalPoint2)
                              {
                        const NOMAD::Eval* eval1 = evalPoint1->getEval(evalType);
                        const NOMAD::Eval* eval2 = evalPoint2->getEval(evalType);
                        return eval1->getFs().lexicographicalCmp(eval2->getFs());
                    });
                }
            }
        }
    }
    return updated;
}

bool NOMAD::DMultiMadsBarrier::updateInfWithPoint(const EvalPoint & evalPoint,
                                               EvalType evalType,
                                               ComputeType computeType,
                                               const bool keepAllPoints,
                                               const bool feasHasBeenUpdated)
{
    
    bool updated = false;
    
    auto eval = evalPoint.getEval(evalType);
    if (!eval->isFeasible(computeType))
    {
        std::string s; // for output info
    
        
        NOMAD::Double h = eval->getH(computeType);
        if (!h.isDefined() || h == NOMAD::INF ||
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
            return false;
        }
        
        if (_xInf.empty())
        {
            // New point is first point
            _xInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
            _xFilterInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
            _currentIncumbentInf = _xInf[0];
            updated = true;
            OUTPUT_DEBUG_START
            s = "New current incumbent infeasible: " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
        }
        // Insertion into the two sets of infeasible non dominated points.
        else
        {
            // Try to insert into _xInfFilter.
            bool insert = true;
            std::vector<bool> isInXinfFilter(_xFilterInf.size(), true);
            int currentInd = 0;
            
            for (const auto& xFilterInf : _xFilterInf)
            {
                auto compFlag = evalPoint.compMO(*xFilterInf, evalType);
                if (compFlag == CompareType::DOMINATED)
                {
                    insert = false;
                    break;
                }
                else if (compFlag == CompareType::DOMINATING)
                {
                    OUTPUT_DEBUG_START
                    s = "xInf dominates the filter element " + xFilterInf->display();
                    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                    OUTPUT_DEBUG_END
                    updated = true;
                    isInXinfFilter[currentInd] = false;
                }
                else if (compFlag == CompareType::EQUAL)
                {
                    if (!keepAllPoints)
                    {
                        insert = false;
                        break;
                    }
                    
                    // If new point is not already there, add it.
                    if (findEvalPoint(_xFilterInf.begin(), _xFilterInf.end(), evalPoint) != _xFilterInf.end())
                    {
                        OUTPUT_DEBUG_START
                        s = "xInf is already here: " + evalPoint.display();
                        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                        OUTPUT_DEBUG_END
                        insert = false;
                    }
                    else
                    {
                        updated = true;
                        break;
                    }
                }
                currentInd++;
            }
            
            if (insert)
            {
                // Remove all dominated elements of _xInfFilter
                currentInd = 0;
                _xFilterInf.erase(std::remove_if(_xFilterInf.begin(), _xFilterInf.end(),
                                                 [&currentInd, &isInXinfFilter](const EvalPointPtr ev)
                                                 {
                    bool isRemoved = !isInXinfFilter[currentInd];
                    currentInd++;
                    return isRemoved;
                }), _xFilterInf.end());
                _xFilterInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
                
                // Order by lexicographic ordering.
                std::sort(_xFilterInf.begin(), _xFilterInf.end(),
                          [&evalType](const EvalPointPtr evalPoint1, const EvalPointPtr evalPoint2)->bool
                          {
                    const NOMAD::Eval* eval1 = evalPoint1->getEval(evalType);
                    const NOMAD::Eval* eval2 = evalPoint2->getEval(evalType);
                    return eval1->getFs().lexicographicalCmp(eval2->getFs());
                });
                
                // Try to insert the point into _xInf.
                insert = true;
                
                currentInd = 0;
                std::vector<bool> isInXinf(_xInf.size(), true);
                
                for (const auto& xInf : _xInf)
                {
                    auto compFlag = evalPoint.compMO(*xInf, evalType, true);
                    if (compFlag == CompareType::DOMINATED)
                    {
                        insert = false;
                        break;
                    }
                    // The points could be only equal according to the f values
                    else if ((compFlag == CompareType::DOMINATING) ||
                             (evalPoint.compMO(*xInf, evalType) == CompareType::DOMINATING))
                    {
                        OUTPUT_DEBUG_START
                        s = "xInf dominates the filter element " + xInf->display();
                        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                        OUTPUT_DEBUG_END
                        updated = true;
                        isInXinf[currentInd] = false;
                    }
                    currentInd++;
                }
                
                if (insert)
                {
                    // Remove all dominated elements of _xInf.
                    currentInd = 0;
                    _xInf.erase(std::remove_if(_xInf.begin(), _xInf.end(),
                                               [&currentInd, &isInXinf](const EvalPointPtr ev)
                                               {
                        bool isRemoved = !isInXinf[currentInd];
                        currentInd++;
                        return isRemoved;
                    }), _xInf.end());
                    updated = true;
                    _xInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
                    
                    
                    std::sort(_xInf.begin(), _xInf.end(),
                              [&evalType](const EvalPointPtr evalPoint1, const EvalPointPtr evalPoint2)->bool
                              {
                        const NOMAD::Eval* eval1 = evalPoint1->getEval(evalType);
                        const NOMAD::Eval* eval2 = evalPoint2->getEval(evalType);
                        return eval1->getFs().lexicographicalCmp(eval2->getFs());
                    });
                    
                }
            }
        }
    }
    return updated;
    
}


// Update the barrier (feas and inf). Once done calls for update the current best feas and inf.
bool NOMAD::DMultiMadsBarrier::updateWithPoints(
                                        const std::vector<EvalPoint>& evalPointList,
                                        EvalType evalType,
                                        ComputeType computeType,
                                        const bool keepAllPoints)
{
    bool updated = false;
    bool updatedFeas = false;
    bool updatedInf = false;
    std::string s; // for output info

    // Temporary infeasible incumbent. For insertion of more than one improving/full success
    NOMAD::EvalPoint xInfTmp;

    OUTPUT_DEBUG_START
    s = "Updating DMultiMadsBarrier (" + std::to_string(_nobj) + " objectives)";
    s += " with " + std::to_string(evalPointList.size()) + " suggested points";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    s = "Current barrier: ";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    std::vector<std::string> vs = display(4);
    for (const auto & s: vs)
    {
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    OUTPUT_DEBUG_END

    // do separate loop on evalPointList
    // First loop update the bestFeasible.
    // If a point is a full success setupdatedFeas = true.
    // This flag is used in the second loop.
    for (const auto & evalPoint : evalPointList)
    {
        // All points must have mesh parameters defined.
        checkMeshParameters(evalPoint);

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
        
        if (computeType == ComputeType::STANDARD && eval->getFs(computeType).size() != _nobj)
        {
            s = "DMultiMadsBarrier update: number of objectives is equal to " + std::to_string(_nobj);
            s += ". Trying to add this point with number of objectives " + std::to_string(eval->getFs(computeType).size());
            s += ": " + evalPoint.display();
            throw NOMAD::Exception(__FILE__, __LINE__, s);
        }
        
        updatedFeas = updateFeasWithPoint(evalPoint,evalType, computeType, keepAllPoints) || updatedFeas ;

    }

    // Do separate loop on evalPointList
    // Second loop update the bestInfeasible.
    // Use the flag oneFeasEvalFullSuccess.
    // If the flag is true hmax will not change. A point improving the best infeasible should not replace it.

    for (const auto & evalPoint : evalPointList)
    {
        auto eval = evalPoint.getEval(evalType);

        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            // Suggested point is not good.
            continue;
        }
        
        updatedInf = updateInfWithPoint(evalPoint, evalType, computeType, keepAllPoints, updatedFeas) || updatedInf ;

    }

    updated = updated || updatedFeas || updatedInf;
    
    if (updated)
    {
        // Set n and check that all points have the same dimension
        setN();
        
        // Update incumbents
        updateCurrentIncumbents();
        
    }

    
    
    OUTPUT_DEBUG_START
    if (updated)
    {
        s = "Updated barrier";
        if (updatedFeas && updatedInf)
        {
            s += " with feasible and infeasible points";
        }
        else if (updatedFeas)
        {
            s += " with feasible points";
        }
        else if (updatedInf)
        {
            s += " with infeasible points";
        }
    }
    else
    {
        s = "Barrier not updated";
    }
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    return updated;
}


void NOMAD::DMultiMadsBarrier::updateXInfAndFilterInfAfterHMaxSet()
{
    if (_xInf.size() == 0)
    {
        return;
    }

    size_t currentInd = 0;

    // Remove all infeasible incumbent solutions below the threshold _hMax.
    std::vector<bool> isInXInf(_xInf.size(), true);
    for (const auto& xInf : _xInf)
    {
        const NOMAD::Eval* eval = xInf->getEval(NOMAD::EvalType::BB);
        NOMAD::Double h = eval->getH();

        if (h > _hMax)
        {
            isInXInf[currentInd] = false;

        }
        currentInd++;
    }

    currentInd = 0;
    _xInf.erase(std::remove_if(_xInf.begin(), _xInf.end(),
                               [&currentInd, &isInXInf](const EvalPointPtr evalPoint)
                               {
                                   bool isRemoved = !isInXInf[currentInd];
                                   ++currentInd;
                                   return isRemoved;
                               }), _xInf.end());

    currentInd = 0;
    // Remove all infeasible non dominated solutions below the threshold _hMax.
    std::vector<bool> isInXFilterInf(_xFilterInf.size(), true);
    for (const auto& xFilterInf : _xFilterInf)
    {
        const NOMAD::Eval* eval = xFilterInf->getEval(NOMAD::EvalType::BB);
        NOMAD::Double h = eval->getH();

        if (h > _hMax)
        {
            isInXFilterInf[currentInd] = false;

        }
        currentInd++;
    }

    currentInd = 0;
    _xFilterInf.erase(std::remove_if(_xFilterInf.begin(), _xFilterInf.end(),
                                     [&currentInd, &isInXFilterInf](const EvalPointPtr evalPoint)
                                     {
                                         bool isRemoved = !isInXFilterInf[currentInd];
                                         ++currentInd;
                                         return isRemoved;
                                     }), _xFilterInf.end());

    std::sort(_xFilterInf.begin(), _xFilterInf.end(),
              [](const EvalPointPtr evalPoint1, const EvalPointPtr evalPoint2) -> bool
              {
                  const NOMAD::Eval* eval1 = evalPoint1->getEval(NOMAD::EvalType::BB);
                  const NOMAD::Eval* eval2 = evalPoint2->getEval(NOMAD::EvalType::BB);
                  return eval1->getFs().lexicographicalCmp(eval2->getFs());
              });

    // And reinsert potential infeasible non dominated points into the set of infeasible
    // solutions.
    currentInd = 0;
    isInXInf = std::vector<bool>(_xFilterInf.size(), false);
    for (const auto& evalPoint : _xFilterInf)
    {
        // If new point is not in _xInf, check if it is not dominated
        if (findEvalPoint(_xInf.begin(), _xInf.end(), *evalPoint) == _xInf.end())
        {
            size_t currentIndTmp = 0;
            bool insert = true;
            for (const auto& evalPointInf : _xFilterInf)
            {
                if (currentIndTmp != currentInd)
                {
                    auto compFlag = evalPoint->compMO(*evalPointInf, NOMAD::EvalType::BB, true);
                    if (compFlag == CompareType::DOMINATED)
                    {
                        insert = false;
                        break;
                    }
                    else if (compFlag == CompareType::DOMINATING)
                    {
                        isInXInf[currentIndTmp] = false;
                    }
                }
                currentIndTmp++;
            }
            isInXInf[currentInd] = insert;
        }
        currentInd++;
    }

    for (size_t i = 0; i < isInXInf.size(); ++i)
    {
        if (isInXInf[i])
        {
            _xInf.push_back(_xFilterInf[i]);
        }
    }

    std::sort(_xInf.begin(), _xInf.end(),
              [](const EvalPointPtr evalPoint1, const EvalPointPtr evalPoint2) -> bool
              {
                  const NOMAD::Eval* eval1 = evalPoint1->getEval(NOMAD::EvalType::BB);
                  const NOMAD::Eval* eval2 = evalPoint2->getEval(NOMAD::EvalType::BB);
                  return eval1->getFs().lexicographicalCmp(eval2->getFs());
              });
    return;
}


// Nice formatting display.
std::vector<std::string> NOMAD::DMultiMadsBarrier::display(const size_t max) const
{
    std::vector<std::string> vs;

    size_t nbXFeas = 0;
    size_t nbXInf = 0;
    size_t nbXFilterInf = 0;

    for (const auto & xFeas : _xFeas)
    {
        vs.push_back("X_FEAS " + xFeas->displayAll());
        nbXFeas++;
        if (nbXFeas >= max && _xFeas.size() > max)
        {
            vs.push_back("... (total " + std::to_string(_xFeas.size()) + ")");
            break;
        }
    }
    for (const auto &xInf : _xInf)
    {
        vs.push_back("X_INF " + xInf->displayAll() );
        nbXInf++;

        if (nbXInf >= max && _xInf.size() > max)
        {
            vs.push_back("... (total " + std::to_string(_xInf.size()) + ")");
            break;
        }
    }
    for (const auto &xFilterInf: _xFilterInf)
    {
        vs.push_back("X_FILTER" + xFilterInf->displayAll());
        nbXFilterInf++;

        if (nbXFilterInf >= max && _xFilterInf.size() > max)
        {
            vs.push_back("... (total " + std::to_string(_xFilterInf.size()) + ")");
            break;
        }
    }

    vs.push_back("H_MAX " + getHMax().display(NOMAD::DISPLAY_PRECISION_FULL));
    vs.push_back("Ref Best Feasible:   " + (_refBestFeas ? _refBestFeas->displayAll() : "NULL") );
    vs.push_back("Ref Best Infeasible: " + (_refBestInf ? _refBestInf->displayAll() : "NULL"));

    return vs;
}



NOMAD::Double NOMAD::DMultiMadsBarrier::getMeshMaxFrameSize(const NOMAD::EvalPointPtr pt) const
{
    NOMAD::Double maxRealVal = -1.0;
    NOMAD::Double maxIntegerVal = -1.0;

    auto mesh = pt->getMesh();
    
    // Detect if mesh is sub dimension and pt are in full dimension.
    bool meshIsInSubDimension = false;
    if (pt->getMesh()->getSize() < pt->size())
    {
        meshIsInSubDimension = true;
    }
    
    size_t shift = 0;
    for (size_t i = 0; i < pt->size(); ++i)
    {
        // Do not use access the frame size for fixed variables.
        if (meshIsInSubDimension && _fixedVariables[i].isDefined())
        {
            shift++;
            continue;
        }
        if (_bbInputsType[i] == NOMAD::BBInputType::CONTINUOUS)
        {
            maxRealVal = std::max(maxRealVal, mesh->getDeltaFrameSize(i-shift));
        }
        else if (_bbInputsType[i] == NOMAD::BBInputType::INTEGER)
        {
            maxIntegerVal = std::max(maxIntegerVal, mesh->getDeltaFrameSize(i-shift));
        }
    }
    if (maxRealVal > 0.0)
    {
        return maxRealVal; // Some values are real: get norm inf on these values only.
    }
    else if (maxIntegerVal > 0.0)
    {
        return maxIntegerVal; // No real value but some integer values: get norm inf
        // on these values only
    }
    else
    {
        return 1.0; // Only binary variables: any elements of the iterate list can be
        // chosen;
    }

}

void NOMAD::DMultiMadsBarrier::updateCurrentIncumbentInf()
{
    _currentIncumbentInf = nullptr;
    if (_xFeas.size() > 0 && _xInf.size() > 0)
    {
        // Get the infeasible solution with maximum dominance move below the _hMax threshold,
        // according to the set of best feasible incumbent solutions.
        size_t currentInd = 0;
        double maxDomMove = -NOMAD::INF;

        for (size_t j = 0; j < _xInf.size(); ++j)
        {
            // Compute dominance move
            // = min \sum_{1}^m max(fi(y) - fi(x), 0)
            //   y \in Fk
            double tmpDomMove = NOMAD::INF;
            const NOMAD::Eval* evalInf = _xInf[j]->getEval(NOMAD::EvalType::BB);
            NOMAD::Double h = evalInf->getH();

            if (h.isDefined() && h <= _hMax)
            {
                for (const auto &xFeas : _xFeas)
                {

                    double sumVal = 0.0;
                    const NOMAD::Eval *evalFeas = xFeas->getEval(NOMAD::EvalType::BB);

                    // Compute \sum_{1}^m max (fi(y) - fi(x), 0)
                    for (size_t i = 0; i < _nobj; i++)
                    {
                        sumVal += std::max(evalFeas->getFs()[i].todouble() -
                                           evalInf->getFs()[i].todouble(),
                                           0.0);
                    }
                    if (tmpDomMove > sumVal)
                    {
                        tmpDomMove = sumVal;
                    }
                }

                // Get the maximum dominance move index
                if (maxDomMove < tmpDomMove) {
                    maxDomMove = tmpDomMove;
                    currentInd = j;
                }
            }
        }

        // In this case, all infeasible solutions are "dominated" in terms of fvalues
        // by at least one element of Fk
        if (NOMAD::Double(maxDomMove) == 0.0)
        {
            // In this case, get the infeasible solution below the _hMax threshold which has
            // minimal dominance move, when considered a maximization problem.
            double minDomMove = NOMAD::INF;
            currentInd = 0;

            for (size_t j = 0; j < _xInf.size(); ++j)
            {
                // Compute dominance move
                // = min \sum_{1}^m max(fi(x) - fi(y), 0)
                //   y \in Fk
                double tmpDomMove = NOMAD::INF;
                const NOMAD::Eval* evalInf = _xInf[j]->getEval(NOMAD::EvalType::BB);

                NOMAD::Double h = evalInf->getH();
                if (h.isDefined() && h <= _hMax)
                {
                    for (const auto& xFeas : _xFeas)
                    {

                        double sumVal = 0.0;
                        const NOMAD::Eval* evalFeas = xFeas->getEval(NOMAD::EvalType::BB);

                        // Compute \sum_{1}^m max (fi(x) - fi(y), 0)
                        for (size_t i = 0; i < _nobj; i++)
                        {
                            sumVal += std::max(evalInf->getFs()[i].todouble() -
                                               evalFeas->getFs()[i].todouble(), 0.0);
                        }
                        if (tmpDomMove > sumVal)
                        {
                            tmpDomMove = sumVal;
                        }
                    }

                    // Get the minimal dominance move index
                    if (minDomMove > tmpDomMove)
                    {
                        minDomMove = tmpDomMove;
                        currentInd = j;
                    }
                }
            }
        }
        _currentIncumbentInf = _xInf[currentInd];
    }
    else
    {
        _currentIncumbentInf = getFirstXInfNoXFeas();
    }

}


void NOMAD::DMultiMadsBarrier::updateCurrentIncumbentFeas()
{
    if (_xFeas.size() == 0)
    {
        _currentIncumbentFeas = nullptr;
        return;
    }
    
    if (_xFeas.size() == 1)
    {
        _currentIncumbentFeas =_xFeas[0];
        return ;
        
    }
    
    // Set max frame size of all elements
    NOMAD::Double maxFrameSizeFeasElts = -1.0;
    for (const auto& xFeas: _xFeas)
    {
        maxFrameSizeFeasElts = std::max(getMeshMaxFrameSize(xFeas), maxFrameSizeFeasElts);
    }
    
    // Select candidates
    std::vector<bool> canBeFrameCenter(_xFeas.size(), false);
    size_t nbSelectedCandidates = 0;
    
    // see article DMultiMads Algorithm 4.
    for (size_t i = 0; i < _xFeas.size(); ++i)
    {
        
        NOMAD::Double maxFrameSizeElt = getMeshMaxFrameSize(_xFeas[i]);
        
        // Casting is required to avoid an integer overflow.
        if ((std::pow(10.0, -(double)_incumbentSelectionParam) * maxFrameSizeFeasElts <= maxFrameSizeElt) )
        {
            canBeFrameCenter[i] = true;
            nbSelectedCandidates += 1;
        }
    }
    
    // Only one point in the barrier.
    if (nbSelectedCandidates == 1)
    {
        auto it = std::find(canBeFrameCenter.begin(), canBeFrameCenter.end(), true);
        if (it == canBeFrameCenter.end())
        {
            std::string s = "Error: DMultiMadsBarrier, should not reach this condition";
            throw NOMAD::Exception(__FILE__,__LINE__,s);
        }
        else
        {
            size_t selectedInd = std::distance(canBeFrameCenter.begin(), it);
            _currentIncumbentFeas =_xFeas[selectedInd];
        }
    }
    // Only two points in the barrier.
    else if ((nbSelectedCandidates == 2) && (_xFeas.size() == 2))
    {
        const NOMAD::Eval* eval1 = _xFeas[0]->getEval(NOMAD::EvalType::BB);
        const NOMAD::Eval* eval2 = _xFeas[1]->getEval(NOMAD::EvalType::BB);
        auto objv1 = eval1->getFs();
        auto objv2 = eval2->getFs();
        if (objv1.abs().max() > objv2.abs().max())
        {
            _currentIncumbentFeas =_xFeas[0];
        }
        else
        {
            _currentIncumbentFeas =_xFeas[1];
        }
    }
    // More than three points in the barrier.
    else
    {
        // First case: biobjective optimization. Points are already ranked by lexicographic order.
        if (_nobj == 2)
        {
            size_t currentBestInd = 0;
            NOMAD::Double maxGap = -1.0;
            NOMAD::Double currentGap;
            for (size_t obj = 0; obj < _nobj; ++obj)
            {
                // Get extreme values value according to one objective
                NOMAD::Double fmin = _xFeas[0]->getEval(NOMAD::EvalType::BB)->getFs()[obj];
                NOMAD::Double fmax = _xFeas[_xFeas.size()-1]->getEval(NOMAD::EvalType::BB)->getFs()[obj];
                
                // In this case, it means all elements of _xFeas are equal (return the first one)
                if (fmin == fmax)
                {
                    break;
                }
                
                // Intermediate points
                for (size_t i = 1; i < _xFeas.size()-1;++i)
                {
                    currentGap = _xFeas[i+1]->getEval(NOMAD::EvalType::BB)->getFs()[obj] -
                    _xFeas[i-1]->getEval(NOMAD::EvalType::BB)->getFs()[obj];
                    currentGap /= (fmax -fmin);
                    if (canBeFrameCenter[i] && currentGap >= maxGap)
                    {
                        maxGap = currentGap;
                        currentBestInd = i;
                    }
                }
                
                // Extreme points
                currentGap = 2 * (_xFeas[_xFeas.size() - 1]->getEval(NOMAD::EvalType::BB) ->getFs()[obj] -
                                  _xFeas[_xFeas.size() - 2] ->getEval(NOMAD::EvalType::BB) ->getFs()[obj]);
                currentGap /= (fmax - fmin);
                if (canBeFrameCenter[_xFeas.size()-1] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = _xFeas.size()-1;
                }
                
                currentGap = 2 * (_xFeas[1]->getEval(NOMAD::EvalType::BB)->getFs()[obj] -
                                  _xFeas[0]->getEval(NOMAD::EvalType::BB)->getFs()[obj]);
                currentGap /= (fmax -fmin);
                if (canBeFrameCenter[0] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = 0;
                }
            }
            _currentIncumbentFeas =_xFeas[currentBestInd];
        }
        // More than 2 objectives
        else
        {
            std::vector<std::pair<EvalPoint, size_t> > tmpXFeasPInd(_xFeas.size());
            
            // Initialize it.
            for (size_t i = 0; i < tmpXFeasPInd.size(); ++i)
            {
                tmpXFeasPInd[i] = std::make_pair(*_xFeas[i], i);
            }

            size_t currentBestInd = 0;
            NOMAD::Double maxGap = -1.0;
            NOMAD::Double currentGap;
            
            for (size_t obj = 0; obj < _nobj; ++obj)
            {
                // Sort elements of tmpXFeasPInd according to objective obj (in ascending order)
                std::sort(tmpXFeasPInd.begin(), tmpXFeasPInd.end(),
                          [obj](const std::pair<EvalPoint, size_t>& t1, const std::pair<EvalPoint, size_t> t2)->bool
                          {
                    const NOMAD::Eval* eval1 = t1.first.getEval(NOMAD::EvalType::BB);
                    const NOMAD::Eval* eval2 = t2.first.getEval(NOMAD::EvalType::BB);
                    return eval1->getFs()[obj] < eval2->getFs()[obj];
                });
                
                // Get extreme values value according to one objective
                NOMAD::Double fmin = tmpXFeasPInd[0].first.getEval(NOMAD::EvalType::BB)->getFs()[obj];
                NOMAD::Double fmax = tmpXFeasPInd[tmpXFeasPInd.size()-1].first.getEval(NOMAD::EvalType::BB)->getFs()[obj];
                
                // Can happen for exemple when we have several minima or for more than three objectives
                if (fmin == fmax)
                {
                    fmin = 0.0;
                    fmax = 1.0;
                }
                
                // Intermediate points
                for (size_t i = 1; i < tmpXFeasPInd.size()-1;++i)
                {
                    currentGap = tmpXFeasPInd[i+1].first.getEval(NOMAD::EvalType::BB)->getFs()[obj] -
                    tmpXFeasPInd[i-1].first.getEval(NOMAD::EvalType::BB)->getFs()[obj];
                    currentGap /= (fmax -fmin);
                    if (canBeFrameCenter[tmpXFeasPInd[i].second] && currentGap >= maxGap)
                    {
                        maxGap = currentGap;
                        currentBestInd = tmpXFeasPInd[i].second;
                    }
                }
                
                // Extreme points
                currentGap = 2 * (tmpXFeasPInd[tmpXFeasPInd.size() - 1].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj] -
                                  tmpXFeasPInd[tmpXFeasPInd.size() - 2].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj]);
                currentGap /= (fmax - fmin);
                if (canBeFrameCenter[tmpXFeasPInd[tmpXFeasPInd.size()-1].second] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = tmpXFeasPInd[tmpXFeasPInd.size()-1].second;
                }
                
                currentGap = 2 * (tmpXFeasPInd[1].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj] -
                                  tmpXFeasPInd[0].first.getEval(NOMAD::EvalType::BB) ->getFs()[obj]);
                currentGap /= (fmax -fmin);
                if (canBeFrameCenter[tmpXFeasPInd[0].second] && currentGap > maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = tmpXFeasPInd[0].second;
                }
            }
            _currentIncumbentFeas = _xFeas[currentBestInd];
        }
        
    }
    
    
}
