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
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Cache/CacheSet.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Eval/MeshBase.hpp"
#include "../../Output/OutputQueue.hpp"

// Initialization from cache
void NOMAD::DMultiMadsBarrier::init(const NOMAD::Point& fixedVariables,
                                    bool barrierInitializedFromCache)
{
    std::vector<NOMAD::EvalPoint> cachePoints;
    
    
    if (_computeType.evalType != NOMAD::EvalType::BB || _computeType.fhComputeTypeS.computeType != NOMAD::ComputeType::STANDARD )
    {
        std::string s = "Error: Eval type must be BB and Compute Type must be standard";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    if (fixedVariables.isEmpty())
    {
        std::string s = "Error: Fixed variable of dimension 0";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
    if (barrierInitializedFromCache)
    {
        checkCache();

        // Order lexicographically a set of non dominated points.
        auto lexicographicallyObjOrder = [](const NOMAD::FHComputeType& computeType,
                                            std::vector<EvalPointPtr>& xSet)
        {
            std::sort(xSet.begin(), xSet.end(),
                      [computeType](const EvalPointPtr& evalPoint1, const EvalPointPtr& evalPoint2)->bool
                      {
                          return evalPoint1->getFs(computeType).lexicographicalCmp(evalPoint2->getFs(computeType));
                      });
        };

        // Get the best feasible and infeasible solutions from cache.
        // Point from cache are in full dimension.
        // Convert them to subproblem dimension.
        // NB: all solutions (and not only the non dominated ones are considered)
        auto cache = CacheBase::getInstance().get();
        if (cache->findBestFeas(cachePoints, fixedVariables, _computeType) > 0)
        {
            for (const auto & evalPoint : cachePoints)
            {
                NOMAD::EvalPointPtr evalPointSub = std::make_shared<EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariables));
                _xFeas.push_back(evalPointSub);
            }
            cachePoints.clear();
            lexicographicallyObjOrder(_computeType, _xFeas);
        }
        if (cache->findFilterInf(cachePoints, _hMax, fixedVariables, _computeType) > 0)
        {
            for (const auto &evalPoint : cachePoints)
            {
                // Consider points with h < INF. That is, points that are not excluded by extreme barrier constraints.
                if (evalPoint.getH(_computeType) < NOMAD::INF)
                {
                    NOMAD::EvalPointPtr evalPointSub = std::make_shared<EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariables));
                    _xInf.push_back(evalPointSub);
                }
            }
            cachePoints.clear();
            lexicographicallyObjOrder(_computeType, _xInf);
        }

        // Get non dominated infeasible points from cache.
        // Points from cache are in full dimension.
        // Convert them to subproblem dimension.
        if (cache->findFilterInf(cachePoints, _hMax, fixedVariables, _computeType) > 0)
        {
            for (const auto & evalPoint : cachePoints)
            {
                NOMAD::EvalPointPtr evalPointSub = std::make_shared<EvalPoint>( evalPoint.makeSubSpacePointFromFixed(fixedVariables));
                _xFilterInf.push_back(evalPointSub);
            }
            cachePoints.clear();
            lexicographicallyObjOrder(_computeType, _xFilterInf);
        }
    }

    if (!_xFeas.empty() || !_xInf.empty())
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
                                    const std::vector<NOMAD::EvalPoint>& evalPointList)
{
    bool updated = updateWithPoints(evalPointList, true); // All points are considered.

    if (! updated)
         return;
    
    // Need to set bb input type
    // Initialize _bbInputsType. By default, it is a vector of real.
    if (_bbInputsType.empty())
    {
        _bbInputsType = BBInputTypeList(_n, NOMAD::BBInputType::CONTINUOUS);
    }
    
    if (_computeType.evalType != NOMAD::EvalType::BB || _computeType.fhComputeTypeS.computeType != NOMAD::ComputeType::STANDARD )
    {
        std::string s = "Error: Eval type must be BB and Compute Type must be standard";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }
    
    if (_bbInputsType.size() != _n)
    {
        std::string s = "Error: Inputs dimensions of DMultiMadsBarrier do not match dimensions of provided input types.";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    // Check: xFeas or xInf could be non-evaluated, but not both.
    if (   (_xFeas.empty() || nullptr == _xFeas[0]->getEval(_computeType.evalType))
        && (_xInf.empty() || nullptr == _xInf[0]->getEval(_computeType.evalType)))
    {
        std::string s = "Barrier constructor: xFeas or xInf must be in the barrier.\n";
        if (!_xFeas.empty())
        {
            s += "There are " + std::to_string(_xFeas.size()) + " xFeas, the first one is:\n";
            s += _xFeas[0]->displayAll(NOMAD::defaultFHComputeTypeS);
        }
        if (!_xInf.empty())
        {
            s += "There are " + std::to_string(_xInf.size()) + " xInf, the first one is:\n";
            s += _xInf[0]->displayAll(NOMAD::defaultFHComputeTypeS);
        }
        if (_xFeas.empty() && _xInf.empty())
        {
            s += "There are no xFeas and no xInf defined.";
        }
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    checkHMax();

    // Update the incumbents used by DMultiMads algo as frame centers and bounds.
    updateCurrentIncumbents();
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
    // Always update the current infeasible incumbents after changing hmax.
    updateCurrentIncumbentInfMaxH();
    updateCurrentIncumbentInf();
    updateCurrentIdealInf();
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
    
    // Update the current incumbents. Both the feasible and infeasible ones.
    updateCurrentIncumbents();
    updateCurrentIdealFeas();
}



void NOMAD::DMultiMadsBarrier::checkXFeasIsFeas(const NOMAD::EvalPoint &xFeas)
{
    const auto evalType = _computeType.evalType;

    // If evalType is UNDEFINED, skip this check.
    if (evalType == NOMAD::EvalType::UNDEFINED)
        return;

    const auto eval = xFeas.getEval(evalType);
    if (eval == nullptr || eval->getEvalStatus() != NOMAD::EvalStatusType::EVAL_OK)
        return;

    const auto computeTypeS = _computeType.Short();
    const NOMAD::Double h = eval->getH(computeTypeS);
    if (!h.isDefined() || 0.0 != h)
    {
        std::string err = "Error: DMultiMadsBarrier: xFeas' h value must be 0.0, got: " + h.display();
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }

    if (computeTypeS.computeType == NOMAD::ComputeType::STANDARD && eval->getFs(computeTypeS).size() != _nobj)
    {
        std::string err = "Error: DMultiMadsBarrier: xFeas' F must be of size " + std::to_string(_nobj);
        err += ", got: F.size() = " + std::to_string(eval->getFs(computeTypeS).size());
        err += " with following F values " + eval->getFs(computeTypeS).display();
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }
}


NOMAD::EvalPointPtr NOMAD::DMultiMadsBarrier::getFirstXIncInfNoXFeas() const
{
    if (_xFilterInf.empty())
    {
        return nullptr;
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
        return _xInf[0];
    }

    if (nbSelectedCandidates == 1)
    {
        auto it = std::find(canBeFrameCenter.begin(), canBeFrameCenter.end(), true);
        if (it == canBeFrameCenter.end())
        {
            std::string s = "Error: DMultiMadsBarrier, should not reach this condition";
            throw NOMAD::Exception(__FILE__,__LINE__,s);
        }
        const size_t selectedInd = std::distance(canBeFrameCenter.begin(), it);
        return _xInf[selectedInd];
    }

    if ((nbSelectedCandidates == 2) && (_xInf.size() == 2))
    {
        const auto& objv1 = _xInf[0]->getFs(_computeType);
        const auto& objv2 = _xInf[1]->getFs(_computeType);
        const EvalPointPtr xInf = (objv1.abs().max() > objv2.abs().max()) ? _xInf[0] : _xInf[1];
        return xInf;
    }

    // More than two points in the barrier.
    // First case: biobjective optimization. Points are already ranked by lexicographic order.
    if (_nobj == 2)
    {
        size_t currentBestInd = 0;
        NOMAD::Double maxGap = -1.0;
        NOMAD::Double currentGap;
        for (size_t obj = 0; obj < _nobj; ++obj)
        {
            // Get extreme values value according to one objective
            const NOMAD::Double fmin = _xInf[0]->getFs(_computeType)[obj];
            const NOMAD::Double fmax = _xInf[_xInf.size()-1]->getFs(_computeType)[obj];

            // In this case, it means all elements of _xInf are equal.
            // We return the first one.
            if (fmin == fmax)
            {
                return _xInf[0];
            }

            // Intermediate points
            for (size_t i = 1; i < _xInf.size()-1;++i)
            {
                if (!canBeFrameCenter[i])
                    continue;

                currentGap = _xInf[i+1]->getFs(_computeType)[obj] -
                    _xInf[i-1]->getFs(_computeType)[obj];
                currentGap /= (fmax -fmin);
                if (currentGap >= maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = i;
                }
            }

            // Extreme points
            currentGap = 2 * (_xInf[_xInf.size() - 1]->getFs(_computeType)[obj] -
                              _xInf[_xInf.size() - 2]->getFs(_computeType)[obj]);
            currentGap /= (fmax - fmin);
            if (canBeFrameCenter[_xInf.size() - 1] && currentGap > maxGap)
            {
                maxGap = currentGap;
                currentBestInd = _xInf.size() - 1;
            }

            currentGap = 2 * (_xInf[1]->getFs(_computeType)[obj] -
                              _xInf[0]->getFs(_computeType)[obj]);
            currentGap /= (fmax -fmin);
            if (canBeFrameCenter[0] && currentGap > maxGap)
            {
                maxGap = currentGap;
                currentBestInd = 0;
            }
        }
        return _xInf[currentBestInd];
    }

    // Case 2 : more than 2 objectives
    std::vector<std::pair<EvalPointPtr, size_t> > tmpXInfPInd(_xInf.size());

    // Initialize it.
    for (size_t i = 0; i < tmpXInfPInd.size(); ++i)
    {
        tmpXInfPInd[i] = std::make_pair(_xInf[i], i);
    }

    size_t currentBestInd = 0;
    NOMAD::Double maxGap = -1.0;
    NOMAD::Double currentGap;

    for (size_t obj = 0; obj < _nobj; ++obj)
    {
        // Sort elements of tmpXInfPInd according to objective obj (in ascending order)
        std::sort(tmpXInfPInd.begin(), tmpXInfPInd.end(),
                  [obj, this](const std::pair<EvalPointPtr, size_t>& t1, const std::pair<EvalPointPtr, size_t>& t2)->bool
                  {
                      return t1.first->getFs(_computeType)[obj] < t2.first->getFs(_computeType)[obj];
                  });

        // Get extreme values value according to one objective
        NOMAD::Double fmin = tmpXInfPInd[0].first->getFs(_computeType)[obj];
        NOMAD::Double fmax = tmpXInfPInd[tmpXInfPInd.size()-1].first->getFs(_computeType)[obj];

        // Can happen for example when we have several minima or for more than three objectives
        if (fmin == fmax)
        {
            fmin = 0.0;
            fmax = 1.0;
        }

        // Intermediate points
        for (size_t i = 1; i < tmpXInfPInd.size()-1;++i)
        {
            if (!canBeFrameCenter[tmpXInfPInd[i].second])
                continue;

            currentGap = tmpXInfPInd[i+1].first->getFs(_computeType)[obj] -
                tmpXInfPInd[i-1].first->getFs(_computeType)[obj];
            currentGap /= (fmax -fmin);
            if (currentGap >= maxGap)
            {
                maxGap = currentGap;
                currentBestInd = tmpXInfPInd[i].second;
            }
        }

        // Extreme points
        currentGap = 2 * (tmpXInfPInd[tmpXInfPInd.size() - 1].first->getFs(_computeType)[obj] -
                          tmpXInfPInd[tmpXInfPInd.size() - 2].first->getFs(_computeType)[obj]);
        currentGap /= (fmax - fmin);
        if (canBeFrameCenter[tmpXInfPInd[tmpXInfPInd.size()-1].second] && currentGap > maxGap)
        {
            maxGap = currentGap;
            currentBestInd = tmpXInfPInd[tmpXInfPInd.size()-1].second;
        }

        currentGap = 2 * (tmpXInfPInd[1].first->getFs(_computeType)[obj] -
                          tmpXInfPInd[0].first->getFs(_computeType)[obj]);
        currentGap /= (fmax -fmin);
        if (canBeFrameCenter[tmpXInfPInd[0].second] && currentGap > maxGap)
        {
            maxGap = currentGap;
            currentBestInd = tmpXInfPInd[0].second;
        }
    }
    return _xInf[currentBestInd];
}


NOMAD::EvalPointPtr NOMAD::DMultiMadsBarrier::getXInfMinH() const
{
    size_t indXInfMinH = 0;
    NOMAD::Double hMinVal = NOMAD::INF;

    for (size_t i = 0; i < _xInf.size(); ++i)
    {
        const NOMAD::Double h = _xInf[i]->getH(_computeType);

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



void NOMAD::DMultiMadsBarrier::clearXInf()
{
    _xInf.clear();
    _xFilterInf.clear();
    
    // Update the current infeasible incumbents.
    updateCurrentIncumbentInfMaxH();
    updateCurrentIncumbentInf();
    updateCurrentIdealInf();
}


void NOMAD::DMultiMadsBarrier::updateRefBests()
{
    throw NOMAD::Exception(__FILE__, __LINE__, "Ref Bests are not used with DMultiMads.");
}


// The code is the very similar to what is in Barrier. Here we use current incumbents.
NOMAD::SuccessType NOMAD::DMultiMadsBarrier::getSuccessTypeOfPoints(const EvalPointPtr xFeas,
                                                                    const EvalPointPtr xInf)
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
        NOMAD::ComputeSuccessType computeSuccess(_computeType);

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


NOMAD::CompareType NOMAD::DMultiMadsBarrier::updateFeasWithPoint(const EvalPoint & evalPoint,
                                                                 const bool keepAllPoints)
{
    // Dealing with Non Standard computing order
    const auto computeTypeS = _computeType.Short();
    if (computeTypeS.computeType != ComputeType::STANDARD)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Update using non standard compute type is not yet implemented");
    }
    
    const auto eval = evalPoint.getEval(_computeType.evalType);
    if (!eval->isFeasible(computeTypeS))
        return NOMAD::CompareType::UNDEFINED;

    std::string s; // for output info
        
    OUTPUT_DEBUG_START
    s = "Point suggested to update DMultiMadsBarrier (feasible): " + evalPoint.display();
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END
        
    if (eval->getFs(computeTypeS).size() != _nobj)
    {
        s = "DMultiMadsBarrier update: number of objectives is equal to " + std::to_string(_nobj);
        s += ". Trying to add this point with number of objectives " + std::to_string(eval->getFs(computeTypeS).size());
        s += ": " + evalPoint.display();
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    if (_xFeas.empty())
    {
        _xFeas.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
        OUTPUT_DEBUG_START
        s = "New dominating xFeas: " + evalPoint.display();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        return NOMAD::CompareType::DOMINATING;
    }

    // Ensure evalPoint is as good as previous points in xFeas
    std::vector<bool> keepInXFeas(_xFeas.size(), true);
    int currentInd = 0;
    auto compFlag = NOMAD::CompareType::INDIFFERENT;
    for (const auto& xFeas: _xFeas)
    {
        const auto currentCompFlag = evalPoint.compMO(*xFeas, _computeType);
        if (currentCompFlag == CompareType::DOMINATED)
        {
            OUTPUT_DEBUG_START
            s = "evalPoint is dominated by " + xFeas->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            s = "evalPoint is rejected";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            return currentCompFlag;
        }
        if (currentCompFlag == CompareType::EQUAL)
        {
            if (!keepAllPoints)
                return currentCompFlag;

            // If new point is not already there, add it
            // One should never add two times the same point into the barrier.
            if (findEvalPoint(_xFeas.begin(), _xFeas.end(), evalPoint) != _xFeas.end())
            {
                OUTPUT_DEBUG_START
                s = "EvalPoint " + evalPoint.display() + "is already in xFeas";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
                return NOMAD::CompareType::UNDEFINED;
            }
            compFlag = NOMAD::CompareType::EQUAL;
        }
        if (currentCompFlag == CompareType::DOMINATING)
        {
            OUTPUT_DEBUG_START
            s = "EvalPoint " + xFeas->display() + "in xFeas is dominated by " + evalPoint.display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            keepInXFeas[currentInd] = false;
            compFlag = NOMAD::CompareType::DOMINATING;
        }
        currentInd++;
    }

    // Remove all dominated elements.
    currentInd = 0;
    const size_t prevNbFeasElements = _xFeas.size();
    _xFeas.erase(std::remove_if(_xFeas.begin(), _xFeas.end(),
                                [&currentInd, &keepInXFeas](const EvalPointPtr& evalPoint)
                                {
                                const bool isRemoved = !keepInXFeas[currentInd];
                                ++currentInd;
                                return isRemoved;
                                }), _xFeas.end());
    OUTPUT_DEBUG_START
    s = "Removing " + std::to_string(std::max((int) prevNbFeasElements - (int)_xFeas.size(), 0));
    s += " dominating feasible eval points";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END
                    
    // Update mesh and insert the new element.
    OUTPUT_DEBUG_START
    s = "Adding new non dominated eval point in xFeas: " + evalPoint.display();
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    const bool updateMesh = [&]() -> bool
    {
        if (compFlag == NOMAD::CompareType::DOMINATING)
            return true;

        // Check if evalPoint extends the current Pareto front approximation.
        bool extends = false;
        for (size_t obj = 0; obj < _nobj; ++obj)
        {
            if (_currentIdealFeas[obj] > evalPoint.getFs(_computeType)[obj])
            {
                extends = true;
                _currentIdealFeas[obj] = evalPoint.getFs(_computeType)[obj];
            }
        }
        return extends;
    }(); // IIFE

    // Update the mesh -> success
    if (updateMesh)
    {
        const auto dir = evalPoint.getDirection();
        if (nullptr != dir)
        {
            OUTPUT_DEBUG_START
            s = "Update the mesh of eval point to be added";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END

            evalPoint.getMesh()->enlargeDeltaFrameSize(*dir);
        }
    }
                    
    // Insert new element
    _xFeas.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
                    
    return compFlag;
}

NOMAD::CompareType NOMAD::DMultiMadsBarrier::updateInfWithPoint(const EvalPoint & evalPoint,
                                                                const bool keepAllPoints)
{
    const auto eval = evalPoint.getEval(_computeType.evalType);
    const auto computeTypeS = _computeType.Short();
    if (eval->isFeasible(computeTypeS))
        return NOMAD::CompareType::UNDEFINED;

    std::string s; // for output info
    
    const NOMAD::Double h = eval->getH(computeTypeS);
    if (!h.isDefined())
    {
        OUTPUT_DEBUG_START
        s = "H is undefined";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        return NOMAD::CompareType::UNDEFINED;
    }
    if (h == NOMAD::INF || ((_hMax < NOMAD::INF) && (h > _hMax)))
    {
        OUTPUT_DEBUG_START
        s = "H is too large: ";
        s += h.display(NOMAD::DISPLAY_PRECISION_FULL) + " > " + _hMax.display(NOMAD::DISPLAY_PRECISION_FULL) + ", continue.";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        return NOMAD::CompareType::UNDEFINED;
    }

    if (_xInf.empty())
    {
        // New point is first point
        _xInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
        _xFilterInf.push_back(std::make_shared<NOMAD::EvalPoint>(evalPoint));
        OUTPUT_DEBUG_START
        s = "New current incumbent infeasible: " + evalPoint.display();
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
        return NOMAD::CompareType::DOMINATING;
    }

    // Insertion into the two sets of infeasible non dominated points.
    // 1- Try to insert into _xInfFilter.
    std::vector<bool> isInXinfFilter(_xFilterInf.size(), true);
    int currentInd = 0;
    for (const auto& xFilterInf: _xFilterInf)
    {
        const auto currentCompFlag = evalPoint.compMO(*xFilterInf, _computeType);
        if (currentCompFlag == CompareType::DOMINATED)
        {
            OUTPUT_DEBUG_START
            s = "evalPoint is dominated by " + xFilterInf->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            s = "evalPoint is rejected";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            return currentCompFlag;
        }
        if (currentCompFlag == CompareType::DOMINATING)
        {
            OUTPUT_DEBUG_START
            s = "xInf dominates the filter element " + xFilterInf->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            isInXinfFilter[currentInd] = false;
        }
        if (currentCompFlag == CompareType::EQUAL)
        {
            if (!keepAllPoints)
                return NOMAD::CompareType::EQUAL;

            // If new point is not already there, add it.
            // One should never insert several times the same point into the barrier.
            if (findEvalPoint(_xFilterInf.begin(), _xFilterInf.end(), evalPoint) != _xFilterInf.end())
            {
                OUTPUT_DEBUG_START
                s = "xInf is already here: " + evalPoint.display();
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
                OUTPUT_DEBUG_END
                return NOMAD::CompareType::UNDEFINED;
            }
        }
        currentInd++;
    }
            
    // Remove all dominated elements of _xFilterInf
    currentInd = 0;
    const size_t prevNbFilterInfElements = _xFilterInf.size();
    _xFilterInf.erase(std::remove_if(_xFilterInf.begin(), _xFilterInf.end(),
                                     [&currentInd, &isInXinfFilter](const EvalPointPtr& ev)
                                     {
                                         bool isRemoved = !isInXinfFilter[currentInd];
                                         currentInd++;
                                         return isRemoved;
                                     }), _xFilterInf.end());
    OUTPUT_DEBUG_START
    s = "Removing " + std::to_string(std::max((int) prevNbFilterInfElements - (int)_xFilterInf.size(), 0));
    s += " dominating filter infeasible eval points";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    OUTPUT_DEBUG_START
    s = "Adding new non dominated eval point in xFilterInf: " + evalPoint.display();
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END
    auto evalPointPtr = std::make_shared<NOMAD::EvalPoint>(evalPoint);
    _xFilterInf.push_back(evalPointPtr);

    // 2- Try to insert the point into _xInf.
    currentInd = 0;
    std::vector<bool> isInXinf(_xInf.size(), true);
    auto compFlag = NOMAD::CompareType::INDIFFERENT;
    for (const auto& xInf: _xInf)
    {
        const auto currentCompFlag = evalPoint.compMO(*xInf, _computeType, true);
        if (currentCompFlag == CompareType::DOMINATED)
        {

            OUTPUT_DEBUG_START
            s = "evalPoint is dominated by " + xInf->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            s = "evalPoint is rejected";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            return currentCompFlag;
        }
        // The points could be only equal according to the f values
        if ((currentCompFlag == CompareType::DOMINATING) ||
            (evalPoint.compMO(*xInf, _computeType) == CompareType::DOMINATING))
        {
            OUTPUT_DEBUG_START
            s = "xInf dominates the filter element " + xInf->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            isInXinf[currentInd] = false;
            compFlag = NOMAD::CompareType::DOMINATING;
        }
        currentInd++;
    }

    // Remove all dominated elements of _xInf.
    currentInd = 0;
    const size_t prevNbInfElements = _xInf.size();
    _xInf.erase(std::remove_if(_xInf.begin(), _xInf.end(),
                               [&currentInd, &isInXinf](const EvalPointPtr& ev)
                               {
                                   bool isRemoved = !isInXinf[currentInd];
                                   currentInd++;
                                   return isRemoved;
                               }), _xInf.end());
    OUTPUT_DEBUG_START
    s = "Removing " + std::to_string(std::max((int) prevNbInfElements - (int)_xInf.size(), 0));
    s += " non dominated infeasible eval points";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    OUTPUT_DEBUG_START
    s = "Adding new non dominated eval point in xInf: " + evalPoint.display();
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    // Update mesh and insert evalPoint into the set of infeasible solutions.
    const bool updateMesh = [&]() -> bool
    {
        if (compFlag == NOMAD::CompareType::DOMINATING)
            return true;

        // Check if evalPoint extends the current Pareto front approximation.
        bool extends = false;
        for (size_t obj = 0; obj < _nobj; ++obj)
        {
            if (_currentIdealInf[obj] > evalPoint.getFs(_computeType)[obj])
            {
                extends = true;
                _currentIdealInf[obj] = evalPoint.getFs(_computeType)[obj];
            }
        }
        return extends;
    }(); // IIFE

    // Update the mesh -> success
    if (updateMesh)
    {
        const auto dir = evalPoint.getDirection();
        if (nullptr != dir)
        {
            OUTPUT_DEBUG_START
            s = "Update the mesh of eval point to be added";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END

            evalPointPtr->getMesh()->enlargeDeltaFrameSize(*dir);
        }
    }

    // Insert new element.
    // NB: It is important that evalPoint be the same into xFilterInf and
    // xInf; most specifically, it should have the same mesh in both sets.
    _xInf.push_back(evalPointPtr);
                    
    return compFlag;
}


// Update the barrier (feas and inf). Once done calls for update the current best feas and inf.
bool NOMAD::DMultiMadsBarrier::updateWithPoints(
                                        const std::vector<EvalPoint>& evalPointList,
                                        const bool keepAllPoints,
                                        const bool updateInfeasibleIncumbentsAndHmax)
{
    bool updatedFeas = false;
    bool updatedInf = false;
    bool updatedIncFeas = false;
    bool updatedIncInf = false;
    bool rejectInf = false;

    const auto evalType = _computeType.evalType;
    const auto computeTypeS = _computeType.Short();

    std::string s; // for output info

    OUTPUT_DEBUG_START
    s = "Updating DMultiMadsBarrier (" + std::to_string(_nobj) + " objectives)";
    s += " with " + std::to_string(evalPointList.size()) + " suggested points";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    s = "Current barrier: ";
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    std::vector<std::string> vs = display(4, false);
    for (const auto& elt: vs)
    {
        NOMAD::OutputQueue::Add(elt, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    OUTPUT_DEBUG_END

    // Order lexicographically a set of non dominated points.
    auto lexicographicallyObjOrder = [](const NOMAD::FHComputeType& computeType,
                                        std::vector<EvalPointPtr>& xSet)
    {
        std::sort(xSet.begin(), xSet.end(),
                  [computeType](const EvalPointPtr& evalPoint1, const EvalPointPtr& evalPoint2)->bool
                  {
                      return evalPoint1->getFs(computeType).lexicographicalCmp(evalPoint2->getFs(computeType));
                  });
    };

    // do separate loop on evalPointList
    // First loop update the bestFeasible.
    // The flag below is set to update the barrier threshold in a second pass.
    NOMAD::SuccessType feasSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
    size_t nbFeasiblePts = 0;
    for (const auto & evalPoint : evalPointList)
    {
        // All points must have mesh parameters defined.
        checkMeshParameters(evalPoint);

        const auto eval = evalPoint.getEval(evalType);
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
        
        if (computeTypeS.computeType == ComputeType::STANDARD && eval->getFs(computeTypeS).size() != _nobj)
        {
            s = "DMultiMadsBarrier update: number of objectives is equal to " + std::to_string(_nobj);
            s += ". Trying to add this point with number of objectives " + std::to_string(eval->getFs(computeTypeS).size());
            s += ": " + evalPoint.display();
            throw NOMAD::Exception(__FILE__, __LINE__, s);
        }

        if (!eval->isFeasible(computeTypeS))
            continue;

        nbFeasiblePts += 1;

        // Add feasible point into the barrier.
        const auto FkCompFlag = updateFeasWithPoint(evalPoint, keepAllPoints);
        const bool insertEqual = (FkCompFlag == NOMAD::CompareType::EQUAL) && !keepAllPoints;
        updatedFeas = ((FkCompFlag != NOMAD::CompareType::UNDEFINED) && !insertEqual) || updatedFeas;

        updatedIncFeas = (FkCompFlag == NOMAD::CompareType::INDIFFERENT) ||
                         (FkCompFlag == NOMAD::CompareType::DOMINATING) ||
                         ((FkCompFlag == NOMAD::CompareType::EQUAL) && keepAllPoints) ||
                         updatedIncFeas;

        // Set the success type according to the current infeasible incumbent.
        if (feasSuccessType == NOMAD::SuccessType::FULL_SUCCESS)
            continue;

        // If the set of feasible incumbents is empty, we follow the same approach as for the
        // progressive barrier. The iteration is then considered as a success.
        if (_currentIncumbentFeas == nullptr)
        {
            feasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
            continue;
        }

        const auto compFlag = evalPoint.compMO(*_currentIncumbentFeas, _computeType);
        if (compFlag == NOMAD::CompareType::DOMINATING)
            feasSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
    }

    // Lexicographically order the set of feasible non-dominated solutions.
    if (updatedIncFeas)
        lexicographicallyObjOrder(_computeType, _xFeas);

    // Do separate loop on evalPointList
    // Second loop update the bestInfeasible.
    // The flag below is set to update the barrier threshold in a second pass.
    NOMAD::SuccessType infSuccessType = NOMAD::SuccessType::UNSUCCESSFUL;
    const bool areAllFeasible = nbFeasiblePts == evalPointList.size();
    for (const auto & evalPoint : evalPointList)
    {
        // No need to continue.
        if (areAllFeasible)
            break;

        const auto eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            // Suggested point is not good.
            continue;
        }

        if (eval->isFeasible(computeTypeS))
            continue;

        const NOMAD::Double h = eval->getH(computeTypeS);
        if (!h.isDefined() || (h == NOMAD::INF) ||
            ((_hMax < NOMAD::INF) && (h > _hMax)))
        {
            OUTPUT_DEBUG_START
            s = !h.isDefined() ? "H is undefined" : "H is too large: ";
            if (h.isDefined())
            {
                s += h.display(NOMAD::DISPLAY_PRECISION_FULL) + " > " + _hMax.display(NOMAD::DISPLAY_PRECISION_FULL) + ", continue.";
            }
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            rejectInf = true;
            continue;
        }

        // Add infeasible point into the barrier.
        const auto IkCompFlag = updateInfWithPoint(evalPoint, keepAllPoints);
        const bool insertEqual = (IkCompFlag == NOMAD::CompareType::EQUAL) && !keepAllPoints;
        updatedInf = ((IkCompFlag != NOMAD::CompareType::UNDEFINED) && !insertEqual) || updatedInf;

        updatedIncInf = (IkCompFlag == NOMAD::CompareType::INDIFFERENT) ||
                        (IkCompFlag == NOMAD::CompareType::DOMINATING) ||
                        ((IkCompFlag == NOMAD::CompareType::EQUAL) && keepAllPoints) ||
                        updatedIncInf;

        // Set the success type according to the current infeasible incumbent.
        if (feasSuccessType == NOMAD::SuccessType::FULL_SUCCESS ||
            infSuccessType == NOMAD::SuccessType::FULL_SUCCESS)
            continue;

        // If the set of infeasible incumbents is empty, we follow the same approach as for the
        // progressive barrier. The iteration is then considered as a success.
        if (_currentIncumbentInf == nullptr)
        {
            infSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
            continue;
        }

        const auto compFlag = evalPoint.compMO(*_currentIncumbentInf, _computeType);
        if (compFlag == NOMAD::CompareType::DOMINATING)
        {
            infSuccessType = NOMAD::SuccessType::FULL_SUCCESS;
            continue;
        }

        const NOMAD::Double hXInf = _currentIncumbentInf->getH(_computeType);
        if (h.isDefined() && h < hXInf)
            infSuccessType = NOMAD::SuccessType::PARTIAL_SUCCESS;
    }

    // Lexicographically order the set of filter points and infeasible points.
    if (updatedInf)
    {
        lexicographicallyObjOrder(_computeType, _xFilterInf);
        lexicographicallyObjOrder(_computeType, _xInf);
    }

    const bool incumbentsAndHMaxUpToDate = ! (updatedFeas || updatedInf);

    // Update hMax: when some infeasible points are rejected, the iteration is considered
    // as a failure. We still update hMax.
    NOMAD::Double hMaxPrev = _hMax;
    NOMAD::Double hMax = _hMax;
    if (updateInfeasibleIncumbentsAndHmax && (!incumbentsAndHMaxUpToDate || rejectInf))
    {
        if (infSuccessType == NOMAD::SuccessType::PARTIAL_SUCCESS &&
            feasSuccessType < NOMAD::SuccessType::FULL_SUCCESS)
        {
            // hMax <- max_{x \in Uk+1} {h(x) : h(x) < h(xinf)}
            // Note : An improving point (with respect to the infeasible incumbent)
            // must have been generated.
            hMax = 0.0;
            const NOMAD::Double hCurrentXInf = _currentIncumbentInf->getH(_computeType);
            for (const auto& xFilterInf: _xFilterInf)
            {
                const auto eval = xFilterInf->getEval(_computeType.evalType);
                const NOMAD::Double h = eval->getH(_computeType.Short());
                if (h < hCurrentXInf)
                    hMax = std::max(hMax, h);
            }

            OUTPUT_DEBUG_START
            s = "Partial success";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
        }
        else
        {
            OUTPUT_DEBUG_START
            s = ((feasSuccessType == NOMAD::SuccessType::FULL_SUCCESS) ||
                 (infSuccessType == NOMAD::SuccessType::FULL_SUCCESS)) ? "Full success" : "no success";
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END

            if (!_xInf.empty())
            {
                // hMax <- max_{x \in Uk+1} {h(x) : h(xinf) <= h(x) < max_{x in Ik} h(x)}
                // When _currentIncumbentInfHmaxH does not exist, set hMax to
                // hMax := max_{x in Ik+1} h(x).
                const NOMAD::Double hMaxLimSup =
                    _currentIncumbentInfMaxH == nullptr ? hMaxPrev
                                                        : _currentIncumbentInfMaxH->getH(_computeType);
                hMax = _currentIncumbentInf == nullptr ? 0.0 : _currentIncumbentInf->getH(_computeType);
                for (const auto& xFilterInf: _xFilterInf)
                {
                const auto eval = xFilterInf->getEval(_computeType.evalType);
                const NOMAD::Double h = eval->getH(_computeType.Short());
                if (h < hMaxLimSup)
                    hMax = std::max(hMax, h);
                }
            }
        }
        // If there is a problem, reset it to hMaxPrev;
        if (hMax == 0)
            hMax = hMaxPrev;

        // Set hMax and remove infeasible points above the threshold.
        if (hMax != hMaxPrev)
        {
            setHMax(hMax);
        }
    }

    // Set n and check that all points have the same dimension
    const bool updatedInc = updatedIncFeas || updatedIncInf;
    if (updatedInc)
    {
        setN();
    }

    const bool updated = updatedFeas || updatedInf;
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


void NOMAD::DMultiMadsBarrier::updateCurrentIncumbentInfMaxH()
{
    _currentIncumbentInfMaxH = nullptr;
    NOMAD::Double currentHMax = 0.0;
    for (const auto& xInf: _xInf)
    {
        const NOMAD::Double h = xInf->getH(_computeType);
        if (h.isDefined() && h > currentHMax)
        {
            currentHMax = h;
            _currentIncumbentInfMaxH = xInf;
        }
    }
}


void NOMAD::DMultiMadsBarrier::updateXInfAndFilterInfAfterHMaxSet()
{
    if (_xInf.empty())
    {
        return;
    }

    size_t currentInd = 0;
    
    // Remove all infeasible incumbent solutions above the threshold _hMax.
    std::vector<bool> isInXInf(_xInf.size(), true);
    for (const auto& xInf: _xInf)
    {
        const NOMAD::Double h = xInf->getH(_computeType);

        if (h > _hMax)
        {
            isInXInf[currentInd] = false;

        }
        currentInd++;
    }

    currentInd = 0;
    _xInf.erase(std::remove_if(_xInf.begin(), _xInf.end(),
                               [&currentInd, &isInXInf](const EvalPointPtr& evalPoint)
                               {
                                   bool isRemoved = !isInXInf[currentInd];
                                   ++currentInd;
                                   return isRemoved;
                               }), _xInf.end());

    // Remove all infeasible non dominated solutions above the threshold _hMax.
    currentInd = 0;
    std::vector<bool> isInXFilterInf(_xFilterInf.size(), true);
    for (const auto& xFilterInf: _xFilterInf)
    {
        const NOMAD::Double h = xFilterInf->getH(_computeType);

        if (h > _hMax)
        {
            isInXFilterInf[currentInd] = false;

        }
        currentInd++;
    }

    currentInd = 0;
    _xFilterInf.erase(std::remove_if(_xFilterInf.begin(), _xFilterInf.end(),
                                     [&currentInd, &isInXFilterInf](const EvalPointPtr& evalPoint)
                                     {
                                         bool isRemoved = !isInXFilterInf[currentInd];
                                         ++currentInd;
                                         return isRemoved;
                                     }), _xFilterInf.end());

    std::sort(_xFilterInf.begin(), _xFilterInf.end(),
              [this](const EvalPointPtr& evalPoint1, const EvalPointPtr& evalPoint2) -> bool
              {
                  return evalPoint1->getFs(_computeType).lexicographicalCmp(evalPoint2->getFs(_computeType));
              });

    // And reinsert potential infeasible non dominated points into the set of infeasible
    // solutions.
    currentInd = 0;
    isInXInf = std::vector<bool>(_xFilterInf.size(), false);
    for (const auto& evalPoint: _xFilterInf)
    {
        if (findEvalPoint(_xInf.begin(), _xInf.end(), *evalPoint) != _xInf.end())
        {
            currentInd++;
            continue;
        }

        // If new point is not in _xInf, check if it is not dominated
        size_t currentIndTmp = 0;
        bool insert = true;
        for (const auto& evalPointInf: _xFilterInf)
        {
            if (currentIndTmp == currentInd)
            {
                currentIndTmp++;
                continue;
            }

            const auto compFlag = evalPoint->compMO(*evalPointInf, _computeType, true);
            if (compFlag == CompareType::DOMINATED)
            {
                insert = false;
                break;
            }
            if (compFlag == CompareType::DOMINATING)
            {
                isInXInf[currentIndTmp] = false;
            }
            currentIndTmp++;
        }
        isInXInf[currentInd] = insert;
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
              [this](const EvalPointPtr& evalPoint1, const EvalPointPtr& evalPoint2) -> bool
              {
                  return evalPoint1->getFs(_computeType).lexicographicalCmp(evalPoint2->getFs(_computeType));
              });
}


// Nice formatting display.
std::vector<std::string> NOMAD::DMultiMadsBarrier::display(const size_t max, const bool displayMeshes) const
{
    std::vector<std::string> vs;

    size_t nbXFeas = 0;
    size_t nbXInf = 0;
    size_t nbXFilterInf = 0;

    for (const auto & xFeas : _xFeas)
    {
        vs.push_back("X_FEAS " + xFeas->displayAll(NOMAD::defaultFHComputeTypeS));
        if (displayMeshes)
        {
            const auto mesh = xFeas->getMesh();
            vs.push_back("delta mesh size = " + mesh->getdeltaMeshSize().display());
            vs.push_back("Delta mesh size = " + mesh->getDeltaFrameSize().display());
        }
        nbXFeas++;
        if (nbXFeas >= max && _xFeas.size() > max)
        {
            vs.push_back("... (total " + std::to_string(_xFeas.size()) + ")");
            break;
        }
    }
    for (const auto &xInf : _xInf)
    {
        vs.push_back("X_INF " + xInf->displayAll(NOMAD::defaultFHComputeTypeS) );
        if (displayMeshes)
        {
            const auto mesh = xInf->getMesh();
            vs.push_back("delta mesh size = " + mesh->getdeltaMeshSize().display());
            vs.push_back("Delta mesh size = " + mesh->getDeltaFrameSize().display());
        }
        nbXInf++;

        if (nbXInf >= max && _xInf.size() > max)
        {
            vs.push_back("... (total " + std::to_string(_xInf.size()) + ")");
            break;
        }
    }
    for (const auto &xFilterInf: _xFilterInf)
    {
        vs.push_back("X_FILTER" + xFilterInf->displayAll(NOMAD::defaultFHComputeTypeS));
        if (displayMeshes)
        {
            const auto mesh = xFilterInf->getMesh();
            vs.push_back("delta mesh size = " + mesh->getdeltaMeshSize().display());
            vs.push_back("Delta mesh size = " + mesh->getDeltaFrameSize().display());
        }
        nbXFilterInf++;

        if (nbXFilterInf >= max && _xFilterInf.size() > max)
        {
            vs.push_back("... (total " + std::to_string(_xFilterInf.size()) + ")");
            break;
        }
    }

    vs.push_back("H_MAX " + getHMax().display(NOMAD::DISPLAY_PRECISION_FULL));
    vs.push_back("Ref Best Feasible:   " + (_refBestFeas ? _refBestFeas->displayAll(NOMAD::defaultFHComputeTypeS) : "NULL") );
    vs.push_back("Ref Best Infeasible: " + (_refBestInf ? _refBestInf->displayAll(NOMAD::defaultFHComputeTypeS) : "NULL"));

    return vs;
}



NOMAD::Double NOMAD::DMultiMadsBarrier::getMeshMaxFrameSize(const NOMAD::EvalPointPtr& pt) const
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
    if (_xFeas.empty() || _xInf.empty())
    {
        _currentIncumbentInf = getFirstXIncInfNoXFeas();
        return;
    }

    const auto evalType = _computeType.evalType;
    const auto computeTypeS = _computeType.fhComputeTypeS;

    auto computeDomMove = [evalType, computeTypeS](const size_t nobj,
                                                   const EvalPointPtr& xInf,
                                                   const std::vector<EvalPointPtr>& xFeasElems,
                                                   const bool domXFeasElems) -> double
    {
        double minDomMove = std::numeric_limits<double>::infinity();
        for (const auto& xFeas: xFeasElems)
        {
            const NOMAD::Eval* evalFeas = xFeas->getEval(evalType);
            const NOMAD::Eval* evalInf = xInf->getEval(evalType);
            double domMove = 0.0;
            for (size_t i = 0; i < nobj; ++i)
            {
                const double fFeasVal = evalFeas->getFs(computeTypeS)[i].todouble();
                const double fInfVal = evalInf->getFs(computeTypeS)[i].todouble();
                const double diffVal = domXFeasElems ? fFeasVal - fInfVal
                                                     : fInfVal - fFeasVal;
                domMove += std::max(diffVal, 0.0);
            }
            minDomMove = std::min(minDomMove, domMove);
        }
        return minDomMove;
    };

    // Get the infeasible solution with maximum dominance move below the _hMax threshold,
    // according to the set of best feasible incumbent solutions.
    size_t currentInd = 0;
    double maxDomMove = -NOMAD::INF;
    for (size_t j = 0; j < _xInf.size(); ++j)
    {
        const NOMAD::Eval* evalInf = _xInf[j]->getEval(evalType);
        const NOMAD::Double h = evalInf->getH(computeTypeS);
        if (!h.isDefined() || h > _hMax)
            continue;

        // Compute dominance move
        // = min \sum_{1}^m max(fi(y) - fi(x), 0)
        //   y \in Fk
        const double domMove = computeDomMove(_nobj, _xInf[j], _xFeas, true /*domXFeasElems*/);

        // Get the maximum dominance move index
        if (maxDomMove < domMove) {
            maxDomMove = domMove;
            currentInd = j;
        }
    }
    if (NOMAD::Double(maxDomMove) > 0.0)
    {
        _currentIncumbentInf = _xInf[currentInd];
        return;
    }

    // In this case, all infeasible solutions are "dominated" in terms of fvalues
    // by at least one element of Fk. Get the infeasible solution below the _hMax
    // threshold which has minimal dominance move, when considered a maximization
    // problem.
    double minDomMove = NOMAD::INF;
    currentInd = 0;
    for (size_t j = 0; j < _xInf.size(); ++j)
    {
        const NOMAD::Eval* evalInf = _xInf[j]->getEval(evalType);
        const NOMAD::Double h = evalInf->getH(computeTypeS);
        if (!h.isDefined() || h > _hMax)
            continue;

        // Compute dominance move
        // = min \sum_{1}^m max(fi(x) - fi(y), 0)
        //   y \in Fk
        const double domMove = computeDomMove(_nobj, _xInf[j], _xFeas, false /*domXFeasElems*/);

        // Get the minimal dominance move index
        if (minDomMove > domMove)
        {
            minDomMove = domMove;
            currentInd = j;
        }
    }
    _currentIncumbentInf = _xInf[currentInd];
}


void NOMAD::DMultiMadsBarrier::updateCurrentIncumbentFeas()
{
    if (_xFeas.empty())
    {
        _currentIncumbentFeas = nullptr;
        return;
    }
    
    if (_xFeas.size() == 1)
    {
        _currentIncumbentFeas =_xFeas[0];
        return;
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
    
    // See article DMultiMads Algorithm 4.
    for (size_t i = 0; i < _xFeas.size(); ++i)
    {
        
        const NOMAD::Double maxFrameSizeElt = getMeshMaxFrameSize(_xFeas[i]);
        
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
            const size_t selectedInd = std::distance(canBeFrameCenter.begin(), it);
            _currentIncumbentFeas =_xFeas[selectedInd];
        }
        return;
    }

    // Only two points in the barrier.
    if ((nbSelectedCandidates == 2) && (_xFeas.size() == 2))
    {
        const auto objv1 = _xFeas[0]->getFs(_computeType);
        const auto objv2 = _xFeas[1]->getFs(_computeType);
        _currentIncumbentFeas = (objv1.abs().max() > objv2.abs().max()) ? _xFeas[0] : _xFeas[1];
        return;
    }

    // More than three points in the barrier.
    // First case: biobjective optimization. Points are already ranked by lexicographic order.
    if (_nobj == 2)
    {
        size_t currentBestInd = 0;
        NOMAD::Double maxGap = -1.0;
        NOMAD::Double currentGap;
        for (size_t obj = 0; obj < _nobj; ++obj)
        {
            // Get extreme values value according to one objective
            const NOMAD::Double fmin = _xFeas[0]->getFs(_computeType)[obj];
            const NOMAD::Double fmax = _xFeas[_xFeas.size()-1]->getFs(_computeType)[obj];
                
            // For biobjective optimization, it means all elements of _xFeas are equal.
            // We return the first one.
            if (fmin == fmax)
            {
                _currentIncumbentFeas =_xFeas[currentBestInd];
                return;
            }
                
            // Intermediate points
            for (size_t i = 1; i < _xFeas.size()-1;++i)
            {
                if (!canBeFrameCenter[i])
                    continue;

                currentGap = _xFeas[i+1]->getFs(_computeType)[obj] -
                    _xFeas[i-1]->getFs(_computeType)[obj];
                currentGap /= (fmax -fmin);
                if (currentGap >= maxGap)
                {
                    maxGap = currentGap;
                    currentBestInd = i;
                }
            }
                
            // Extreme points
            currentGap = 2 * (_xFeas[_xFeas.size() - 1]->getFs(_computeType)[obj] -
                              _xFeas[_xFeas.size() - 2]->getFs(_computeType)[obj]);
            currentGap /= (fmax - fmin);
            if (canBeFrameCenter[_xFeas.size()-1] && currentGap > maxGap)
            {
                maxGap = currentGap;
                currentBestInd = _xFeas.size()-1;
            }
                
            currentGap = 2 * (_xFeas[1]->getFs(_computeType)[obj] -
                              _xFeas[0]->getFs(_computeType)[obj]);
            currentGap /= (fmax -fmin);
            if (canBeFrameCenter[0] && currentGap > maxGap)
            {
                maxGap = currentGap;
                currentBestInd = 0;
            }
        }
        _currentIncumbentFeas =_xFeas[currentBestInd];
        return;
    }

    // Second case: more than 2 objectives
    std::vector<std::pair<EvalPointPtr, size_t> > tmpXFeasPInd(_xFeas.size());
            
    // Initialize it.
    for (size_t i = 0; i < tmpXFeasPInd.size(); ++i)
    {
        tmpXFeasPInd[i] = std::make_pair(_xFeas[i], i);
    }

    size_t currentBestInd = 0;
    NOMAD::Double maxGap = -1.0;
    NOMAD::Double currentGap;
            
    for (size_t obj = 0; obj < _nobj; ++obj)
    {
        // Sort elements of tmpXFeasPInd according to objective obj (in ascending order)
        std::sort(tmpXFeasPInd.begin(), tmpXFeasPInd.end(),
                  [obj,this](const std::pair<EvalPointPtr, size_t>& t1, const std::pair<EvalPointPtr, size_t>& t2)->bool
                  {
                      return t1.first->getFs(_computeType)[obj] < t2.first->getFs(_computeType)[obj];
                  });
                
        // Get extreme values value according to one objective
        NOMAD::Double fmin = tmpXFeasPInd[0].first->getFs(_computeType)[obj];
        NOMAD::Double fmax = tmpXFeasPInd[tmpXFeasPInd.size()-1].first->getFs(_computeType)[obj];
                
        // Can happen for example when we have several minima or for more than three objectives
        if (fmin == fmax)
        {
            fmin = 0.0;
            fmax = 1.0;
        }
                
        // Intermediate points
        for (size_t i = 1; i < tmpXFeasPInd.size()-1;++i)
        {
            if (!canBeFrameCenter[tmpXFeasPInd[i].second])
                continue;

            currentGap = tmpXFeasPInd[i+1].first->getFs(_computeType)[obj] -
                tmpXFeasPInd[i-1].first->getFs(_computeType)[obj];
            currentGap /= (fmax -fmin);
            if (currentGap >= maxGap)
            {
                maxGap = currentGap;
                currentBestInd = tmpXFeasPInd[i].second;
            }
        }
                
        // Extreme points
        currentGap = 2 * (tmpXFeasPInd[tmpXFeasPInd.size() - 1].first->getFs(_computeType)[obj] -
                          tmpXFeasPInd[tmpXFeasPInd.size() - 2].first->getFs(_computeType)[obj]);
        currentGap /= (fmax - fmin);
        if (canBeFrameCenter[tmpXFeasPInd[tmpXFeasPInd.size()-1].second] && currentGap > maxGap)
        {
            maxGap = currentGap;
            currentBestInd = tmpXFeasPInd[tmpXFeasPInd.size()-1].second;
        }
                
        currentGap = 2 * (tmpXFeasPInd[1].first->getFs(_computeType)[obj] -
                          tmpXFeasPInd[0].first->getFs(_computeType)[obj]);
        currentGap /= (fmax -fmin);
        if (canBeFrameCenter[tmpXFeasPInd[0].second] && currentGap > maxGap)
        {
            maxGap = currentGap;
            currentBestInd = tmpXFeasPInd[0].second;
        }
    }
    _currentIncumbentFeas = _xFeas[currentBestInd];
}


void NOMAD::DMultiMadsBarrier::updateCurrentIdealFeas()
{
    _currentIdealFeas = std::vector<NOMAD::Double>(_nobj, NOMAD::INF);
    for (const auto& xFeas: _xFeas)
    {
        for (size_t obj = 0; obj < _nobj; ++obj)
        {
            if (_currentIdealFeas[obj] > xFeas->getFs(_computeType)[obj])
                _currentIdealFeas[obj] = xFeas->getFs(_computeType)[obj];
        }
    }
}


void NOMAD::DMultiMadsBarrier::updateCurrentIdealInf()
{
    _currentIdealInf = std::vector<NOMAD::Double>(_nobj, NOMAD::INF);
    for (const auto& xInf: _xInf)
    {
        const auto& fvaluesXinf = xInf->getFs(_computeType);
        for (size_t obj = 0; obj < _nobj; ++obj)
        {
            if (_currentIdealInf[obj] > fvaluesXinf[obj])
                _currentIdealInf[obj] = xInf->getFs(_computeType)[obj];
        }
    }
}
