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
/**
 \file   CacheSet.cpp
 \brief  Simple implementation of Cache derived from CacheBase (implementation)
 \author Viviane Rochon Montplaisir
 \date   January 2018
 \see    CacheSet.hpp
 */
#include "../Cache/CacheSet.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Util/fileutils.hpp"
#include "../Util/MicroSleep.hpp"
#include "BBOutputType.hpp"
#include "CompareType.hpp"
#include "ComputeType.hpp"
#include "Eval.hpp"
#include "EvalPoint.hpp"

#include <fstream>
#include <iostream>
#include <iterator>

// Init static members
NOMAD::BBOutputTypeList NOMAD::CacheSet::_bbOutputType = NOMAD::BBOutputTypeList();
NOMAD::ArrayOfDouble NOMAD::CacheSet::_bbEvalFormat = NOMAD::ArrayOfDouble();
std::unique_ptr<NOMAD::CacheBase> NOMAD::CacheBase::_single = nullptr;

std::atomic<size_t> NOMAD::CacheBase::_nbCacheHits;

#ifdef _OPENMP
omp_lock_t NOMAD::CacheSet::_cacheLock;
#endif // _OPENMP


// Initialize CacheSet class.
// To be called by the Constructor.
void NOMAD::CacheSet::init()
{
    if (_cacheParams->toBeChecked())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "CacheParameters::checkAndComply() needs to be called before constructing a CacheSet.");
    }
}


// Terminate CacheSet class.
// To be called by the Destructor.
void NOMAD::CacheSet::destroy()
{
    // Clear the _cache directly.
    // No need to set lock, assuming there is only one cache and
    // that now it is the end of the run, and we are calling its destructor.
    _cache.clear();

#ifdef _OPENMP
    omp_destroy_lock(&_cacheLock);
#endif // _OPENMP
}

void NOMAD::CacheSet::setInstance(const std::shared_ptr<NOMAD::CacheParameters>& cacheParams,
                                  const NOMAD::BBOutputTypeList& bbOutputType,
                                  const NOMAD::ArrayOfDouble& bbEvalFormat)
{
#ifdef _OPENMP
    #pragma omp critical(initCacheLock)
    {
#endif // _OPENMP
        if (nullptr == _single)
        {
#ifdef _OPENMP
            omp_init_lock(&_cacheLock);
#endif // _OPENMP
            _single = std::unique_ptr<NOMAD::CacheSet>(new CacheSet(cacheParams)) ;
        }
        else if (_single->size() != 0)
        {
            std::string err = "Cache is not empty while calling NOMAD::CacheSet::setInstance more than ONCE. Need to reset the cache." ;
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }
#ifdef _OPENMP
    }   // end of critical section
#endif // _OPENMP

    _bbOutputType = bbOutputType;
    _bbEvalFormat = bbEvalFormat;

    // As long as the cache file exists, it is read.
    getInstance()->read();
}


void NOMAD::CacheSet::verifyPointComplete(const NOMAD::Point& point) const
{
    if (!point.isComplete())
    {
        std::string err = "Error: Cache does not support incomplete points.";
        err += " Got point: " + point.display();
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::CacheSet::verifyPointSize(const NOMAD::Point& point) const
{
    if (!_cache.empty() && _n != point.size())
    {
        std::string err = "Error: Cache method called with a point of size ";
        err += std::to_string(point.size());
        err += ": " + point.display();
        err += ". Cache needs points of size " + std::to_string(_n);
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::CacheSet::verifyPointComplete(const NOMAD::EvalPoint& evalPoint) const
{
    verifyPointComplete(*(evalPoint.getX()));
}


void NOMAD::CacheSet::verifyPointSize(const NOMAD::EvalPoint& evalPoint) const
{
    verifyPointSize(*(evalPoint.getX()));
}


// Add a new eval point to the cache.
// Return true if insertion worked, false if not (EvalPoint was already there).
// Note: This is not an optimal implementation. This method is not used
// in the nomad code. It is kept for compliance with the CacheBase interface
// and also because unit tests use it extensively.
bool NOMAD::CacheSet::insert(const NOMAD::EvalPoint &evalPoint)
{
    NOMAD::EvalPoint evalPointFound;
    // Return true if point is not found: Expecting it to be inserted.
    bool inserted = (0 == find(evalPoint, evalPointFound));
    // Ignore smartInsert's return value
    smartInsert(evalPoint, NOMAD::INF_SHORT, NOMAD::EvalType::BB);

    return inserted;
}


// Get EvalPoint evalPoint at Point x from the cache
// Returns the number of EvalPoints found
size_t NOMAD::CacheSet::find(const NOMAD::Point& x, NOMAD::EvalPoint &evalPoint,
                             const NOMAD::EvalType evalType,
                             bool waitIfNotYetAvailable ) const
{
    size_t nbFound = 0;

    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    it = _cache.find(NOMAD::EvalPoint(x));
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    if (it != _cache.end())
    {        
#ifdef _OPENMP
        // Wait for evaluation:
        // If using OpenMP, the EvalPoint may be updated by another thread.
        // Otherwise, do not wait.
        if ( waitIfNotYetAvailable && NOMAD::EvalType::UNDEFINED != evalType)
        {
            
            auto evalStatus = it->getEvalStatus(evalType);
            if ( NOMAD::EvalStatusType::EVAL_IN_PROGRESS == evalStatus )
            {
                OUTPUT_INFO_START
                std::string s = "Start waiting for point ";
                s += x.display() + " to complete.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
                OUTPUT_INFO_END
            }
            while (!_stopWaiting
                   && (NOMAD::EvalStatusType::EVAL_IN_PROGRESS == evalStatus
                       || NOMAD::EvalStatusType::EVAL_NOT_STARTED == evalStatus
                       || NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED == evalStatus))
            {
                usleep(10);
                evalStatus = it->getEvalStatus(evalType);
            }
            if (_stopWaiting && NOMAD::EvalStatusType::EVAL_IN_PROGRESS == evalStatus)
            {
                OUTPUT_INFO_START
                std::string s = "Force stop waiting for point ";
                s += x.display() + " to complete.";
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
                OUTPUT_INFO_END
            }
        }
#endif // _OPENMP
        nbFound = 1;
        evalPoint = *it;
    }
    return nbFound;
}

// Get EvalPoint evalPoint at Point x from the cache for rerun.
bool NOMAD::CacheSet::findInCacheForRerun(const NOMAD::Point& x, NOMAD::EvalPoint &evalPoint) const
{

    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    it = _cacheForRerun.find(NOMAD::EvalPoint(x));
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    if (it != _cacheForRerun.end())
    {
        evalPoint = *it;
        return true;
    }
    return false;
}


// Insert evalPoint in cache.
// Return a boolean indicating if we should eval this point.
// If insertion worked, the point was not in the cache before. Return true.
// If insertion did not work, the point was in the cache before.
// Depending on its EvalStatus, return true if it should be evaluated again,
// false otherwise.
bool NOMAD::CacheSet::smartInsert(const NOMAD::EvalPoint &evalPoint,
                                  short maxNumberEval,
                                  NOMAD::EvalType evalType)
{
    verifyPointComplete(evalPoint);
    verifyPointSize(evalPoint);

    // First insert sets n (even if insert fails)
    if (_cache.empty())
    {
        _n = evalPoint.size();
    }

    bool inserted = false;
    std::pair<NOMAD::EvalPointSet::iterator,bool> ret;   // Return of the insert()
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    ret = _cache.insert(evalPoint);
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    inserted = ret.second;
    bool canEval = (*ret.first).toEval(maxNumberEval, evalType);
    bool doEval = canEval;
    
    if (-1 == evalPoint.getTag())
    {
        throw NOMAD::Exception(__FILE__, __LINE__," Eval point should have its tag set before smart insert.");
    }

    if (inserted && canEval)
    {
        doEval = true;
    }
    else if (nullptr == (ret.first)->getEval(evalType))
    {
        // Point already in cache, but not evaluated.
        // NOTE: We do not know if this point is in the evaluation queue yet, or not.
        // We do not know here if the evaluation queue is cleared between runs.
        // If doEval is set to true, the point could be evaluated twice.
        // If doEval is set to false, there might be cases where it is not evaluated at all.

        // Only warn when in blackbox context.
        if (NOMAD::EvalType::BB == evalType)
        {

            // Upate the tag of point already in cache but not evaluated
            ret.first->setTag(evalPoint.getTag());
            
            OUTPUT_INFO_START
            std::string s = "Point already in cache (but not BB evaluated): ";
            s += ret.first->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
            OUTPUT_INFO_END

            // Avoid re-evaluating BB.
            doEval = canEval;
        }
        else if (NOMAD::EvalType::MODEL == evalType)
        {
            // It is ok to re-evaluate MODEL points.
            doEval = true;
        }
        else if (NOMAD::EvalType::SURROGATE == evalType)
        {
            doEval = canEval;
        }
    }
    else
    {
        // Cache hit.
        doEval = canEval;

        // Only count as cache hit when using Blackbox Eval.
        if (!inserted && NOMAD::EvalType::BB == evalType)
        {
            _nbCacheHits++;
            OUTPUT_INFO_START
            std::string s = "Cache hit: ";
            s += ret.first->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
            OUTPUT_INFO_END
        }
        if (doEval)
        {
            std::cout << "Warning: CacheSet: smartInsert: New evaluation of point found in cache " << (*ret.first).display() << std::endl;
        }
    }

    return doEval;
}


size_t NOMAD::CacheSet::find(const NOMAD::Point& x, std::vector<NOMAD::EvalPoint> &evalPointList) const
{
    verifyPointComplete(x);
    verifyPointSize(x);
    evalPointList.clear();

    NOMAD::EvalPoint evalPoint;
    size_t found = find(x, evalPoint);
    if (found > 0)
    {
        evalPointList.push_back(evalPoint);
    }
    return found;
}


size_t NOMAD::CacheSet::find(const NOMAD::Eval &refeval,
                             std::function<bool(const NOMAD::Eval&, const NOMAD::Eval&, const NOMAD::FHComputeTypeS&)> comp,
                             std::vector<NOMAD::EvalPoint> &evalPointList,
                             const NOMAD::FHComputeType& computeType) const
{
    evalPointList.clear();
    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        const NOMAD::Eval* eval = it->getEval(computeType.evalType);
        if (nullptr == eval)
        {
            continue;
        }
        if (comp(*eval, refeval, computeType.Short()))
        {
            const NOMAD::EvalPoint& evalPoint(*it);
            evalPointList.push_back(evalPoint);
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP

    return evalPointList.size();
}


// Get best eval points, using comp()
size_t NOMAD::CacheSet::findBest(std::function<bool(const NOMAD::Eval&, const NOMAD::Eval&, const NOMAD::FHComputeTypeS&)> comp,
                     std::vector<NOMAD::EvalPoint> &evalPointList,
                     const bool findFeas,
                     const NOMAD::Double& hMax,
                     const NOMAD::Point& fixedVariable,
                     const NOMAD::FHComputeType& computeType) const
{
    evalPointList.clear();
    NOMAD::EvalPointSet::const_iterator it;
    NOMAD::Eval refeval;
    
    auto evalType = computeType.evalType;
    auto compactComputeType = computeType.Short();

#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        const NOMAD::EvalPoint& evalPoint(*it);
        const NOMAD::Eval* eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (findFeas != eval->isFeasible(compactComputeType))
        {
            continue;
        }
        NOMAD::Double h = eval->getH(compactComputeType);
        if (! h.isDefined())
        {
            continue;
        }
        // If hMax == INF all infeasible points (PB and EB) are considered. Otherwise, only h <=hMax are considered
        if ( hMax < NOMAD::INF && h > hMax )
        {
            continue;
        }
        // Must be in the subspace defined by fixedVariable
        if (!evalPoint.hasFixed(fixedVariable))
        {
            continue;
        }

        if (refeval.getEvalStatus()==NOMAD::EvalStatusType::EVAL_STATUS_UNDEFINED)
        {
            // Found first point
            refeval = *eval;
            evalPointList.push_back(evalPoint);
        }
        else if (*eval == refeval)
        {
            // Found first point
            // Found a point with eval == refeval
            evalPointList.push_back(evalPoint);
        }
        else if (comp(*eval, refeval, compactComputeType))
        {
            // Found a better point
            refeval = *eval;
            // Reset list with new best
            evalPointList.clear();
            evalPointList.push_back(evalPoint);
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP

    return evalPointList.size();
}

bool NOMAD::CacheSet::hasFeas(const NOMAD::FHComputeType& computeType) const
{
    bool ret = false;

#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (const auto& it : _cache)
    {
        const NOMAD::Eval* eval = it.getEval(computeType.evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (eval->isFeasible(computeType.Short()))
        {
            ret = true;
            break;
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP

    return ret;
}


bool NOMAD::CacheSet::hasInfeas(const NOMAD::FHComputeType& computeType) const
{
    bool ret = false;

#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (const auto& it : _cache)
    {
        const NOMAD::Eval* eval = it.getEval(computeType.evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (!eval->isFeasible(computeType.Short()))
        {
            ret = true;
            break;
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP

    return ret;
}



size_t NOMAD::CacheSet::find(const NOMAD::Point & X,
                             std::function<bool(const NOMAD::Point&, const NOMAD::EvalPoint &)> crit,
                             std::vector<NOMAD::EvalPoint> &evalPointList,
                             int maxEvalPoints) const
{
    verifyPointComplete(X);
    verifyPointSize(X);
    evalPointList.clear();

    bool stopWhenMaxFound = (maxEvalPoints > 0);
    bool errSizeDisplayed = false;  // Error about size to be displayed only once.
    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        if (X.size() != it->size())
        {
            if (!errSizeDisplayed)
            {
                std::string err = "CacheSet: find: Looking for a point of size ";
                err += NOMAD::itos(X.size());
                err += " but the cache points are of size ";
                err += NOMAD::itos(it->size());
                std::cout << "Warning: CacheSet: find: Looking for a point of size " << X.size() << " but found cache point of size " << it->size() << std::endl;
                errSizeDisplayed = true;
            }
            continue; // Points are in different dimensions -skip.
        }

        if (crit(X, *it))
        {
            const NOMAD::EvalPoint& evalPoint(*it);
            evalPointList.push_back(evalPoint);
            if (stopWhenMaxFound && evalPointList.size() >= (size_t)maxEvalPoints)
            {
                break;
            }
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    return evalPointList.size();
}


size_t NOMAD::CacheSet::find(std::function<bool(const NOMAD::EvalPoint&)> crit,
                             std::vector<NOMAD::EvalPoint> &evalPointList) const
{
    evalPointList.clear();
    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        const NOMAD::EvalPoint& evalPoint(*it);
        if (crit(evalPoint))
        {
            evalPointList.push_back(evalPoint);
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP

    return evalPointList.size();
}


void NOMAD::CacheSet::browse(std::function<void(const NOMAD::EvalPoint&)> crit) const
{

    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        const NOMAD::EvalPoint& evalPoint(*it);
        crit(evalPoint);
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
}

size_t NOMAD::CacheSet::find(std::function<bool(const NOMAD::EvalPoint&)> crit1,
                             std::function<bool(const NOMAD::EvalPoint&)> crit2,
                             std::vector<NOMAD::EvalPoint> &evalPointList) const
{
    evalPointList.clear();

    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        if ( crit1(*it) && crit2(*it) )
        {
            const NOMAD::EvalPoint& evalPoint(*it);
            evalPointList.push_back(evalPoint);
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    return evalPointList.size();
}


size_t NOMAD::CacheSet::findBestFeas(std::vector<EvalPoint> &evalPointList,
                                     const Point& fixedVariable,
                                     const FHComputeType & completeComputeType) const
{
    evalPointList.clear();
    auto compactComputeType = completeComputeType.Short();
    auto computeType = compactComputeType.computeType;
    auto evalType = completeComputeType.evalType;

    size_t nobj = 0;
    for (const auto & bbo: getBbOutputType())
    {
        if (bbo.isObjective())
        {
            nobj += 1;
        }
    }

    if (((computeType == ComputeType::STANDARD) && (nobj == 1)) ||
        ((computeType == ComputeType::DMULTI_COMBINE_F) && (nobj > 1 )) ||
        computeType == ComputeType::PHASE_ONE ||
        computeType == ComputeType::UNDEFINED ||
        computeType == ComputeType::USER      )
    {
        findBest(NOMAD::Eval::compEvalFindBest, evalPointList, true, 0,
                 fixedVariable, completeComputeType);
        return evalPointList.size();
    }
    
    
    std::list<NOMAD::EvalPoint> tmpEvalPointList;
    NOMAD::EvalPointSet::const_iterator itCache;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (itCache = _cache.begin(); itCache != _cache.end(); ++itCache)
    {
        const NOMAD::EvalPoint& evalPoint(*itCache);
        const NOMAD::Eval* eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (!eval->isFeasible(compactComputeType))
        {
            continue;
        }
        // Must be in the subspace defined byFixedVariable
        if (!evalPoint.hasFixed(fixedVariable))
        {
            continue;
        }
        // For robustness, be sure the cache picks up points which
        // have the same number of objectives
        size_t nobjEval = 0;
        for (const auto & bbo: eval->getBBOutputTypeList())
        {
            if (bbo.isObjective())
            {
                nobjEval += 1;
            }
        }
        if (nobjEval != nobj)
        {
            continue;
        }

        // Found first point
        if (tmpEvalPointList.empty())
        {
            tmpEvalPointList.push_back(evalPoint);
        }
        else
        {
            // Two cases:
            // 1- biobjective: points are ordered by lexicographic order.
            // Finding and removing dominated points is extremely efficient.
            //
            // See Algorithm 2 of
            //
            // A. Jaszkiewicz and T. Lust,
            // "ND-Tree-Based Update: A Fast Algorithm for the Dynamic Nondominance Problem,"
            // IEEE Transactions on Evolutionary Computation,
            // vol. 22, no. 5, pp. 778-791, Oct. 2018,
            // doi: 10.1109/TEVC.2018.2799684.
            //
            // One could also simply order the points by lexicographic order with one pass to get
            // all non dominated ones.
            //
            if (nobj == 2)
            {
                bool insert = false;
                auto isBelowf1Eval = [&evalType, &compactComputeType, eval](const EvalPoint& ev)
                {
                    return ev.getEval(evalType)->getFs(compactComputeType)[0] <= eval->getFs(compactComputeType)[0];
                };
                // Find the last element of the list which satisfies the condition
                auto itPfreverse = std::find_if(tmpEvalPointList.rbegin(), tmpEvalPointList.rend(), isBelowf1Eval);
                std::list<EvalPoint>::iterator itPfforward;

                if (itPfreverse == tmpEvalPointList.rend())
                {
                    // In this case, evalPoint has the smallest f1 value of the list
                    // and can be inserted at the beginning.
                    insert = true;
                }
                else
                {
                    // Check that evalPoint is non dominated
                    if (eval->getFs(compactComputeType)[1] < itPfreverse->getFs(completeComputeType)[1])
                    {
                        insert = true;
                        // Two subcases
                        // 1- evalPoint dominates itPfreverse element: will be inserted before
                        // all (potential) equal elements with itPfreverse values.
                        if (eval->getFs(compactComputeType)[0] == itPfreverse->getFs(completeComputeType)[0])
                        {
                            NOMAD::EvalPoint tmpEvalPoint(*itPfreverse);
                            // DO NOT UNDERSTAND: why when I do not create an EvalPoint, do I have a user rejected status ?
                            const NOMAD::Eval* evalTmp = tmpEvalPoint.getEval(evalType);

                            // Skip all equal elements.
                            itPfreverse++;
                            while (itPfreverse != tmpEvalPointList.rend())
                            {
                                NOMAD::EvalPoint tmp2EvalPoint(*itPfreverse);
                                const NOMAD::Eval* evalTmp2 = tmp2EvalPoint.getEval(evalType);
                                if ((evalTmp->getFs(compactComputeType)[0] != evalTmp2->getFs(compactComputeType)[0]) ||
                                    (evalTmp->getFs(compactComputeType)[1] != evalTmp2->getFs(compactComputeType)[1]))
                                {
                                    break;
                                }
                                itPfreverse++;
                            }
                        }
                        // 2- evalPoint is non dominated: will be inserted after itPfreverse element.
                    }
                    // or equal
                    else if ((eval->getFs(compactComputeType)[0] == itPfreverse->getFs(completeComputeType)[0]) &&
                             (eval->getFs(compactComputeType)[1] == itPfreverse->getFs(completeComputeType)[1]))
                    {
                        // evalPoint will be inserted after itPfreverse element
                        insert = true;
                    }
                }
                if (insert)
                {
                    // Add new evalPoint
                    tmpEvalPointList.insert(itPfreverse.base(), evalPoint);

                    // Remove points after evalPoint
                    itPfforward = itPfreverse.base();
                    while (itPfforward != tmpEvalPointList.end())
                    {
                        // evalj element is dominated.
                        const NOMAD::Eval* evalj = itPfforward->getEval(evalType);
                        if (eval->getFs(compactComputeType)[1] <= evalj->getFs(compactComputeType)[1])
                        {
                            tmpEvalPointList.erase(itPfforward++);
                            continue;
                        }
                        itPfforward++;
                    }
                }
            }
            // 2- More than two objectives. In this case, no order structure is exploitable.
            else
            {
                bool insert = true;
                auto itPf = tmpEvalPointList.begin();
                while (itPf != tmpEvalPointList.end())
                {
                    auto compFlag = evalPoint.compMO(*itPf, completeComputeType);
                    if (compFlag == CompareType::DOMINATED)
                    {
                        insert = false;
                        break;
                    }
                    if (compFlag == CompareType::DOMINATING)
                    {
                        tmpEvalPointList.erase(itPf++);
                        continue;
                    }
                    itPf++;
                }
                if (insert)
                {
                    tmpEvalPointList.push_front(evalPoint);
                }
            }
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    std::copy(tmpEvalPointList.begin(), tmpEvalPointList.end(), std::back_inserter(evalPointList));
    return evalPointList.size();
}

// Find best infeasible points with h<hmax:  least infeasible points with smallest f -> index 0 (and above), best f points with smallest h-> last index (and below). 
// All best f points have the same bboutputs. Idem for the least infeasible points.
size_t NOMAD::CacheSet::findBestInf(std::vector<EvalPoint> &evalPointList,
                                    const Double& hMax,
                                    const Point& fixedVariable,
                                    const FHComputeType& completeComputeType) const
{
    evalPointList.clear();
    
    auto evalType = completeComputeType.evalType;
    auto compactComputeType = completeComputeType.Short();

    size_t nobj = 0;
    for (const auto & bbo: getBbOutputType())
    {
        if (bbo.isObjective())
        {
            nobj += 1;
        }
    }

    NOMAD::EvalPointSet::const_iterator itCache;
    
    // Refs values (f and h) for both bestF and leastInf
    NOMAD::ArrayOfDouble bestFRefFs(nobj,NOMAD::INF);
    NOMAD::Double bestFRefH(NOMAD::INF);
    NOMAD::Double leastInfRefH(NOMAD::INF);
    NOMAD::ArrayOfDouble leastInfRefFs(nobj,NOMAD::INF);
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (itCache = _cache.begin(); itCache != _cache.end(); ++itCache)
    {
        const NOMAD::EvalPoint& evalPoint(*itCache);
        const NOMAD::Eval* eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (eval->isFeasible(compactComputeType))
        {
            continue;
        }
        NOMAD::Double h = eval->getH(compactComputeType);
        if (!h.isDefined() || h > hMax || h == NOMAD::INF)
        {
            continue;
        }
        // Must be in the subspace defined byFixedVariable
        if (!evalPoint.hasFixed(fixedVariable))
        {
            continue;
        }
        // For robustness, be sure the cache picks up points which
        // have the same number of objectives
        size_t nobjEval = 0;
        for (const auto &bbo: eval->getBBOutputTypeList())
        {
            if (bbo.isObjective())
            {
                nobjEval += 1;
            }
        }
        if (nobjEval != nobj)
        {
            continue;
        }
        NOMAD::ArrayOfDouble fs = eval->getFs(compactComputeType);
        
        // Two types of best inf but no duplication of points. If leastInf and bestF are the same we put single point in the list (see below in the second step).
        
        // Better f (still infeasible though)
        // For multiobjective, compare all objectives in the arrayOfDouble (no dominance).
        if (fs.isComplete() && fs < bestFRefFs )
        {
            bestFRefFs = fs;
            bestFRefH = h;
        }
        
        // lower infeas (do not care about f)
        if (h < leastInfRefH)
        {
            leastInfRefH = h;
            leastInfRefFs = fs;
        }
    }
    
    // Create the list with bestF (last index and below if multiple point) and leastInf (index 0 and above if multiple points)
    for (itCache = _cache.begin(); itCache != _cache.end(); ++itCache)
    {
        // Must be eval ok
        const NOMAD::Eval* eval = itCache->getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        // Must be in the subspace defined byFixedVariable
        if (!itCache->hasFixed(fixedVariable))
        {
            continue;
        }

        NOMAD::ArrayOfDouble fs = eval->getFs(compactComputeType);
        NOMAD::Double h = eval->getH(compactComputeType);
        if (fs == bestFRefFs && h == bestFRefH)
        {
            evalPointList.push_back(*itCache);
            continue;
        }
        if (h == leastInfRefH && fs == leastInfRefFs)
        {
            evalPointList.insert(evalPointList.begin(),*itCache);
        }
    }

#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    return evalPointList.size();
}


size_t NOMAD::CacheSet::findFilterInf(std::vector<NOMAD::EvalPoint> &evalPointList,
                                      const Double& hMax,
                                      const Point& fixedVariable,
                                      const FHComputeType & completeComputeType) const
{
    auto evalType = completeComputeType.evalType;
    auto compactComputeType = completeComputeType.Short();
    auto computeType = compactComputeType.computeType;
    
    evalPointList.clear();

    size_t nobj = 0;
    for (const auto & bbo: getBbOutputType())
    {
        if (bbo.isObjective())
        {
            nobj += 1;
        }
    }

    // Notion of dominance is not defined for those compute types.
    // Let's find some best infeasible points with h<hmax: least infeasible point -> index 0 and best f point -> last index
    if (computeType == ComputeType::PHASE_ONE ||
        computeType == ComputeType::UNDEFINED ||
        computeType == ComputeType::USER)
    {
        // NB: not the most efficient version... vector->list!!!
        findBestInf(evalPointList, hMax, fixedVariable, completeComputeType);
        return evalPointList.size();
    }

    std::list<NOMAD::EvalPoint> tmpEvalPointList;
    NOMAD::EvalPointSet::const_iterator itCache;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (itCache = _cache.begin(); itCache != _cache.end(); ++itCache)
    {
        const NOMAD::EvalPoint& evalPoint(*itCache);
        const NOMAD::Eval* eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (eval->isFeasible(compactComputeType)){
            continue;
        }
        NOMAD::Double h = eval->getH(compactComputeType);
        if (!h.isDefined() || h > hMax || h == NOMAD::INF)
        {
            continue;
        }
        // Must be in the subspace defined byFixedVariable
        if (!evalPoint.hasFixed(fixedVariable))
        {
            continue;
        }
        // For robustness, be sure the cache picks up points which
        // have the same number of objectives
        size_t nobjEval = 0;
        for (const auto & bbo: eval->getBBOutputTypeList())
        {
            if (bbo.isObjective())
            {
                nobjEval += 1;
            }
        }
        if (nobjEval != nobj)
        {
            continue;
        }
        // The set of non dominated points is empty, so insert it.
        if (tmpEvalPointList.empty())
        {
            tmpEvalPointList.push_back(evalPoint);
        }
        else
        {
            // Insertion into a non-empty set.
            bool insert = true;
            auto itInfPf = tmpEvalPointList.begin();
            while (itInfPf != tmpEvalPointList.end())
            {
                auto compFlag = evalPoint.compMO(*itInfPf, completeComputeType, false);
                if (compFlag == NOMAD::CompareType::DOMINATED)
                {
                    insert = false;
                    break;
                }
                else if (compFlag == NOMAD::CompareType::DOMINATING)
                {
                    tmpEvalPointList.erase(itInfPf++);
                    continue;
                }
                itInfPf++;
            }
            if (insert)
            {
                tmpEvalPointList.insert(tmpEvalPointList.begin(),evalPoint);
            }
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    std::copy(tmpEvalPointList.begin(), tmpEvalPointList.end(), std::back_inserter(evalPointList));
    return evalPointList.size();
}





// Update EvalPoint in cache.
// Look for Point and update the Eval part.
// If the point is not found, it is a serious issue.
// Don't throw an exception, to avoid a crash, but write a warning.
// Returns true if update succeeded, false if there was an error.
bool NOMAD::CacheSet::update(const NOMAD::EvalPoint& evalPoint, NOMAD::EvalType  evalType, const NOMAD::MeshBasePtr mesh)
{
    bool updateOk = false;

    if (nullptr == evalPoint.getEval(evalType))
    {
        // Cannot update to a null Eval. Warn the user.
        std::string err = "Warning: CacheSet: Update: Cannot update to a NULL Eval for Point ";
        err += evalPoint.displayAll();
        std::cout << err << std::endl;
        return false;
    }

    NOMAD::EvalPointSet::const_iterator it;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    it = _cache.find(evalPoint);
    if (it == _cache.end())
    {
        std::string err = "Warning: CacheSet: Update: Did not find EvalPoint to update in cache: " + evalPoint.displayAll();
        std::cout << err << std::endl;
        NOMAD::OutputQueue::Add(err, NOMAD::OutputLevel::LEVEL_WARNING);
    }
    else
    {
        // Update EvalPoint in cache directly.
        // Since we are not changing the Point part, which is the only part
        // used for sorting, the cache should remain coherent.
        auto cacheEvalPoint = const_cast<NOMAD::EvalPoint*>(&*it);
        cacheEvalPoint->setEval(*evalPoint.getEval(evalType), evalType);
        if (NOMAD::EvalType::BB == evalType)
        {
            cacheEvalPoint->setNumberBBEval(evalPoint.getNumberBBEval());
        }
        if (nullptr != mesh)
        {
            cacheEvalPoint->setMesh(mesh);
        }
        
        // Update revealing status of the point (DiscoMads algorithm)
        cacheEvalPoint->setRevealingStatus(evalPoint.getRevealingStatus());
        
        // Update user fail eval check flag of the point (DiscoMads algorithm)
        cacheEvalPoint->setUserFailEvalCheck(evalPoint.getUserFailEvalCheck());

        updateOk = true;
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP

    return updateOk;
}


// Empty the cache and reset number of cache hits
bool NOMAD::CacheSet::clear()
{
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    _cache.clear();
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP

    // Note: We might not want to reset - in that case, remove this line.
    resetNbCacheHits();

    // Reset size of points in cache
    _n = 0;

    return true;
}


// Clear all quad and sgtelib model evaluations from the cache
void NOMAD::CacheSet::clearModelEval(const int mainThreadNum)
{
    processOnAllPoints(NOMAD::EvalPoint::clearModelEval, mainThreadNum);
}


// Purge the cache for space.
//
// Current implementation:
// For now, we only keep EvalPoints which have an f under the mean f.
// We should also keep points that are cache hits, or that were recent cache hits.
// It could be related to the iteration where the point was added,
// or to the number of purges it survived, etc.
// We might need a CachePoint class to keep that information.
// Etc.
//
// General idea: We want to keep points that might be hit again, to
// avoid having to recompute them.
// We also want to keep points that are good enough to be interesting to the
// user.
//
// Note June 2021: We are now ignoring points for which eval status is not EVAL_OK.
void NOMAD::CacheSet::purge()
{
    std::cout << "Warning: Calling Cache purge. Size is " << _cache.size() << " max is " << _maxSize << ". Some points will be removed from the cache." << std::endl;
    if ( _maxSize== NOMAD::INF_SIZE_T || _cache.size() < _maxSize)
    {
        // Do nothing
        return;
    }
    size_t nbRemovedLast = 1;

#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP

    while (_cache.size() >= _maxSize)
    {
        NOMAD::EvalPointSet tmpCache;
        NOMAD::Double meanF;
        size_t nbElemWithF = computeMeanF(meanF);
        //std::cout << "Debug: purge: meanF = " << meanF << " nb elem = " << nbElemWithF << std::endl;

        if (nbElemWithF > 0 && nbRemovedLast > 0)
        {
            // Remove all EvalPoints for which f is over or equal to the mean.
            // For this, use a temporary set/cache, because we
            // cannot iterate over a set and erase items at the same time.
            NOMAD::EvalPointSet::const_iterator it;
            for (it = _cache.begin(); it != _cache.end(); ++it)
            {
                if (NOMAD::EvalStatusType::EVAL_OK != it->getEvalStatus(NOMAD::EvalType::BB))
                {
                    continue;
                }
                if (!it->getF(defaultFHComputeType).isDefined())
                {
                    continue;
                }
                if (it->getF(defaultFHComputeType) < meanF)
                {
                    //std::cout << "Debug: purge: insert EvalPoint with f = " << it->getF(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD) << " to tmpCache" << std::endl;
                    tmpCache.insert(*it);
                }
                else
                {
                    //std::cout << "Debug: purge: Do not insert EvalPoint with f = " << it->getF(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD) << " to tmpCache" << std::endl;
                }
            }
        }
        else
        {
            // Remove arbitrary half the elements of cache.
            // Keep the first half.
            size_t i = 0;
            NOMAD::EvalPointSet::const_iterator it;
            for (it = _cache.begin(); i < _cache.size() / 2; ++it, i++)
            {
                tmpCache.insert(*it);
            }
        }

        // If tmpCache is empty, set nbRemovedLast to 0 and the next loop
        // will go to the "else" case.
        if (tmpCache.empty())
        {
            nbRemovedLast = 0;
        }
        else
        {
            //std::cout << "Debug: Cache has     " << _cache.size() << " elements." << std::endl;
            //std::cout << "Debug: Tmp cache has " << tmpCache.size() << " elements." << std::endl;
            nbRemovedLast = _cache.size() - tmpCache.size();
            _cache.clear();
            _cache = std::move(tmpCache);
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
}


// Naive way to compute the mean f for all points in the cache.
size_t NOMAD::CacheSet::computeMeanF(NOMAD::Double &mean) const
{
    size_t nbElem = 0;
    NOMAD::Double total = 0;
    mean.reset();
    NOMAD::EvalPointSet::const_iterator it;
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        if (NOMAD::EvalStatusType::EVAL_OK != it->getEvalStatus(NOMAD::EvalType::BB))
        {
            continue;
        }
        NOMAD::Double f = it->getF(defaultFHComputeType);
        if (f.isDefined())
        {
            total += f;
            nbElem++;
        }
    }
    if (nbElem > 0)
    {
        mean = total / (double)nbElem;
    }

    return nbElem;
}


// Call function func on all points generated by mainThreadNum
void NOMAD::CacheSet::processOnAllPoints(void (*func)(NOMAD::EvalPoint&), const int mainThreadNum)
{
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (const auto& it : _cache)
    {
        auto evalPoint = const_cast<NOMAD::EvalPoint*>(&it);
        if (   -1 == mainThreadNum 
            || evalPoint->getThreadAlgo() == mainThreadNum)
        {
            func(*evalPoint);
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
}


void NOMAD::CacheSet::deleteModelEvalOnly(const int mainThreadNum)
{
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    for (auto it = _cache.begin(); it != _cache.end();)
    {
        if (mainThreadNum != it->getThreadAlgo())
        {
            it++;
        }
        else
        {
            bool foundOtherEval = false;
            for (size_t i = 0; (i < (size_t)NOMAD::EvalType::LAST && !foundOtherEval); i++)
            {
                auto evalType = NOMAD::EvalType(i);
                if (NOMAD::EvalType::MODEL != evalType && nullptr != it->getEval(evalType))
                {
                    foundOtherEval = true;
                }
            }
            if (foundOtherEval)
            {
                it++;
            }
            else
            {
                // Only MODEL evaluation, or no evaluation, for this point.
                it = _cache.erase(it);
            }
        }
    }
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
}


// Write cache to file _filename
// This function will use operator<< defined below
bool NOMAD::CacheSet::write() const
{
    OUTPUT_INFO_START
    std::string s = "Write cache file " + _filename;
    NOMAD::OutputQueue::Add(s);
    OUTPUT_INFO_END
    return NOMAD::write(*this, _filename);
}


// Read _filename as written by write(), and add the points to the cache.
// This function will use operator>> defined below.
bool NOMAD::CacheSet::read()
{
    bool fileRead = false;
    if (NOMAD::checkReadFile(_filename))
    {
        OUTPUT_INFO_START
        std::string s = "Read cache file " + _filename;
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_NORMAL);
        OUTPUT_INFO_END
        fileRead = NOMAD::read(*this, _filename);
    }
    return fileRead;
}


// Display all points in cache
// Useful mostly for debugging purposes
std::string NOMAD::CacheSet::displayAll() const
{
    std::string retStr;
    for (const auto& evalPoint : _cache)
    {
        retStr += evalPoint.displayAll() + "\n";
    }

    return retStr;
}


// Display only EvalPoints that have a BB or SURROGATE eval that is good. Only eval status is checked.
// This method is used to write points to cache.
std::ostream& NOMAD::CacheSet::displayPointsWithEval(std::ostream& os) const
{
    for (const auto& evalPoint : _cache)
    {
        if ( (nullptr != evalPoint.getEval(NOMAD::EvalType::BB) && evalPoint.getEval(NOMAD::EvalType::BB)->goodForCacheFile() ) ||
            (nullptr != evalPoint.getEval(NOMAD::EvalType::SURROGATE) && evalPoint.getEval(NOMAD::EvalType::SURROGATE)->goodForCacheFile() ) )
        {
            os << evalPoint.displayForCache(_bbEvalFormat) << std::endl;
        }
    }

    return os;
}

void NOMAD::CacheSet::moveEvalPointToCacheForRerun()
{
    _cacheForRerun = _cache;
    _cache.clear();
}

// Display only EvalPoints that have an eval.
std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::CacheSet& cache)
{
    os << "CACHE_HITS " << cache.getNbCacheHits() << std::endl;
    os << "BB_OUTPUT_TYPE " << cache.getBbOutputType() << std::endl;
    cache.displayPointsWithEval(os);

    return os;
}


// Get these EvalPoints from stream
std::istream& NOMAD::operator>>(std::istream& is, NOMAD::CacheSet& cache)
{
    std::string s;
    NOMAD::BBOutputTypeList bbOutputTypes;

    is >> s;
    if ("CACHE_HITS" == s)
    {
        size_t cacheHits;
        is >> cacheHits;
        cache.setNbCacheHits(cacheHits);
    }
    else
    {
        // Put back s to istream.
        for (unsigned i = 0; i < s.size(); i++)
        {
            is.unget();
        }
    }

    is >> s;
    if ("BB_OUTPUT_TYPE" == s)
    {
        while (is >> s && is.good() && !is.eof())
        {
            if (NOMAD::ArrayOfDouble::pStart == s)
            {
                is.unget();
                break;
            }
            else
            {
                bbOutputTypes.emplace_back(s);
            }
        }

        cache.setBBOutputType(bbOutputTypes);
    }


    NOMAD::EvalPoint evalPoint;
    while (is >> evalPoint && is.good() && !is.eof())
    {
        evalPoint.setBBOutputType(bbOutputTypes);
        evalPoint.updateTag();
        evalPoint.setEvalIsFromCacheFile(true);
        cache.insert(evalPoint);
    }

    return is;
}


