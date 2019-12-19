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
/**
 \file   CacheSet.cpp
 \brief  Simple implementation of Cache derived from CacheBase (implementation)
 \author Viviane Rochon Montplaisir
 \date   January 2018
 \see    CacheSet.hpp
 */
#include "../Cache/CacheSet.hpp"
#include "../Math/Point.hpp"
#include "../Output/OutputQueue.hpp"

#include <fstream>
#include <iostream>

// Init static members
NOMAD::BBOutputTypeList NOMAD::CacheSet::_bbOutputType = NOMAD::BBOutputTypeList();
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
#ifdef _OPENMP
    omp_init_lock(&_cacheLock);
#endif // _OPENMP
}


// Terminate CacheSet class.
// To be called by the Destructor.
void NOMAD::CacheSet::destroy()
{
    // Clear the _cache directly.
    // No need to set lock, assuming there is only one cache and
    // that now it is the end of the run and we are calling its destructor.
    _cache.clear();

#ifdef _OPENMP
    omp_destroy_lock(&_cacheLock);
#endif // _OPENMP
}

void NOMAD::CacheSet::setInstance(const std::shared_ptr<NOMAD::CacheParameters>& cacheParams,
                                  const BBOutputTypeList& bbOutputType)
{
#ifdef _OPENMP
    // Lock cache set before creating instance
    omp_set_lock(&_cacheLock);
#endif

    if ( _single == nullptr )
    {
        _single = std::unique_ptr<NOMAD::CacheSet>(new CacheSet(cacheParams)) ;
    }
    else
    {
        std::string err = "Cannot get instance. NOMAD::CacheSet::setInstance must be called only ONCE before calling NOMAD::CacheBase::getInstance()" ;
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

#ifdef _OPENMP
    omp_unset_lock(& _cacheLock);
#endif // _OPENMP

    _bbOutputType = bbOutputType;

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
    if (_cache.size() > 0 && _n != point.size())
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
bool NOMAD::CacheSet::insert(const NOMAD::EvalPoint &evalPoint)
{
    verifyPointComplete(evalPoint);
    verifyPointSize(evalPoint);
    // First insert sets n (even if insert fails)
    if (0 == _cache.size())
    {
        _n = evalPoint.size();
    }

    std::pair<EvalPointSet::iterator,bool> ret;
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    ret = _cache.insert(evalPoint);
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    bool inserted = ret.second;
    if (!inserted)
    {
        _nbCacheHits++;
    }

    return inserted;
}


// Get EvalPoint evalPoint at Point x from the cache
// Returns the number of EvalPoints found
size_t NOMAD::CacheSet::find(const NOMAD::Point& x, NOMAD::EvalPoint &evalPoint) const
{
    size_t nbFound = 0;
    EvalPointSet::const_iterator it;
    it = _cache.find(NOMAD::EvalPoint(x));
    if (it != _cache.end())
    {
        evalPoint = *it;
        nbFound = 1;
    }
    return nbFound;
}


// Insert evalPoint in cache.
// Return a boolean indicating if we should eval this point.
// If insertion worked, the point was not in the cache before. Return true.
// If insertion did not work, the point was in the cache before.
// Depending on its EvalStatus, return true if it should be evaluated again,
// false otherwise.
bool NOMAD::CacheSet::smartInsert(const NOMAD::EvalPoint &evalPoint,
                                  const short maxNumberEval,
                                  const EvalType& evalType)
{
    verifyPointComplete(evalPoint);
    verifyPointSize(evalPoint);

    // First insert sets n (even if insert fails)
    if (0 == _cache.size())
    {
        _n = evalPoint.size();
    }

    bool inserted = false;
    bool doEval;
    std::pair<EvalPointSet::iterator,bool> ret;   // Return of the insert()
#ifdef _OPENMP
    omp_set_lock(&_cacheLock);
#endif // _OPENMP
    ret = _cache.insert(evalPoint);
#ifdef _OPENMP
    omp_unset_lock(&_cacheLock);
#endif // _OPENMP
    inserted = ret.second;
    bool canEval = (*ret.first).toEval(maxNumberEval, evalType);
    if (inserted && canEval)
    {
        doEval = true;
    }
    else if (nullptr == (ret.first)->getEval(evalType))
    {
        // Point already inserted, but not evaluated.
        // NOTE: We do not know if this point is in the evaluation queue yet, or not.
        // We do not know here if the evaluation queue is cleared between runs.
        // TODO: Clear this out. For now, set doEval to true, but this could
        // mean that the point is evaluated twice.

        // Only warn outside of sgte context.
        if (NOMAD::EvalType::BB == evalType)
        {
            std::string s = "Point already inserted in cache, but not evaluated: ";
            s += ret.first->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
        }

        doEval = true;
    }
    else
    {
        // Cache hit.
        doEval = canEval;

        // Only count as cache hit when using Blackbox Eval.
        if (NOMAD::EvalType::BB == evalType)
        {
            _nbCacheHits++;
            std::string s = "Cache hit: ";
            s += ret.first->display();
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
        }
        if (doEval)
        {
            std::cerr << "Warning: CacheSet: smartInsert: New evaluation of point found in cache " << (*ret.first).display() << std::endl;
        }
    }

    return doEval;
}


size_t NOMAD::CacheSet::find(const NOMAD::Point x, std::vector<NOMAD::EvalPoint> &evalPointList) const
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
                             bool (*comp)(const NOMAD::Eval&, const NOMAD::Eval&),
                             std::vector<NOMAD::EvalPoint> &evalPointList,
                             const EvalType& evalType) const
{
    evalPointList.clear();
    EvalPointSet::const_iterator it;
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        const NOMAD::Eval* eval = it->getEval(evalType);
        if (nullptr == eval)
        {
            continue;
        }
        if (comp(*eval, refeval))
        {
            NOMAD::EvalPoint evalPoint(*it);
            evalPointList.push_back(evalPoint);
        }
    }

    return evalPointList.size();
}


// Get best eval points, using comp()
size_t NOMAD::CacheSet::findBest(bool (*comp)(const NOMAD::Eval&, const NOMAD::Eval&),
                     std::vector<NOMAD::EvalPoint> &evalPointList,
                     const bool findFeas,
                     const NOMAD::Double& hMax,
                     const NOMAD::Point& fixedVariable,
                     const EvalType& evalType) const
{
    evalPointList.clear();
    EvalPointSet::const_iterator it;
    auto refeval = std::shared_ptr<NOMAD::Eval>(nullptr);

    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        NOMAD::EvalPoint evalPoint(*it);
        const NOMAD::Eval* eval = evalPoint.getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (findFeas != eval->isFeasible())
        {
            continue;
        }
        if (eval->getH() > hMax)
        {
            continue;
        }
        // Must be in the sub-space defined by fixedVariable
        if (!evalPoint.hasFixed(fixedVariable))
        {
            continue;
        }

        if (nullptr == refeval)
        {
            // Found first point
            refeval = std::make_shared<NOMAD::Eval>(*eval);
            evalPointList.push_back(evalPoint);
        }
        else if (*eval == *refeval)
        {
            // Found another point with eval == refeval
            evalPointList.push_back(evalPoint);
        }
        else if (comp(*eval, *refeval))
        {
            // Found a better point
            *refeval = *eval;
            // Reset list with new best
            evalPointList.clear();
            evalPointList.push_back(evalPoint);
        }
    }

    return evalPointList.size();
}


size_t NOMAD::CacheSet::findBestFeas(std::vector<NOMAD::EvalPoint> &evalPointList,
                                     const NOMAD::Point& fixedVariable,
                                     const EvalType& evalType) const
{
    findBest(NOMAD::Eval::compEvalFindBest, evalPointList, true, 0, fixedVariable, evalType);
    return evalPointList.size();
}


bool NOMAD::CacheSet::hasFeas(const EvalType& evalType) const
{
    bool ret = false;

    for (auto it = _cache.begin(); it != _cache.end(); ++it)
    {
        const NOMAD::Eval* eval = (*it).getEval(evalType);
        if (nullptr == eval || NOMAD::EvalStatusType::EVAL_OK != eval->getEvalStatus())
        {
            continue;
        }
        if (eval->isFeasible())
        {
            ret = true;
            break;
        }
    }

    return ret;
}


size_t NOMAD::CacheSet::findBestInf(std::vector<NOMAD::EvalPoint> &evalPointList,
                                    const NOMAD::Double& hMax,
                                    const NOMAD::Point& fixedVariable,
                                    const EvalType& evalType) const
{
    findBest(NOMAD::Eval::compEvalFindBest, evalPointList, false, hMax, fixedVariable, evalType);

    return evalPointList.size();
}



size_t NOMAD::CacheSet::find(NOMAD::Point X,
                 NOMAD::Double distance,
                 std::vector<NOMAD::EvalPoint> &evalPointList,
                 int maxEvalPoints) const
{
    verifyPointComplete(X);
    verifyPointSize(X);
    evalPointList.clear();

    bool stopWhenMaxFound = (maxEvalPoints > 0);
    bool errSizeDisplayed = false;  // Error about size to be displayed only once.
    EvalPointSet::const_iterator it;
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        if (X.size() != it->size())
        {
            if (!errSizeDisplayed)
            {
                std::string err = "CacheSet: find: Looking for a point of size ";
                err += X.size();
                err += " but the cache points are of size ";
                err += it->size();
                std::cerr << "Warning: CacheSet: find: Looking for a point of size " << X.size() << " but found cache point of size " << it->size() << std::endl;
                errSizeDisplayed = true;
            }
            continue; // Points are in different dimensions -skip.
        }

        if (NOMAD::Point::dist(X, *it) <= distance)
        {
            NOMAD::EvalPoint evalPoint(*it);
            evalPointList.push_back(evalPoint);
            if (stopWhenMaxFound && evalPointList.size() >= (uint)maxEvalPoints)
            {
                break;
            }
        }
    }
    return evalPointList.size();
}


size_t NOMAD::CacheSet::find(bool (*crit)(const NOMAD::EvalPoint&),
                     std::vector<NOMAD::EvalPoint> &evalPointList) const
{
    evalPointList.clear();
    EvalPointSet::const_iterator it;
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        NOMAD::EvalPoint evalPoint(*it);
        if (crit(evalPoint))
        {
            evalPointList.push_back(evalPoint);
        }
    }

    return evalPointList.size();
}


// Update EvalPoint in cache.
// Look for Point and update the Eval part.
// If the point is not found, it is a serious issue.
// Don't throw an exception, to avoid a crash, but write a warning.
// Returns true if update succeded, false if there was an error.
bool NOMAD::CacheSet::update(const NOMAD::EvalPoint& evalPoint, const EvalType& evalType)
{
    bool updateOk = false;

    if (nullptr == evalPoint.getEval(evalType))
    {
        // Cannot update to a null Eval. Warn the user.
        std::string err = "Warning: CacheSet: Update: Cannot update to a NULL Eval for Point ";
        err += evalPoint.displayAll();
        std::cerr << err << std::endl;
        return false;
    }

    EvalPointSet::const_iterator it;
    it = _cache.find(evalPoint);
    if (it == _cache.end())
    {
        std::string err = "Warning: CacheSet: Update: Did not find EvalPoint to update in cache: " + evalPoint.displayAll();
        NOMAD::OutputQueue::Add(err, NOMAD::OutputLevel::LEVEL_WARNING);
    }
    else
    {
        // Update EvalPoint in cache directly.
        // Since we are not changing the Point part, which is the only part
        // used for sorting, the cache should remain coherent.
        auto cacheEvalPoint = const_cast<NOMAD::EvalPoint*>(&*it);
        cacheEvalPoint->setEval(*evalPoint.getEval(evalType), evalType);
        cacheEvalPoint->setNumberEval(evalPoint.getNumberEval());
        updateOk = true;
    }

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

    return true;
}


// Clear all sgte evaluations from the cache
void NOMAD::CacheSet::clearSgte()
{
    processOnAllPoints(NOMAD::EvalPoint::clearEvalSgte);
}


// Purge the cache for space.
//
// Current implementation:
// For now, we only keep EvalPoints which have an f under the mean f.
// We should also keep points that are cache hits, or that were recent cache hits.
// TODO: There should be a concept of time / oldness of points.
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
void NOMAD::CacheSet::purge()
{
    std::cerr << "Warning: Calling Cache purge. Size is " << _cache.size() << " max is " << _maxSize << ". Some points will be removed from the cache." << std::endl;
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
        EvalPointSet tmpCache;
        NOMAD::Double meanF;
        size_t nbElemWithF = computeMeanF(meanF);
        //std::cout << "Debug: purge: meanF = " << meanF << " nb elem = " << nbElemWithF << std::endl;

        if (nbElemWithF > 0 && nbRemovedLast > 0)
        {
            // Remove all EvalPoints for which f is over or equal to the mean.
            // For this, use a temporary set/cache, because we
            // cannot iterate over a set and erase items at the same time.
            EvalPointSet::const_iterator it;
            for (it = _cache.begin(); it != _cache.end(); ++it)
            {
                if (!it->getF(NOMAD::EvalType::BB).isDefined())
                {
                    continue;
                }
                if (it->getF(NOMAD::EvalType::BB) < meanF)
                {
                    //std::cout << "Debug: purge: insert EvalPoint with f = " << it->getF(NOMAD::EvalType::BB) << " to tmpCache" << std::endl;
                    tmpCache.insert(*it);
                }
                else
                {
                    //std::cout << "Debug: purge: Do not insert EvalPoint with f = " << it->getF(NOMAD::EvalType::BB) << " to tmpCache" << std::endl;
                }
            }
        }
        else
        {
            // Remove arbitrarly half the elements of cache.
            // Keep the first half.
            // TODO Review how is the set ordered, could we use that
            // order to better fill tmpCache (in both f and non-f cases) ?
            size_t i = 0;
            EvalPointSet::const_iterator it;
            for (it = _cache.begin(); i < _cache.size() / 2; ++it, i++)
            {
                tmpCache.insert(*it);
            }
        }

        // If tmpCache is empty, set nbRemovedLast to 0 and the next loop
        // will go to the "else" case.
        if (0 == tmpCache.size())
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
    EvalPointSet::const_iterator it;
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        NOMAD::Double f = it->getF(NOMAD::EvalType::BB);
        if (f.isDefined())
        {
            total += f;
            nbElem++;
        }
    }
    if (nbElem > 0)
    {
        mean = total / nbElem;
    }

    return nbElem;
}


// Call function func on all points in cache.
void NOMAD::CacheSet::processOnAllPoints(void (*func)(NOMAD::EvalPoint&))
{
    for (auto it = _cache.begin(); it != _cache.end(); ++it)
    {
        auto evalPoint = const_cast<NOMAD::EvalPoint*>(&*it);
        func(*evalPoint);
    }
}


// Write cache to file _filename
// This function will use operator<< defined below
bool NOMAD::CacheSet::write() const
{
    std::string s = "Write cache file " + _filename;
    NOMAD::OutputQueue::Add(s);
    return NOMAD::write(*this, _filename);
}


// Read _filename as written by write(), and add the points to the cache.
// This function will use operator>> defined below.
bool NOMAD::CacheSet::read()
{
    bool fileRead = false;
    if (NOMAD::checkReadFile(_filename))
    {
        std::string s = "Read cache file " + _filename;
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_NORMAL);
        fileRead = NOMAD::read(*this, _filename);
    }
    return fileRead;
}


// Display all points in cache
// Useful mostly for debugging purposes
std::ostream& NOMAD::CacheSet::displayAll(std::ostream& os) const
{
    EvalPointSet::const_iterator it;
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        os << *it;
    }

    return os;
}


// Display only EvalPoints that have a BB eval that is good. Only eval status
// is checked.
// This method is used to write points to cache.
std::ostream& NOMAD::CacheSet::displayPointsWithEval(std::ostream& os) const
{
    EvalPointSet::const_iterator it;
    for (it = _cache.begin(); it != _cache.end(); ++it)
    {
        NOMAD::EvalPoint evalPoint = *it;
        if (nullptr != evalPoint.getEval(NOMAD::EvalType::BB) && evalPoint.getEval(NOMAD::EvalType::BB)->goodForCacheFile())
        {
            os << evalPoint << std::endl;
        }
    }

    return os;
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
        NOMAD::BBOutputTypeList bbOutputTypes;

        while (is >> s && is.good() && !is.eof())
        {
            if (NOMAD::ArrayOfDouble::pStart == s)
            {
                is.unget();
                break;
            }
            else
            {
                bbOutputTypes.push_back(NOMAD::stringToBBOutputType(s));
            }
        }

        cache.setBBOutputType(bbOutputTypes);
    }


    NOMAD::EvalPoint evalPoint;
    while (is >> evalPoint && is.good() && !is.eof())
    {
        cache.insert(evalPoint);
    }

    // Need to recompute F and H on all cache points
    NOMAD::CacheBase::getInstance()->processOnAllPoints(NOMAD::CacheSet::recomputeFH);

    return is;
}


// Recompute f and h, using BB eval.
void NOMAD::CacheSet::recomputeFH(NOMAD::EvalPoint& evalPoint)
{
    auto eval = evalPoint.getEval(NOMAD::EvalType::BB);
    if (eval)
    {
        eval->setBBOutputAndRecompute(eval->getBBOutput(), _bbOutputType);
    }
}
