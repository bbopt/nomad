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
 * \file   CacheSet.hpp
 * \brief  Simple implementation of Cache derived from CacheBase, using set or unordered_set.
 * \author Viviane Rochon Montplaisir
 * \date   January 2018
 * \see    CacheSet.cpp
 */

#ifndef __NOMAD_4_5_CACHESET__
#define __NOMAD_4_5_CACHESET__

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP
#include "../Cache/CacheBase.hpp"
#include "../Eval/EvalPoint.hpp"

#include "../nomad_platform.hpp"
#include "../nomad_nsbegin.hpp"


/// Class implementing the abstract class \b CacheBase
/**
* Uses a set or unordered set of EvalPoint for the cache.
*/
class DLL_EVAL_API CacheSet : public CacheBase {

private:

    /// Lock for multithreading
#ifdef _OPENMP
    static omp_lock_t  _cacheLock;
#endif // _OPENMP

    static BBOutputTypeList    _bbOutputType;  ///< Corresponds to parameter BB_OUTPUT_TYPE used for this cache
    static ArrayOfDouble       _bbEvalFormat;  ///< Used to write cache correctly

    EvalPointSet _cache;  ///< The set of points that constitutes the cache.
    EvalPointSet _cacheForRerun;  ///< The set of points that constitutes the cache used for rerun only (empty if not in rerun mode). Filled with points from a cache file. Used for evaluation, not for "cache hit". 

    

    /// Constructor
    /**
     \param cacheParams     The parameters for cache -- \b IN.
     */
    explicit CacheSet(const std::shared_ptr<CacheParameters>& cacheParams)
      : CacheBase(cacheParams),
        _cache()
    {
        init();
    }


public:
    /*---------------*/
    /* Class Methods */
    /*---------------*/

    /// Destructor
    virtual ~CacheSet()
    {
        destroy();
    }

    /// Set the singleton (done once)
    /**
     \param cacheParams     The cache parameters -- \b IN.
     \param bbOutputType    List of the blackbox output type -- \b IN.
     \param bbEvalFormat    Format to write cache correctly -- \b IN.
     */
    static void setInstance(const std::shared_ptr<CacheParameters>& cacheParams,
                            const BBOutputTypeList& bbOutputType,
                            const ArrayOfDouble& bbEvalFormat = ArrayOfDouble());

    /// Get the list of blackbox output types
    static BBOutputTypeList getBbOutputType() { return _bbOutputType; }

    /// Set the list of blackbox output type
    /**
     \param bbOutputType    The list to use in cache -- \b IN.
     */
    static void setBBOutputType(const BBOutputTypeList& bbOutputType) { _bbOutputType = bbOutputType; }

    /// Add a new EvalPoint to the cache
    /**
     \param evalPoint   The eval point to insert in the cache     -- \b IN.
     \return            \c true if the insertion succeeded and \c false if the EvalPoint was already in the cache.
     */
    bool insert(const EvalPoint &evalPoint) override;

    /// Get eval point at point x from the cache (there can be only one in CacheSet).
    /**
     \param x           The point to find                   -- \b IN.
     \param evalPoint   The copy of the data in cache       -- \b OUT.
     \param evalType    If not UNDEFINED, wait for the point to be evaluated for this EvalType. -- \b IN.
     \param waitIfNotYetAvailable    Flag to control if we wait for the point to have an evaluation for this evaltype. -- \b IN.
     \return            An integer: 1 if found, 0 otherwise
     */
    size_t find(const Point & x, EvalPoint &evalPoint,
                const EvalType evalType = EvalType::UNDEFINED,
                bool waitIfNotYetAvailable = true ) const override;
    
    
    /// Get eval point at point x from the cache for rerun (there can be only one in CacheSet).
    /**
     \param x           The point to find                   -- \b IN.
     \param evalPoint   The returned eval point that matches x  -- \b IN/OUT.
     \return true if the evalPoint  found in cache for rerun, false otherwise.
     */
    bool findInCacheForRerun(const Point & x,
                             NOMAD::EvalPoint &evalPoint ) const override;

    /// Insert evalPoint in cache.
    /**
     * evalPoint's tag (mutable) is updated.
     * Return a boolean indicating if we should eval this point. \n
     * If insertion worked, the point was not in the cache before. Return true. \n
     * If insertion did not work, the point was in the cache before. \n
     * Depending on its EvalStatus, return true if it should be evaluated again, false otherwise.

     \param evalPoint       The point to insert                         -- \b IN.
     \param maxNumberEval   The max number of evaluations of the point  -- \b IN.
     \param evalType        Which eval of the EvalPoint to look at -- \b IN.
     \return                A boolean indicating if we should eval this point.
     */
    bool smartInsert(const EvalPoint &evalPoint,
                     const short maxNumberEval,
                     EvalType evalType ) override;

    /// Get eval point at point x from the cache and return it in a list.
    /**
     \param x                The point to find                              -- \b IN.
     \param evalPointList    The eval point corresponding to x in a list    -- \b OUT.
     \return                 1 if found, 0 otherwise
     */
    size_t find(const Point& x,
                std::vector<EvalPoint> &evalPointList) const override;

    /**
     * \brief Get all eval points for which comp(refeval) returns true.
     *
     * All eval points for which eval is inferior to refeval.

     \param refeval                The point to find                                              -- \b IN.
     \param comp                       The comparison function                                        -- \b IN.
     \param evalPointList   The eval points that verify comp()==true returned in a list    -- \b OUT.
     \param computeType   Which type of computation (eval type, compute type and h norm type)  -- \b IN.
     \return                 The number of eval points found.
     */
    size_t find(const Eval &refeval,
                std::function<bool(const Eval&, const Eval&, const FHComputeTypeS&)> comp,
                std::vector<EvalPoint> &evalPointList,
                const FHComputeType& computeType) const override;

    

    /**
     * \brief Get best eval points, using comp().
     *
     * Only the points with eval status EVAL_OK are considered.
     \param comp                        The comparison function                                    -- \b IN.
     \param evalPointList    The best eval points that verify comp()==true in a list    -- \b OUT.
     \param findFeas               Flag to find feasible points                               -- \b IN.
     \param hMax                        Maximum acceptable value for h, when findFeas is false     -- \b IN.
     \param fixedVariable    Searching for a subproblem defined by this point -- \b IN.
     \param computeType   Which type of computation (eval type, compute type and h norm type)  -- \b IN.
     \return                 The number of eval points found.
     */
    virtual size_t findBest(std::function<bool(const Eval&, const Eval&,const FHComputeTypeS&)> comp,
                            std::vector<EvalPoint> &evalPointList,
                            const bool findFeas,
                            const Double& hMax,
                            const Point& fixedVariable,
                            const FHComputeType& computeType) const override;


    /// Test if cache contains a feasible point.
    /**
     \param computeType   Which type of computation (eval type, compute type and h norm type)  -- \b IN.
     \return \c true if the cache contains at least one feasible point, \c false otherwise.
     */
    bool hasFeas(const FHComputeType & computeType) const override;
    
    /// Test if cache contains an infeasible point.
    /**
     \return \c true if the cache contains at least one infeasible point, \c false otherwise.
     */
    bool hasInfeas(const FHComputeType & computeType) const override;


    /// Get all eval points verifying a criterion function with respect to point X. The criterion function can be a measure of distance to X.
    /**

     \param X                               The point of reference                              -- \b IN.
     \param crit                        The criterion function                              -- \b IN.
     \param evalPointList    The eval points within the prescribed distance of X -- \b OUT.
     \param maxEvalPoints    The maximum number of points to select              -- \b IN.
     \return                 The number of eval points found.
     */
    virtual size_t find(const Point & X,
                        std::function<bool(const Point&, const EvalPoint &)> crit,
                        std::vector<EvalPoint> &evalPointList,
                        int maxEvalPoints = 0) const override;

    /// \brief Find using custom criteria.
    /**
     All the points for which crit() return true are put in evalPointList.

     \param crit                        The criteria function                               -- \b IN.
     \param evalPointList    The eval points within verifying the criteria -- \b OUT.
     \return                 The number of eval points found.
     */
    virtual size_t find(std::function<bool(const EvalPoint&)> crit,
                        std::vector<EvalPoint> &evalPointList) const override;
    
    /// Browse cache using criteria. The function can have access to remote info using the lambda
    /// function capture by reference.
    /**
    \param crit            The criteria function                               -- \b IN.
    */
    virtual void browse(std::function<void(const EvalPoint&)> crit) const override;


    /// \brief Find using custom criteria  and distance to a point.
    /**
     All the points for which the two crit() functions return true are put in evalPointList.

     \param crit1                      The first criteria function                                    -- \b IN.
     \param crit2                      The second criteria function                               -- \b IN.
     \param evalPointList    The eval points verifying the criteria -- \b OUT.
     \return                 The number of eval points found.
     */
    virtual size_t find(std::function<bool(const EvalPoint&)> crit1,
                        std::function<bool(const EvalPoint&)> crit2,
                        std::vector<EvalPoint> &evalPointList) const override;


    /// Get all non dominated (or equal) best feasible eval points using dominance criterion
    /// Used for multiobjective optimization
    /// NB: To use with precaution, computationally costly (n log n for two objectives).
    /**
     \param evalPointList   The best non dominated feasible eval points in a list  -- \b OUT.
     \param fixedVariable   Searching for a subproblem defined by this point -- \b IN.
     \param computeType   Which type of computation (eval type, compute type and h norm type)  -- \b IN.
     \return                The number of eval points found.
     */
    virtual size_t findBestFeas(std::vector<EvalPoint> &evalPointList,
                                const Point& fixedVariable = Point(),
                                const FHComputeType & computeType = defaultFHComputeType) const override;


    /// Find best infeasible points with h<=hmax:
    ///  -> index 0 and above if duplicates, least infeasible point with smallest f
    ///  -> last index and below if duplicates, best f with smallest h
    /// All best f points have the same blackbox outputs. Idem for the least infeasible points.
    /// Works also for multiobjective optimization
    /**
     \param evalPointList   The best non dominated feasible eval points in a list  -- \b OUT.
     \param fixedVariable   Searching for a subproblem defined by this point -- \b IN.
     \param hMax            Select a point if h <= hMax                                                 -- \b IN.
     \param computeType   Which type of computation (eval type, compute type and h norm type)  -- \b IN.
     \return                The number of eval points found.
     */
    virtual size_t findBestInf(std::vector<EvalPoint> &evalPointList,
                               const Double& hMax = INF,
                               const Point& fixedVariable = Point(),
                               const FHComputeType & computeType = defaultFHComputeType) const override;


    /// Get all non dominated (or equal) infeasible eval points using dominance criterion
    /// Works for single and multiobjective optimization
    /// NB: To use with precaution, computationally costly (O(n^2 m) where n is the number
    /// of points in the cache and m the number of objectives)
    /**
     \param evalPointList   The non dominated infeasible eval points  -- \b OUT.
     \param fixedVariable   Searching for a subproblem defined by this point -- \b IN.
     \param hMax            Select a point if h <= hMax                                                 -- \b IN.
     \param computeType   Which type of computation (eval type, compute type and h norm type)  -- \b IN.
     \return                The number of eval points found.
     */
    virtual size_t findFilterInf(std::vector<NOMAD::EvalPoint> &evalPointList,
                                 const Double& hMax,
                                 const Point& fixedVariable,
                                 const FHComputeType & computeType) const override;


    /// \brief Update EvalPoint in cache.
    /**
     * Look for Point and update the Eval part.
     * Eval is assumed non-NULL.
     * If the point is not found, throw an exception.

     \param evalPoint       The eval point to update  -- \b IN.
     \param evalType         Which eval of the EvalPoint to look at -- \b IN.
     \param mesh                  Update the eval point with a mesh -- \b IN.
     \return            A boolean indicating if update succeeded (\c true), \c false if there was an error.
     */
    bool update(const EvalPoint& evalPoint, EvalType  evalType, const MeshBasePtr mesh) override;

    /// Return number of eval points in the cache.
    size_t size() const override { return _cache.size(); }

    /// Empty the cache.
    bool clear() override;

    /// Clear all model (sgtelib) evaluations from the cache
    void clearModelEval(const int mainThreadNum) override;

    /** Purge the cache to get under CACHE_SIZE_MAX.
     */
    void purge() override;

    /// Write cache to file _filename.
    bool write() const override;

    /// Read file given by _filename.
    bool read() override;

    /// Display all points in cache.
    std::string displayAll() const override;

    /** Display only EvalPoints that have an eval.
     * This method is used to write the cache file.
     * \note the EvalPoint's Eval must satisfy method Eval::goodForCacheFile().
     */
    std::ostream& displayPointsWithEval(std::ostream& os) const;

    /// Compute the mean f.
    /**
     \param mean       The eval point to update -- \b OUT.
     \return        The number of EvalPoints for which f is defined.
     */
    size_t computeMeanF(Double &mean) const override;

    /// Call function func() on all EvalPoint in cache for Evals that were generated by mainThreadNum.
    void processOnAllPoints(void (*func)(EvalPoint&), const int mainThreadNum = -1) override;

    void deleteModelEvalOnly(const int mainThreadNum) override;
    
    /// Move eval points from cache set to cache set for rerun
    void moveEvalPointToCacheForRerun() override;

private:
    /// Private initialization function for internal use by constructor.
    void init();

    /// Private function for internal use by destructor.
    void destroy();

    /// Helper function for find and insertion.
    /**
     Throw exception if error. Do nothing otherwise.

     \param point       The point to verify  -- \b IN.
     */
    void verifyPointComplete(const Point& point) const;

    /// Helper function for find and insertion.
    /**
     Throw exception if error. Do nothing otherwise.

     \param point       The point to verify  -- \b IN.
     */
    void verifyPointSize(const Point& point) const;

    /// Helper function for find and insertion.
    /**
     Throw exception if error. Do nothing otherwise.

     \param evalPoint       The point to verify  -- \b IN.
     */
    void verifyPointComplete(const EvalPoint& evalPoint) const;

    /// Helper function for find and insertion.
    /**
     Throw exception if error. Do nothing otherwise.

     \param evalPoint       The point to verify  -- \b IN.
     */
    void verifyPointSize(const EvalPoint& evalPoint) const;
};

/// Display only EvalPoints that have an eval.
DLL_EVAL_API std::ostream& operator<<(std::ostream& os, const CacheSet& cache);

/// Get these EvalPoints from stream
DLL_EVAL_API std::istream& operator>>(std::istream& is, CacheSet& cache);

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_5_CACHESET__
