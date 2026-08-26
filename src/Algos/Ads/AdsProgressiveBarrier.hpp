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

#ifndef __NOMAD_4_6_ADSPROGRESSIVEBARRIER__
#define __NOMAD_4_6_ADSPROGRESSIVEBARRIER__

#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Eval/EvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for ADS progressive barrier.
/**
 * This barrier extends ProgressiveBarrier with functionality specific to ADS algorithm.
 * It maintains a list of cache point tags that can be used for punctured space checks.
 */
class DLL_ALGO_API AdsProgressiveBarrier : public ProgressiveBarrier
{
private:
    /// List of cache point tags to use for punctured space checks
    std::vector<size_t> _cachePointTags;

public:
    /// Constructor
    /**
     * hMax will be updated during optimization.
     \param hMax            The max of h to keep a point in the barrier -- \b IN.
     \param fixedVariable   The fixed variables have a fixed value -- \b IN.
     \param evalType        Type of evaluation (BB or MODEL) -- \b IN.
     \param computeType  Type of function computation (standard, phase-one or user) -- \b IN.
     \param evalPointList   Additional points to consider in building the barrier -- \b IN.
     \param barrierInitializedFromCache Flag to initialize the barrier from cache -- \b IN.
     */
    AdsProgressiveBarrier(const Double& hMax = INF,
            const Point& fixedVariable = Point(),
            EvalType evalType = EvalType::BB,
            FHComputeTypeS computeType = defaultFHComputeTypeS,
            const std::vector<EvalPoint>& evalPointList = std::vector<EvalPoint>(),
            bool barrierInitializedFromCache = true)
      : ProgressiveBarrier(hMax, fixedVariable, evalType, computeType, evalPointList, barrierInitializedFromCache)
    {
    }
    
    /// Copy constructor
    AdsProgressiveBarrier(const AdsProgressiveBarrier & b) : ProgressiveBarrier(b)
    {
        _cachePointTags = b._cachePointTags;
    }

    std::shared_ptr<BarrierBase> clone() const override
    {
        return std::make_shared<AdsProgressiveBarrier>(*this);
    }

    /*-----------------*/
    /* Cache point tags management */
    /*-----------------*/

    /// Add a cache point tag to the list
    /**
     \param tag    The tag to add -- \b IN.
     */
    void addCachePointTag(size_t tag);

    /// Add multiple cache point tags to the list
    /**
     \param tags   The tags to add -- \b IN.
     */
    void addCachePointTags(const std::vector<size_t>& tags);

    /// Remove a cache point tag from the list
    /**
     \param tag    The tag to remove -- \b IN.
     */
    void removeCachePointTag(size_t tag);

    /// Clear all cache point tags
    void clearCachePointTags();

    /// Get all cache point tags
    /**
     \return   The list of cache point tags.
     */
    const std::vector<size_t>& getCachePointTags() const { return _cachePointTags; }

    /// Check if a trial point is in the punctured space using the stored cache point tags
    /**
     * This method checks if a trial point is in the punctured space by comparing
     * it only with the cache points whose tags are in _cachePointTags.
     * The points are sorted according to the barrier ordering (best points first).
     *
     \param trialPoint      The trial point to check -- \b IN.
     \param minDelta        The minimum mesh size (delta) -- \b IN.
     \param iterationTag    The tag of the current iteration -- \b IN.
     \return                True if the point is in the punctured space (at most one point found), false otherwise.
     */
    bool isInPuncturedSpace(const Point& trialPoint,
                           const Double& minDelta,
                           size_t iterationTag) const;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_6_ADSPROGRESSIVEBARRIER__

