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

#ifndef __NOMAD_4_0_CACHEINTERFACE__
#define __NOMAD_4_0_CACHEINTERFACE__

#include "Step.hpp"

#include "../nomad_nsbegin.hpp"


// Class interface for the cache.
/**
 * Used by algorithm and step classes.
 * The CacheInterface takes care of converting points from subproblems
   to full dimension before adding them to the cache, and from
   full dimension to subproblems when retrieving them from the cache.
 */
class CacheInterface
{
private:

    const Step* _step;      ///< Step that uses the Cache
    Point _fixedVariable;   ///< Full dimension point including fixed variables

public:
    /// Constructor
    /**
     \param step            The step using this CacheInterface
     */
    explicit CacheInterface(const Step* step)
      : _step(step)
    {
        init();
    }

    /// Find best feasible point(s) in cache
    /**
     \param evalPointList  The found evaluation points -- \b OUT.
     \param evalType             Criteria for EvalType -- \b IN.
     \param refeval               The point of reference                                      -- \b IN.
     \return                Number of points found
     */
    size_t findBestFeas(std::vector<EvalPoint> &evalPointList,
                        const EvalType& evalType,
                        const ComputeType& computeType,
                        const Eval* refeval) const;

    /// Find best infeasible point(s) in cache
    /**
     \param evalPointList   The found evaluation points -- \b OUT.
     \param hMax                       Points' h value must be under this value -- \b IN.
     \param evalType              Points' EvalType to look at -- \b IN.
     \param refeval                The point of reference                   -- \b IN.
     \return                 Number of points found
     */
    size_t findBestInf(std::vector<EvalPoint> &evalPointList,
                       const Double& hMax,
                       const EvalType& evalType,
                       const ComputeType& computeType,
                       const Eval* refeval) const;

    /// Interface for CacheBase::smartInsert.
    /**
     The full dimension point is reconstructed from step fixed variables information.
     \param evalPoint     The point to insert -- \b IN.
     \param maxNumberEval Maximun number of times this point may be evaluated -- \b IN.
     \param evalType      Criteria for EvalType -- \b IN.
     \return              \c True if the point may be sent for evaluation, \c false otherwise
     */
    bool smartInsert(const EvalPoint &evalPoint,
                     const short maxNumberEval = 1,
                     const EvalType& evalType = EvalType::BB);

    /// Interface for CacheBase::smartFind.
    /**
     Transform the point into full space before looking into cache.
     \param x           The point to find -- \b IN.
     \param evalPoint   The found EvalPoint -- \b OUT.
     \param evalType    If not UNDEFINED, wait for Eval of this EvalType to be completed. -- \b IN.
     */
    size_t find(const Point& x, EvalPoint &evalPoint,
                const EvalType& evalType = EvalType::UNDEFINED);

    /// Find points in the cache fulfilling a criteria
    /**
     \param crit                        The criteria function (function of EvalPoint) -- \b IN.
     \param evalPointList    The vector of EvalPoints found -- \b OUT.
     \param findInSubspace  The flag to find in subspace -- \b IN.
     \return              The number of points found
    */
    size_t find(std::function<bool(const EvalPoint&)> crit,
                std::vector<EvalPoint> &evalPointList,
                bool findInSubspace = false ) const;


    /// Get all points from the cache
    /**
     \param evalPointList The vector of EvalPoints -- \b OUT
     \return              The number of points
    */
    size_t getAllPoints(std::vector<EvalPoint> &evalPointList) const;

private:

    /// Helper for constructor
    void init();

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_CACHEINTERFACE__
