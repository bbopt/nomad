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

#ifndef __NOMAD400_CACHEINTERFACE__
#define __NOMAD400_CACHEINTERFACE__

#include "../Cache/CacheBase.hpp"
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

    const Step* _step; ///< Step that uses the Cache

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
     \param evalPointList     The found evaluation points -- \b OUT.
     \param fixedVariable    Fixed variables have a fixed value  -- \b IN.
     \param evalType          Criteria for EvalType -- \b IN.
     \return                  Number of points found
     */
    size_t findBestFeas(std::vector<EvalPoint> &evalPointList,
                        const Point& fixedVariable,
                        const EvalType& evalType) const;

    /// Find best infeasible point(s) in cache
    /**
     \param evalPointList   The found evaluation points -- \b OUT.
     \param hMax            Points' h value must be under this value -- \b IN
     \param fixedVariable  Variables whose values must be set -- \b IN
     \param evalType        Points' EvalType to look at -- \b IN.
     \return                Number of points found
     */
    size_t findBestInf(std::vector<EvalPoint> &evalPointList,
                       const Double& hMax,
                       const Point& fixedVariable,
                       const EvalType& evalType) const;

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
     */
    size_t find(const Point x, EvalPoint &evalPoint);

    /// Find points in the cache fulfilling a criteria
    /**
     \param crit          The criteria function (function of EvalPoint) -- \b IN
     \param evalPointList The vector of EvalPoints found -- \b OUT
     \return              The number of points found
    */
    size_t find(bool (*crit)(const EvalPoint&),
                std::vector<EvalPoint> &evalPointList) const;

    /// Get all points from the cache
    /**
     \param evalPointList The vector of EvalPoints -- \b OUT
     \return              The number of points
    */
    size_t getAllPoints(std::vector<EvalPoint> &evalPointList) const;

private:

    /// Helper for constructor
    void init();

    /// Helper for converting points to subspace.
    static void convertPointListToSub(std::vector<EvalPoint>& evalPointList,
                                      const Point& fixedVariable);

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_CACHEINTERFACE__
