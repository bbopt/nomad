
#ifndef __NOMAD400_CACHEINTERFACE__
#define __NOMAD400_CACHEINTERFACE__

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
     \param evalPointList     The found evaluation points -- \b OUT.
     \param evalType          Criteria for EvalType -- \b IN.
     \return                  Number of points found
     */
    size_t findBestFeas(std::vector<EvalPoint> &evalPointList,
                        const EvalType& evalType,
                        const Eval* refeval) const;

    /// Find best infeasible point(s) in cache
    /**
     \param evalPointList   The found evaluation points -- \b OUT.
     \param hMax            Points' h value must be under this value -- \b IN
     \param evalType        Points' EvalType to look at -- \b IN.
     \return                Number of points found
     */
    size_t findBestInf(std::vector<EvalPoint> &evalPointList,
                       const Double& hMax,
                       const EvalType& evalType,
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
     */
    size_t find(const Point x, EvalPoint &evalPoint);

    /// Find points in the cache fulfilling a criteria
    /**
     \param crit          The criteria function (function of EvalPoint) -- \b IN
     \param evalPointList The vector of EvalPoints found -- \b OUT
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

#endif // __NOMAD400_CACHEINTERFACE__
