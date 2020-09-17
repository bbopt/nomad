
#ifndef __NOMAD400_LHS__
#define __NOMAD400_LHS__

#include <vector>
#include "../Math/Point.hpp"
using namespace std;

#include "../nomad_nsbegin.hpp"


/// \brief Latin Hypercube Sampling class.
/**
 * Input:
 *  n dimension
 *  p number of desired samples
 *  lowerBound and upperBound of type ArrayOfDouble indicating lower and upper bounds. Must be completely defined.
 *  seed (optional)
 *
 * Output:
 *  p points of dimension n, distributed in the l x u hyper-rectangle of R^n
 */
class LHS
{
private:
    size_t _n;  ///< dimension
    size_t _p;  ///< number of samples
    ArrayOfDouble    _lowerBound; ///< lower bounds
    ArrayOfDouble    _upperBound; ///< upper bounds

public:
    /// Constructor
    /**
     \param n               Dimension -- \b IN.
     \param p               Number of samples -- \b IN.
     \param lowerBound      Lower bounds -- \b IN.
     \param upperBound      Upper bounds -- \b IN.
     */
    explicit LHS(const size_t n,
                 const size_t p,
                 const ArrayOfDouble lowerBound,
                 const ArrayOfDouble upperBound);

    /// Get lower bound
    /**
     \return Lower bound as \c ArrayOfDouble.
     */
    ArrayOfDouble getLowerBound() const                  { return _lowerBound; }

    /// Set lower bound
    /**
     \param lowerBound  An \c ArrayOfDouble for lower bound -- \b IN.
     */
    void setLowerBound(const ArrayOfDouble lowerBound)   { _lowerBound = lowerBound; }

    /// Get upper bound
    /**
     \return Upper bound as \c ArrayOfDouble.
     */
    ArrayOfDouble    getUpperBound(void) const           { return _upperBound; }

    /// Set upper bound
    /**
     \param upperBound  An \c ArrayOfDouble for upper bound -- \b IN.
     */
    void setUpperBound(const ArrayOfDouble upperBound)   { _upperBound = upperBound; }

    /// Do the sampling
    /**
     \return A vector \c Points.
     */
    std::vector<Point> Sample() const;

    /// Random permutation of the vector (1, 2, .., p)
    /**
     \param p   Number of positive integers elements in series -- \b IN.
     \return    Vector of positive integers
     */
    static std::vector<size_t> Permutation(const size_t p);
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_LHS__
