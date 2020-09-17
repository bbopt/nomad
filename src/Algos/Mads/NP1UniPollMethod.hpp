#ifndef __NOMAD400_NP1UNIPOLLMETHOD__
#define __NOMAD400_NP1UNIPOLLMETHOD__

#include "../../Algos/Mads/PollMethodBase.hpp"
#include "../../nomad_nsbegin.hpp"

/// Class to perform N+1 Uniform Poll .
class NP1UniPollMethod final : public PollMethodBase
{
public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit NP1UniPollMethod(const Step* parentStep)
      : PollMethodBase(parentStep)
    {
        init();
    }

private:

    /// Helper for constructor.
    void init();

    ///Generate n+1 poll directions uniformly distributed on a unit n-sphere
    /**
     - 2n directions on a unit n-sphere are computed (same as Ortho2nPollMethod).
     - These directions are transformed into n+1 directions on a unit n-sphere.
         -
     \param directions  The directions obtained for this poll -- \b OUT.
     \param n                      The dimension of the variable space  -- \b IN.
      */
     void generateUnitPollDirections(std::list<Direction> &directions, size_t n) const override  ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_ORTHONP1UNIPOLLMETHOD__
