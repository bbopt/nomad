#ifndef __NOMAD400_ORTHO2NPOLLMETHOD__
#define __NOMAD400_ORTHO2NPOLLMETHOD__

#include "../../Algos/Mads/PollMethodBase.hpp"
#include "../../nomad_nsbegin.hpp"

/// Class to perform Orth 2N Poll .
class Ortho2NPollMethod final : public PollMethodBase
{
public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit Ortho2NPollMethod(const Step* parentStep)
      : PollMethodBase(parentStep)
    {
        init();
    }

private:

    /// Helper for constructor.
    /**
     Test if the Ortho2N Poll is enabled or not.
     */
    void init();

    ///Generate 2n polls direction on a unit N-Sphere (no evaluation)
    /**
     - A single direction on unit n-sphere is computed (Poll::computeDirOnUnitSphere).
     - This direction is transformed into 2n directions on a unit n-sphere using the householder transformation.
     \param directions  The directions obtained for this poll -- \b OUT.
     \param n                      The dimension of the variable space  -- \b IN.
      */
     void generateUnitPollDirections(std::list<Direction> &directions, size_t n) const override  ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_ORTHO2NPOLLMETHOD__
