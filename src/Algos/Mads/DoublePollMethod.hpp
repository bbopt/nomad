#ifndef __NOMAD400_DOUBLEPOLLMETHOD__
#define __NOMAD400_DOUBLEPOLLMETHOD__

#include "../../Algos/Mads/PollMethodBase.hpp"
#include "../../nomad_nsbegin.hpp"

/// Class to perform two poll directions (a single direction and its opposite).
class DoublePollMethod final : public PollMethodBase
{
public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit DoublePollMethod(const Step* parentStep )
      : PollMethodBase(parentStep)
    {
        init();
    }

private:

    /// Helper for constructor.
    void init();

    ///Generate 2 poll directions on a unit n-sphere
    /**
     - A single direction on unit n-sphere is computed (Poll::computeDirOnUnitSphere).
     - The negative direction is added.
     \param directions  The directions obtained for this poll -- \b OUT.
     \param n                      The dimension of the variable space -- \b IN.
      */
     void generateUnitPollDirections(std::list<Direction> &directions, size_t n) const override  ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_DOUBLEPOLLMETHOD__
