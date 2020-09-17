#ifndef __NOMAD400_LHSEARCHMETHOD__
#define __NOMAD400_LHSEARCHMETHOD__

#include "../../Algos/Mads/SearchMethodSimple.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the search method using Latin Hypercube.
/**
 The trial points a produced using a LH sampling. The bounds for LH are the problem bounds if available or else are obtained from the frame center and the frame size of MADS.
 */
class LHSearchMethod final : public SearchMethodSimple
{
public:
    /// Constructor
    explicit LHSearchMethod(const Step* parentStep )
      : SearchMethodSimple(parentStep )
    {
        init();
    }

    /**
     \copydoc SearchMethodSimple::generateTrialPointsImp \n
     \note For the LH search method for Mads, the generation of points uses LHS with bounds determined from the frame size and frame center of the currrent iteration of Mads.
     */
    void generateTrialPointsImp() override;

private:
    /// Helper for constructor
    void init();


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_LHSEARCHMETHOD__

