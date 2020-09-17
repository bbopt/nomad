#ifndef __NOMAD400_MADSUPDATE__
#define __NOMAD400_MADSUPDATE__

#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for Step 1. of MADS algorithm: parameter update.
/**
 The update is performed when calling the MadsUpdate::run function.
 */
class MadsUpdate: public Step
{
public:
    // Constructor
    explicit MadsUpdate(const Step* parentStep)
      : Step(parentStep)
    {
        init();
    }


private:

    /// Helper for constructor to check for valid ancestor.
    void init();

    /// No implementation is required for start.
    virtual void    startImp() override {}

    /// Implementation of the run tasks.
    /**
     Gets the best feasible point (xFeas) and best infeasible point (xInf)
     from the cache, and updates the MegaIteration's Barrier member with it.
     Compares new values of xFeas and xInf with previous ones
     - i.e., compute success or failure.
     Enlarges or shrinks the delta (mesh size) and Delta (frame size)
     accordingly.
     */
    virtual bool    runImp()   override;

    /// No implementation is required for start.
    virtual void    endImp()   override {}

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_MADSUPDATE__
