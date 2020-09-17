#ifndef __NOMAD400_INITIALIZATION__
#define __NOMAD400_INITIALIZATION__

#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for initialization (step 0) of an Algorithm
/**
 This an abstract class, each algorithm should probably implement an initialization.
 */
class Initialization: public Step
{
protected:
    std::shared_ptr<Barrier> _barrier;   ///< Barrier constructed from evaluated X0s

public:
    /// Constructor
    /*
     \param parentStep      The parent of this step -- \b IN.
     */
    explicit Initialization(const Step* parentStep)
      : Step(parentStep),
        _barrier(nullptr)
    {
        init();
    }

    /// Destructor
    /**
     Upon destruction, print all that is in the output queue.
     */
    virtual ~Initialization();

    std::shared_ptr<Barrier> getBarrier() const { return _barrier; }

private:
    /// Helper for constructor
    void init();

public:

    virtual void startImp()    override {}
    virtual bool runImp()      override = 0;
    virtual void endImp()      override {}

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_INITIALIZATION__
