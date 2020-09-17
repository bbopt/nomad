#ifndef __NOMAD400_QUAD_MODEL_INITIALIZATION__
#define __NOMAD400_QUAD_MODEL_INITIALIZATION__

#include "../../Algos/Initialization.hpp"
#include "../../Algos/IterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"


class QuadModelInitialization: public Initialization, public IterationUtils
{
private:
    std::shared_ptr<AlgoStopReasons<ModelStopType>> _qmStopReason;
public:
    /// Constructor
    explicit QuadModelInitialization(const Step* parentStep)
      : Initialization(parentStep),
        IterationUtils(parentStep)
    {
        init();
    }

    virtual ~QuadModelInitialization();


private:
    void init();

    virtual void startImp() override;
    virtual bool runImp() override;
    void endImp() override {};

    /// Insert X0s for evaluation or (exclusive) check cache
    void generateTrialPoints() override;

    /// Eval X0s, using blackbox.
    bool eval_x0s();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_INITIALIZATION__

