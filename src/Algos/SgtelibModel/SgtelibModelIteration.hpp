#ifndef __NOMAD400_SGTELIB_MODEL_ITERATION__
#define __NOMAD400_SGTELIB_MODEL_ITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelOptimize.hpp"

#include "../../nomad_nsbegin.hpp"

/// \class Iteration (Step)
class SgtelibModelIteration: public Iteration
{
private:
    // Optimizer for model on sgte function
    std::shared_ptr<SgtelibModelOptimize> _optimize;

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param k                  The iteration number -- \b IN.
     */
    explicit SgtelibModelIteration(const Step *parentStep,
                        const size_t k)
      : Iteration(parentStep, k),
        _optimize(nullptr)
    {
        init();
    }


    // Get/Set
    /// Return oracle points found by the Optimizer
    const EvalPointSet& getOraclePoints() const;


    virtual void startImp() override;
    virtual bool runImp() override;


private:
    void init();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIB_MODEL_ITERATION__
