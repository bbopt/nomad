#ifndef __NOMAD400_SGTELIB_MODEL_MEGAITERATION__
#define __NOMAD400_SGTELIB_MODEL_MEGAITERATION__

// Manager for SgtelibModel iterations.
// Steps:
// - Generate points, using sgte model
// - Evaluate points
// - Post-processing

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/MegaIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// class SgtelibModelMegaIteration (Step)
class SgtelibModelMegaIteration: public MegaIteration, public IterationUtils
{
public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit SgtelibModelMegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<Barrier> barrier,
                              SuccessType success)
      : MegaIteration(parentStep, k, barrier, success),
        IterationUtils(parentStep)
    {
        init();
    }

    virtual ~SgtelibModelMegaIteration();

    /// Generate new points to evaluate
    void generateTrialPoints() override;


private:
    void init();

    virtual void startImp() override ;
    virtual bool runImp() override;
    virtual void endImp() override;

    // Helper for generateTrialPoints()
    void generateIterations();
    void runIterationsAndSetTrialPoints();
    void filterCache();

};

/** Display useful values so that a new MegaIteration could be constructed using these values.
 */
std::ostream& operator<<(std::ostream& os, const SgtelibModelMegaIteration& megaIteration);

/// Get an MegaIteration values from a stream
std::istream& operator>>(std::istream& is, SgtelibModelMegaIteration& megaIteration);

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIB_MODEL_MEGAITERATION__
