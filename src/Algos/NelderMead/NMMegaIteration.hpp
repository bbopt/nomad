#ifndef __NOMAD400_NMMEGAITERATION__
#define __NOMAD400_NMMEGAITERATION__

#include "../../Algos/MegaIteration.hpp"
#include "../../Algos/NelderMead/NMIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for Nelder Mead mega iteration.
/**
 * Manager for Nelder Mead iterations starts, runs and ends.
 * Steps:
 * - Generate points over simplices (start).
 * - Evaluate points (run)
 * - Post-processing (end)
 */
class NMMegaIteration: public MegaIteration
{
protected:

    std::shared_ptr<NMIteration> _nmIteration;

private:

    void init();

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit NMMegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<Barrier> barrier,
                              SuccessType success)
      : MegaIteration(parentStep, k,barrier,success), _nmIteration(nullptr)
    {
        init();
    }
    // No Destructor needed - keep defaults.


    void read(  std::istream& is ) override;
    void display(  std::ostream& os ) const override ;

private:

    /// Implementation of start task.
    /**
     Create a Nelder Mead iteration for a current incumbent: use xFeas or xInf from the barrier. \n
     \note Running the algorithm requires a single iteration object with several start, run, end for the various iterations of the algorithm. This allows to easily maintain a proper simplex in a single NMIteration object.
     */
    virtual void startImp() override ;

    /// Implementation of run task.
    /**
     The algorithm iterations are started, ran and ended sequentially until a stop reason to terminate is obtained. \n
     We have a success if either a better xFeas or
     a dominating or partial success for xInf was found.
     See Algorithm 12.2 from DFBO.
     */
    virtual bool runImp() override;

};

/**
 Display useful values so that a new MegaIteration could be constructed using these values.
 */
std::ostream& operator<<(std::ostream& os, const NMMegaIteration& megaIteration);

/// Get an MegaIteration values from a stream
std::istream& operator>>(std::istream& is, NMMegaIteration& megaIteration);

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_NMMEGAITERATION__
