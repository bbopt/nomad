#ifndef __NOMAD400_NMINITIALIZATION__
#define __NOMAD400_NMINITIALIZATION__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Initialization.hpp"
#include "../../Algos/NelderMead/NMIterationUtils.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for Nelder Mead initialization
/**
 * For Step 0 of NM algorithm: Create and evaluate (if needed) trial points that will be used to form simplex. The points are put in cache.
 */
class NMInitialization: public Initialization, public NMIterationUtils
{
private:

    std::shared_ptr<AlgoStopReasons<NMStopType>> _nmStopReason;

public:
    /// Constructor
    /*
     \param parentStep      The parent of this step -- \b IN.
     */
    explicit NMInitialization(const Step* parentStep)
      : Initialization(parentStep),
        NMIterationUtils(parentStep)
    {
        init();
    }

    /// Destructor
    virtual ~NMInitialization() {}


private:
    /// Helper for constructor
    void init();

    /// Implementation of start task
    /**
     If needed, generate trial points and put them in cache to form simplex.
     For a standalone optimization (NM_OPTIMIZATION true), initial trial points must be generated to form a valid simplex around x0. Otherwise, the cache will be used to construct the simplex.
     */
    virtual void startImp() override ;

    /// Implementation of run task
    /**
     For a standalone NM, evaluate the trial points generated during start (simplex is created later)
     Otherwise, there are no trial points available and a failed stop reason is set.
     */
    virtual bool runImp() override ;

    // Update _evalPointList member with evaluated trial points for future use
    void endImp() override;

    /// Generate new points to form simplex
    void generateTrialPoints() override;

    /// Helper for start
    bool checkCacheCanFormSimplex ( void );

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_NMINITIALIZATION__
