#ifndef __NOMAD400_NMSHRINK__
#define __NOMAD400_NMSHRINK__

#include "../../Algos/NelderMead/NMIterationUtils.hpp"
#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for shrink step for NM algorithm.
/**
 The SHRINK step is executed when NM is called as a standalone algorithm and when all reflective steps (REFLECT, EXPAND, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION) fail to obtain improvement. \n
 The class manages the creation and evaluation of the shrunk simplex (start, run, end).
 In the start function, before shrinking the simplex, an update of the main barrier must be performed by calling NMUpdate start, run and end. The shrunk simplex is obtained in start task. \n
 Evaluation of the new simplex is performed in run.
 */
class NMShrink: public Step , public NMIterationUtils
{
private:
    Double _gamma;
public:
    /// Constructor
    /**
     \param parentStep The parent of this NM step
     */
    explicit NMShrink(const Step* parentStep )
      : Step( parentStep ) ,
        NMIterationUtils ( parentStep )
    {
        init();
    }
    virtual ~NMShrink() {}

    /// Implementation of the start tasks for simplex shrink.
    /**
     - update the barrier
     - call NMShrink::generateTrialPoints
     */
    virtual void    startImp() override ;

    /// Implementation of the run task for simplex shrink.
    /**
     Evaluate the trial points.
     */
    virtual bool    runImp() override ;

    /// Implementation of the end task for simplex shrink.
     /**
      Call default IterationUtils::postProcessing.
      */
    virtual void    endImp() override ;

    /// Generate new points to evaluate
    /**
     The new shrunk simplex is obtained with the formula y[k] = y0[k] + _gamma*(yi[k]-y0[k]). Where y0 is the frame center and yi are elements of the previous simplex. Gamma must be a parameter in ]0;1]
     */
    void generateTrialPoints() override;


private:

    /// Helper for constructor
    void init();

    /*---------------------------------*/
    /* Private methods used by Shrink */
    /*---------------------------------*/


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_NMSHRINK__
