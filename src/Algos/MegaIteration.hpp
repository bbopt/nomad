#ifndef __NOMAD400_MEGAITERATION__
#define __NOMAD400_MEGAITERATION__

#include "../Algos/Iteration.hpp"
#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

/// Class to manage the iterations.
/**
 A mega iteration of an algorithm:
 *  - generates a lot of points over multiple strategies (ex., Poll and Search for Mads).
 *  - evaluates points
 *  - performs post-processing

 The run and start functions of mega iteration must be implemented for each optimizer
 that has several phases of points creation that we want to combine before launching evaluations (for example, ::MadsMegaIteration, ::NMMegaIteration, ...).

  \note As an hypothesis, the time load is taken by the evaluation,
  which is parallelized over all evaluations simultaneously.
  The iteration generation, including trial points generation,
  has little time load, so they do not need to be parallelized.
  It is also preferable to keep parallelization to the only place where
  it matters the most to avoid errors. \n
  There is no parallelization at the algorithmic level.
  The algorithms are run in main thread(s) only; other threads are available for evaluations.
*/
class MegaIteration: public Step
{


protected:
    std::vector<std::shared_ptr<Iteration>> _iterList; ///< A collection of additional iterations: Help generate more eval points.

    /**
     The barrier holds xFeas, xInf and hMax.
    xFeas and xInf will be used as frame centers.
    \note This barrier is in subspace.
     */
    std::shared_ptr<Barrier> _barrier;


    size_t _k;          ///< Main counter
    size_t _nbIterRun;  ///< Number of iterations run within this MegaIteration

    /**
     Success type of this step. Initialized with the run of the previous
     mega iteration, so that the update step knows what to do
     (for example, enlarge or reduce the mesh).
     At the end of run step during the mega iteration of an algorithm,
     the success type is updated (MegaIteration::setSuccessType) with the latest
     success type.
     */
    SuccessType _megaIterationSuccess;



public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling. Belongs to subproblem space. -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit MegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<Barrier> barrier,
                              SuccessType success);

    // No Destructor needed - keep defaults.


    /*---------*/
    /* Get/Set */
    /*---------*/

    size_t getNbIterations() const                              { return _iterList.size(); }
    std::shared_ptr<Iteration> getIter(size_t i) const          { return _iterList[i]; }

    size_t getK() const                                         { return _k; }
    void setK(const size_t k)                                   { _k = k; }
    size_t getNextK() const;

    // Barrier
    const std::shared_ptr<Barrier> getBarrier() const           { return _barrier; }
    void setBarrier(const std::shared_ptr<Barrier> &barrier)    { _barrier = barrier; }


    /**
     Success is initialized with the success of the previous
     MegaIteration. This is useful to know for the Update step, and some Search steps.
     At the end of the MegaIteration, Success is updated with the latest SuccessType.
     */
    const SuccessType& getSuccessType() const       { return _megaIterationSuccess; }
    void setSuccessType(const SuccessType& success) { _megaIterationSuccess = success; }

    /// Compute the number of xFeas and xInf points to use to create iterations
    void computeMaxXFeasXInf(size_t &maxXFeas, size_t &maxXInf);

    /*---------*/
    /* Others  */
    /*---------*/

    virtual void read(std::istream& is);
    virtual void display(std::ostream& os) const ;

private:
    /// Helper for constructor
    void init();

protected:

    virtual void startImp()    override = 0;
    virtual bool runImp()      override = 0;

    /// Implementation of the end task.
    /**
     Perform callback check for user termination and clear the iteration list.
     */
    virtual void endImp()      override;

};


/**
 Display useful values so that a new MegaIteration could be constructed using these values.
 */
std::ostream& operator<<(std::ostream& os, const MegaIteration& megaIteration);

/// Get an MegaIteration values from a stream
std::istream& operator>>(std::istream& is, MegaIteration& megaIteration);



#include "../nomad_nsend.hpp"

#endif // __NOMAD400_MEGAITERATION__
