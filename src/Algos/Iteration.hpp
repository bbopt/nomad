#ifndef __NOMAD400_ITERATION__
#define __NOMAD400_ITERATION__

#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for iteration of an Algorithm.
/**
 This is an abstract class, each algorithm must implement its own iteration.
 */
class Iteration: public Step
{
protected:

    size_t _k; ///< Iteration number

    void init(); ///< Utility for constructor

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param k                  The iteration number -- \b IN.
     */
    explicit Iteration(const Step *parentStep,
                       const size_t k)
      : Step( parentStep ),
        _k(k)
    {
        init();
    }

    /// Destructor
    /**
     When iteration is done, Flush prints output queue.
     */
    virtual ~Iteration();

    // Get/Set

    /// Get iteration number
    /**
     Iteration number is incremented when calling the default Iteration::start().
     */
    size_t getK() const { return _k; }

    /// Increment iteration number by one
    /// To be used only when a single Iteration is used over and over, e.g. Nelder Mead
    void incK() { _k++; }

    /**
     \return \c nullptr for algorithms that do not use a mesh. Otherwise, this function must be reimplemented in algorithm specific iteration (for example, MadsIteration, NMIteration).
     */
    virtual const std::shared_ptr<MeshBase> getMesh() const { return nullptr; }

    /**
     \return \c nullptr for algorithms that do not use a frame center. Otherwise, this function must be reimplemented in algorithm specific iteration (for example, MadsIteration, NMIteration, ...).
     */
    virtual const std::shared_ptr<EvalPoint> getFrameCenter() const { return nullptr; }

protected:

    /**
     This must be implemented when an algorithm has its own iteration.
     */
    virtual void startImp()    override = 0;

    /**
     This must be implemented when an algorithm has its own iteration.
     */
    virtual bool runImp()      override = 0;

    /**
     The default implement for end function displays the stop reason and calls the customized end function if provided by the user. \n
     If an end implementation function specific to an algorithm is required, it is convenient to call this function for default task.
     */
    virtual void endImp()      override;

    /**
     \return \c true if this Iteration is considered the principal iteration of its parent MegaIteration.
    By default, return true.
    */
    virtual bool isMainIteration() const { return true; }

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_ITERATION__
