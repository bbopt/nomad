#ifndef __NOMAD400_TERMINATION__
#define __NOMAD400_TERMINATION__

#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

///  Class for termination of an algorithm.
/**
 The terminate function checks for termination criterions such as MAX_ITERATIONS, MAX_TIME, STOP_IF_FEASIBLE and set the stop reason.
 */
class Termination: public Step
{
public:
    /// Constructor
    explicit Termination(const Step* parentStep,
                         const std::shared_ptr<RunParameters>& runParams = nullptr,
                         const std::shared_ptr<PbParameters>& pbParams = nullptr)
      : Step(parentStep, runParams, pbParams)
    {
        init();
    }

    /// Destructor
    virtual ~Termination() {}

    /**
     The terminate function is called when algorithm are performing iterations during a run. At each iteration, we test if a stop criterion is reached.
     */
    virtual bool terminate(size_t iteration);

    virtual void    startImp() override; ///< Will update the step name

    /// Implementation for run task of algorithm Termination.
    /**
     \return \c true is a stop reason requires termination of an algorithm, \c false otherwise.
     */
    virtual bool    runImp()   override;

    /// Implementation for end tasks of algorithm Termination.
    /**
     Upon completing an algorithm run, this end function is called to display termination info.
     */
    virtual void    endImp()   override;

private:

    /// Helper for constructor
    void init();

    /// Helper for end
    bool solHasFeas() const;

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_TERMINATION__
