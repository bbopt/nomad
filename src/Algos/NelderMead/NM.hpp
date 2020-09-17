#ifndef __NOMAD400_NM__
#define __NOMAD400_NM__


#include "../../Algos/Algorithm.hpp"
#include "../../Algos/AlgoStopReasons.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class implementing Nelder Mead non-linear simplex algorithm for constrained problems.
/**
 See the NM-Mads paper: https://link.springer.com/article/10.1007/s10589-018-0016-0 for details.
 */
class NM: public Algorithm
{
public:
    /// Constructor
    /**
     \param parentStep          The parent of this Step -- \b IN.
     \param stopReasons         The stop reasons for NM -- \b IN.
     \param runParams           The run parameters that control NM -- \b IN.
     \param pbParams            The problem parameters that control NM -- \b IN.
     */
    explicit NM(const Step* parentStep,
                std::shared_ptr<AlgoStopReasons<NMStopType>> stopReasons,
                const std::shared_ptr<RunParameters>& runParams,
                const std::shared_ptr<PbParameters>& pbParams)
      : Algorithm(parentStep, stopReasons, runParams, pbParams)
    {
        init();
    }

    /// Destructor
    virtual ~NM() {}

    virtual void readInformationForHotRestart() override ;

private:
    /// Helper for constructor
    void init();

    /// Implementation for run tasks.
    /**
     - Algorithm execution for single-objective.
     - Loop on NMMegaIteration (start, run, end) until a stop reason to terminate is obtained.
     - Update the succes type
     - Perform Termination tasks (start, run, end)
     - Update the SearchMethod success type with best success found.
     \return \c true
     */
    virtual bool runImp() override;

    /// Implementation for start tasks.
    /**
     - Set the stop reason to STARTED
     - Reset sub-algorithm counter
     - Perform Initialization tasks (start, run, end)
     */
    virtual void startImp() override;



};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_NM__
