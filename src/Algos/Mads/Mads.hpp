
#ifndef __NOMAD400_MADS__
#define __NOMAD400_MADS__

#include "../../Algos/Algorithm.hpp"
#include "../../Algos/AlgoStopReasons.hpp"

#include "../../nomad_nsbegin.hpp"


/// The (M)esh (A)daptive (D)irect (S)earch algorithm.
/**
\note AllParameters and EvaluatorControl are held by MainStep.
Cache is a singleton all by itself.
MegaIteration holds the algorithm-related structures: Mesh, Barrier.
 */
class Mads: public Algorithm
{
public:
    /// Constructor
    /**
     \param parentStep          The parent of this step -- \b IN.
     \param stopReasons         The stop reasons for MADS -- \b IN.
     \param runParams           The run parameters that control MADS -- \b IN.
     \param pbParams            The problem parameters that control MADS -- \b IN.
     */
    explicit Mads(const Step* parentStep,
                  std::shared_ptr<AlgoStopReasons<MadsStopType>> stopReasons,
                  const std::shared_ptr<RunParameters>& runParams,
                  const std::shared_ptr<PbParameters>& pbParams)
      : Algorithm(parentStep, stopReasons, runParams, pbParams)
    {
        init();
    }


private:
    ///  Initialization of class, to be used by Constructor.
    void init();

    /// Algorithm execution for single-objective.
    /**
     Overrides the default algorithm's run
     \return \c true
     */
    virtual bool runImp() override;

    /// Helper for start()
    void readInformationForHotRestart() override ;



};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_MADS__
